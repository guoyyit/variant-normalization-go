// Package rebuild re-imports the original reported VCF files into a fresh
// shadow generation of the store and atomically publishes it. Because the new
// generation is produced by the same parser/normalizer/split/merge pipeline
// as a clean import, every rebuilt record is identical to re-importing the
// raw data from scratch; records normalized by a previous, incorrect
// normalizer disappear when the generation pointer switches.
package rebuild

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/importer"
	"github.com/guoyyit/variant-normalization-go/pkg/normalize"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// Source is one original input file that contributes to a generation.
type Source struct {
	Name string
	Path string
}

// Result summarizes one rebuild run.
type Result struct {
	Generation int
}

// Option customizes a Rebuilder.
type Option func(*Rebuilder)

// WithSyncBatch bounds how many accepted records accumulate between durable
// checkpoints during the rebuild, keeping replay bounded on crash.
func WithSyncBatch(n int) Option {
	return func(b *Rebuilder) {
		if n > 0 {
			b.syncBatch = n
		}
	}
}

// shadowState records which shadow generation an in-progress rebuild is
// filling. It lets an interrupted run resume the same shadow (keeping its
// checkpoint cursors) while a fresh run starts a new shadow with zero
// cursors.
type shadowState struct {
	Generation int    `json:"generation"`
	Dir        string `json:"dir"`
}

// Rebuilder streams the raw sources into a shadow generation.
type Rebuilder struct {
	router    *store.Router
	ref       normalize.Reference
	cfg       config.Config
	fold      func(existing, incoming variant.Record) variant.Record
	syncBatch int
}

func statePath(checkpointPath string) string { return checkpointPath + ".rebuild" }

// New creates a Rebuilder.
func New(
	router *store.Router,
	ref normalize.Reference,
	cfg config.Config,
	fold func(existing, incoming variant.Record) variant.Record,
	opts ...Option,
) (*Rebuilder, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	b := &Rebuilder{router: router, ref: ref, cfg: cfg, fold: fold, syncBatch: 1000}
	for _, opt := range opts {
		opt(b)
	}
	return b, nil
}

// Rebuild imports the sources into a shadow generation and publishes it.
//
// Memory is bounded by the shadow store's memtable, progress is durable per
// source, and an interrupted run resumes into the same shadow generation from
// its checkpoint; provenance-based folding makes replay idempotent (no
// duplicates, no omissions). Throughout the run the active generation keeps
// serving reads, and the switch is one atomic manifest publish plus pointer
// swap, so readers see either the complete old or complete new generation.
func (b *Rebuilder) Rebuild(ctx context.Context, sources []Source, checkpointPath string) (result Result, err error) {
	if len(sources) == 0 {
		return Result{}, fmt.Errorf("rebuild requires at least one source file")
	}

	shadow, dir, resumed, err := b.openShadow(checkpointPath)
	if err != nil {
		return Result{}, err
	}

	committed := false
	defer func() {
		if !committed {
			// Keep the shadow and its state on interruption/failure so a
			// follow-up Rebuild can resume; only close our handle.
			_ = shadow.Close()
		}
	}()

	im, err := importer.New(shadow, b.ref, b.cfg, checkpointPath)
	if err != nil {
		return Result{}, err
	}
	im.WithSyncBatch(b.syncBatch)
	if !resumed {
		im.ResetProgress()
		rel, rerr := filepath.Rel(b.router.RootDir(), dir)
		if rerr != nil {
			return Result{}, rerr
		}
		if werr := writeShadowState(statePath(checkpointPath), shadowState{
			Generation: b.router.Generation() + 1, Dir: rel,
		}); werr != nil {
			return Result{}, werr
		}
	}

	for _, src := range sources {
		f, oerr := os.Open(src.Path)
		if oerr != nil {
			return Result{}, fmt.Errorf("rebuild source %q: %w", src.Path, oerr)
		}
		_, ierr := im.Resume(ctx, src.Name, filepath.Base(src.Path), f)
		cerr := f.Close()
		if ierr != nil {
			return Result{}, fmt.Errorf("rebuild source %q: %w", src.Path, ierr)
		}
		if cerr != nil {
			return Result{}, fmt.Errorf("rebuild source %q: %w", src.Path, cerr)
		}
	}

	if cerr := b.router.CommitRebuild(shadow, dir); cerr != nil {
		return Result{}, cerr
	}
	committed = true
	if rerr := os.Remove(statePath(checkpointPath)); rerr != nil && !os.IsNotExist(rerr) {
		return Result{}, rerr
	}
	return Result{Generation: b.router.Generation()}, nil
}

// openShadow resumes an unfinished shadow when its persisted state matches
// the next generation; otherwise it creates a fresh shadow, first removing
// any stale shadow directory from an abandoned attempt.
func (b *Rebuilder) openShadow(checkpointPath string) (*store.Store, string, bool, error) {
	pending, err := readShadowState(statePath(checkpointPath))
	if err != nil {
		return nil, "", false, err
	}
	if pending != nil {
		full := filepath.Join(b.router.RootDir(), pending.Dir)
		if pending.Generation == b.router.Generation()+1 && dirExists(full) {
			st, oerr := store.Open(full, b.fold)
			if oerr != nil {
				return nil, "", false, oerr
			}
			return st, full, true, nil
		}
	}

	// Abandoned or superseded shadow state: drop it and start fresh.
	if pending != nil {
		stale := filepath.Join(b.router.RootDir(), pending.Dir)
		_ = os.RemoveAll(stale)
		_ = os.Remove(statePath(checkpointPath))
	}
	st, dir, err := b.router.BeginRebuild(b.fold)
	if err != nil {
		return nil, "", false, err
	}
	return st, dir, false, nil
}

func readShadowState(path string) (*shadowState, error) {
	raw, err := os.ReadFile(path)
	if errors.Is(err, os.ErrNotExist) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	var st shadowState
	if err := json.Unmarshal(raw, &st); err != nil {
		return nil, fmt.Errorf("corrupt rebuild state %s: %w", path, err)
	}
	return &st, nil
}

func writeShadowState(path string, st shadowState) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(path), ".rebuild-*")
	if err != nil {
		return err
	}
	name := tmp.Name()
	encErr := json.NewEncoder(tmp).Encode(st)
	syncErr := tmp.Sync()
	closeErr := tmp.Close()
	if encErr != nil {
		os.Remove(name)
		return encErr
	}
	if syncErr != nil {
		os.Remove(name)
		return syncErr
	}
	if closeErr != nil {
		os.Remove(name)
		return closeErr
	}
	return os.Rename(name, path)
}

func dirExists(p string) bool {
	info, err := os.Stat(p)
	return err == nil && info.IsDir()
}
