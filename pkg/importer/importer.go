package importer

import (
	"context"
	"errors"
	"fmt"
	"io"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/merge"
	"github.com/guoyyit/variant-normalization-go/pkg/normalize"
	"github.com/guoyyit/variant-normalization-go/pkg/split"
	"github.com/guoyyit/variant-normalization-go/pkg/stats"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
	"github.com/guoyyit/variant-normalization-go/pkg/vcf"
)

// Fold is the deterministic read-modify-write function applied by the store
// under its write lock. Replays carrying provenance already present in the
// stored record return the stored state unchanged (idempotent resume);
// otherwise records are merged and annotations are recomputed.
func Fold(cfg config.Config) func(existing, incoming variant.Record) variant.Record {
	return func(existing, incoming variant.Record) variant.Record {
		if len(existing.Sources) == 0 {
			stats.Recompute(&incoming, cfg)
			incoming.ConfigID = cfg.ID
			return incoming
		}
		if len(incoming.Sources) > 0 && containsSource(existing.Sources, incoming.Sources[0]) &&
			hasAlt(existing, incoming) {
			return existing
		}
		merged := merge.Merge(existing, incoming, cfg)
		stats.Recompute(&merged, cfg)
		return merged
	}
}

func containsSource(sources []variant.SourceRef, want variant.SourceRef) bool {
	for _, src := range sources {
		if src == want {
			return true
		}
	}
	return false
}

// hasAlt reports whether every incoming alt is already present in the
// stored group, used together with provenance equality to distinguish a true
// replay from a distinct alt of the same multi-allelic source line.
func hasAlt(existing, incoming variant.Record) bool {
	for _, alt := range incoming.Alt {
		found := false
		for _, have := range existing.Alt {
			if have == alt {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return len(incoming.Alt) > 0
}

// Rejection describes one record that could not be imported.
type Rejection struct {
	Line   int64
	Reason string
}

// Report summarizes one import run.
type Report struct {
	Source   string
	Lines    int64
	Imported int
	Rejected []Rejection
}

// Importer streams VCF text into the store with checkpointed progress.
type Importer struct {
	store      *store.Store
	ref        normalize.Reference
	cfg        config.Config
	checkpoint *checkpoint
	syncBatch  int
}

// New creates an Importer backed by st.
func New(st *store.Store, ref normalize.Reference, cfg config.Config, checkpointPath string) (*Importer, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	cp, err := openCheckpoint(checkpointPath)
	if err != nil {
		return nil, err
	}
	return &Importer{store: st, ref: ref, cfg: cfg, checkpoint: cp, syncBatch: 1000}, nil
}

// WithSyncBatch sets how many imported records accumulate between durable
// checkpoints. Larger batches improve throughput; a crash replays at most one
// batch, with provenance dedup keeping replay idempotent.
func (im *Importer) WithSyncBatch(n int) *Importer {
	if n > 0 {
		im.syncBatch = n
	}
	return im
}

// ResetProgress clears the checkpoint cursors. Rebuild uses this when it
// starts a new empty shadow generation so stale cursors from a previously
// completed generation cannot cause it to skip records.
func (im *Importer) ResetProgress() {
	im.checkpoint.reset()
}

// Import streams every record from r. Invalid lines are reported via
// Rejected but never abort the stream. Replaying the same source is safe:
// the (source,file,line) provenance tuple makes each upsert idempotent.
func (im *Importer) Import(ctx context.Context, sourceName, fileName string, r io.Reader) (Report, error) {
	return im.run(ctx, sourceName, fileName, r, 0)
}

// Resume continues r from the checkpointed byte offset. r must support
// seeking (a file), so the caller can reopen the original input.
func (im *Importer) Resume(ctx context.Context, sourceName, fileName string, r io.ReadSeeker) (Report, error) {
	prog := im.checkpoint.progress(sourceName)
	if prog.Offset > 0 {
		if _, err := r.Seek(prog.Offset, io.SeekStart); err != nil {
			return Report{}, err
		}
	}
	return im.run(ctx, sourceName, fileName, r, prog.Line)
}

func (im *Importer) run(ctx context.Context, sourceName, fileName string, r io.Reader, startLine int64) (Report, error) {
	report := Report{Source: sourceName}

	parser := vcf.NewParser(r)
	parser.SetLine(startLine)
	batch := 0
	lastOffset := parser.ByteOffset()
	lastLine := startLine

	for {
		select {
		case <-ctx.Done():
			return report, ctx.Err()
		default:
		}

		rec, err := parser.Next()
		if errors.Is(err, io.EOF) {
			if batch > 0 {
				if err := im.checkpointBatch(sourceName, lastOffset, lastLine); err != nil {
					return report, err
				}
			}
			report.Lines = parser.Line()
			return report, nil
		}
		if err != nil {
			var perr *vcf.ParseError
			if errors.As(err, &perr) {
				report.Rejected = append(report.Rejected, Rejection{
					Line: perr.Line, Reason: perr.Reason})
				continue
			}
			return report, err
		}

		rec.Sources = []variant.SourceRef{{
			Source: sourceName, File: fileName, Line: parser.Line(),
		}}

		normalized, err := im.normalizeRecord(rec)
		if err != nil {
			report.Rejected = append(report.Rejected, Rejection{
				Line: parser.Line(), Reason: err.Error()})
			continue
		}

		if err := im.ingest(normalized); err != nil {
			return report, fmt.Errorf("line %d: %w", parser.Line(), err)
		}
		report.Imported += len(normalized)
		batch += len(normalized)
		lastOffset = parser.ByteOffset()
		lastLine = parser.Line()
		if batch >= im.syncBatch {
			if err := im.checkpointBatch(sourceName, lastOffset, lastLine); err != nil {
				return report, err
			}
			batch = 0
		}
	}
}

// checkpointBatch durably flushes the WAL, then persists the cursor. Ordering
// guarantees a crash replays at most this batch; source provenance makes the
// replay idempotent (no duplicates, no omissions).
func (im *Importer) checkpointBatch(source string, offset, line int64) error {
	if err := im.store.Sync(); err != nil {
		return err
	}
	return im.checkpoint.advance(source, offset, line)
}

// normalizeRecord splits the raw multi-allelic record and normalizes each
// single-alt piece independently.
func (im *Importer) normalizeRecord(rec variant.Record) ([]variant.Record, error) {
	parts := split.Split(rec)
	out := make([]variant.Record, 0, len(parts))
	for _, part := range parts {
		norm, err := normalize.Normalize(part, im.ref)
		if err != nil {
			return nil, err
		}
		out = append(out, norm)
	}
	return out, nil
}

func (im *Importer) ingest(recs []variant.Record) error {
	for _, rec := range recs {
		// The store folds incoming records atomically under its write lock
		// (provenance dedup -> merge -> recompute), so the importer only
		// upserts the normalized single-alt record.
		if err := im.store.Upsert(rec); err != nil {
			return err
		}
	}
	return nil
}
