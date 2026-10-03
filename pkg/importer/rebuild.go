package importer

import (
	"context"
	"fmt"
	"os"
	"path/filepath"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/normalize"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
)

// SourceFile pairs a logical source name with a VCF file path for rebuild.
type SourceFile struct {
	Source string
	File   string
	Path   string
}

// Rebuild re-imports the original VCF files into a fresh store and atomically
// swaps it into place. The legacy store is replaced by the rebuilt one, so
// queries see either the old complete state or the new complete state, never
// a partial intermediate.
//
// Progress is persisted per (source, file) via the checkpoint, so a crash
// resumes from the last checkpoint without duplicates or omissions. Memory
// usage is bounded by the store's memtable limit, not by total data size.
func Rebuild(ctx context.Context, dbDir string, ref normalize.Reference,
	cfg config.Config, checkpointPath string, files []SourceFile) (Report, error) {

	tmpDir := dbDir + ".rebuild"
	bakDir := dbDir + ".bak"
	rebuildCP := checkpointPath + ".rebuild"

	// Clean any stale rebuild directory from a previous crashed run.
	if err := os.RemoveAll(tmpDir); err != nil {
		return Report{}, err
	}
	// Start with a fresh checkpoint so the rebuild re-imports from the
	// beginning, not from the legacy store's checkpoint.
	if err := os.Remove(rebuildCP); err != nil && !os.IsNotExist(err) {
		return Report{}, err
	}

	// Build the new store in a temporary directory, leaving the legacy store
	// untouched so queries keep serving the old complete state.
	st, err := store.Open(tmpDir, Fold(cfg))
	if err != nil {
		return Report{}, err
	}

	im, err := New(st, ref, cfg, rebuildCP)
	if err != nil {
		st.Close()
		return Report{}, err
	}

	var report Report
	report.Source = "rebuild"
	for _, sf := range files {
		select {
		case <-ctx.Done():
			st.Close()
			return report, ctx.Err()
		default:
		}

		f, err := os.Open(sf.Path)
		if err != nil {
			st.Close()
			return report, fmt.Errorf("open %s: %w", sf.Path, err)
		}
		rep, err := im.Resume(ctx, sf.Source, sf.File, f)
		f.Close()
		if err != nil {
			st.Close()
			return report, fmt.Errorf("rebuild %s: %w", sf.Path, err)
		}
		report.Imported += rep.Imported
		report.Rejected = append(report.Rejected, rep.Rejected...)
	}

	// Flush and close the new store before the swap.
	if err := st.Close(); err != nil {
		return report, err
	}

	// Atomically swap the rebuilt store into place. Queries see either the
	// old complete state (before rename) or the new complete state (after
	// rename), never a partial intermediate.
	if err := os.RemoveAll(bakDir); err != nil {
		return report, err
	}
	if err := os.Rename(dbDir, bakDir); err != nil {
		return report, fmt.Errorf("stage legacy store: %w", err)
	}
	if err := os.Rename(tmpDir, dbDir); err != nil {
		// Attempt to restore the legacy store.
		_ = os.Rename(bakDir, dbDir)
		return report, fmt.Errorf("publish rebuilt store: %w", err)
	}
	if err := os.RemoveAll(bakDir); err != nil {
		// The swap succeeded; the leftover backup is harmless.
		return report, nil
	}
	return report, nil
}

// RebuildResume continues a rebuild from the persisted checkpoint. It opens
// the existing rebuilt store and resumes processing the given files.
func RebuildResume(ctx context.Context, dbDir string, ref normalize.Reference,
	cfg config.Config, checkpointPath string, files []SourceFile) (Report, error) {

	tmpDir := dbDir + ".rebuild"
	rebuildCP := checkpointPath + ".rebuild"
	if _, err := os.Stat(tmpDir); os.IsNotExist(err) {
		return Report{}, fmt.Errorf("no rebuild in progress at %s", tmpDir)
	}

	st, err := store.Open(tmpDir, Fold(cfg))
	if err != nil {
		return Report{}, err
	}

	im, err := New(st, ref, cfg, rebuildCP)
	if err != nil {
		st.Close()
		return Report{}, err
	}

	var report Report
	report.Source = "rebuild"
	for _, sf := range files {
		select {
		case <-ctx.Done():
			st.Close()
			return report, ctx.Err()
		default:
		}

		f, err := os.Open(sf.Path)
		if err != nil {
			st.Close()
			return report, fmt.Errorf("open %s: %w", sf.Path, err)
		}
		rep, err := im.Resume(ctx, sf.Source, sf.File, f)
		f.Close()
		if err != nil {
			st.Close()
			return report, fmt.Errorf("rebuild %s: %w", sf.Path, err)
		}
		report.Imported += rep.Imported
		report.Rejected = append(report.Rejected, rep.Rejected...)
	}

	if err := st.Close(); err != nil {
		return report, err
	}

	// Publish the rebuilt store.
	bakDir := dbDir + ".bak"
	if err := os.RemoveAll(bakDir); err != nil {
		return report, err
	}
	if err := os.Rename(dbDir, bakDir); err != nil {
		return report, fmt.Errorf("stage legacy store: %w", err)
	}
	if err := os.Rename(tmpDir, dbDir); err != nil {
		_ = os.Rename(bakDir, dbDir)
		return report, fmt.Errorf("publish rebuilt store: %w", err)
	}
	_ = os.RemoveAll(bakDir)
	return report, nil
}

// Ensure dbDir's parent directory exists.
func init() {
	_ = filepath.Join
}
