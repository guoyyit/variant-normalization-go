package importer

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// sourceProgress is the durable resume cursor for one input source.
type sourceProgress struct {
	// Offset is the byte position in the input immediately after the last
	// fully processed record.
	Offset int64 `json:"offset"`
	// Line is the 1-based line number of that record.
	Line int64 `json:"line"`
	// File is the input file name this progress refers to.
	File string `json:"file"`
}

// checkpointFile is the on-disk progress manifest.
type checkpointFile struct {
	Sources map[string]sourceProgress `json:"sources"`
}

// checkpoint persists per-source byte offsets. Every Advance is written to a
// temp file, fsynced and atomically renamed, so a crash leaves either the old
// or the new complete manifest - never a torn one.
type checkpoint struct {
	path string
	data checkpointFile
}

func openCheckpoint(path string) (*checkpoint, error) {
	cp := &checkpoint{path: path, data: checkpointFile{Sources: map[string]sourceProgress{}}}
	raw, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return cp, nil
		}
		return nil, err
	}
	if err := json.Unmarshal(raw, &cp.data); err != nil {
		return nil, fmt.Errorf("corrupt checkpoint %s: %w", path, err)
	}
	if cp.data.Sources == nil {
		cp.data.Sources = map[string]sourceProgress{}
	}
	return cp, nil
}

// checkpointKey returns the composite key for a (source, file) pair.
func checkpointKey(source, file string) string {
	return source + "\x00" + file
}

func (c *checkpoint) progress(source, file string) sourceProgress {
	return c.data.Sources[checkpointKey(source, file)]
}

// advance moves the cursor and atomically persists the manifest.
func (c *checkpoint) advance(source, file string, offset, line int64) error {
	c.data.Sources[checkpointKey(source, file)] = sourceProgress{
		Offset: offset, Line: line, File: file,
	}
	return c.persist()
}

func (c *checkpoint) persist() error {
	if err := os.MkdirAll(filepath.Dir(c.path), 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(c.path), ".checkpoint-*")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	enc := json.NewEncoder(tmp)
	if err := enc.Encode(c.data); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return err
	}
	if err := tmp.Sync(); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmpName)
		return err
	}
	if err := os.Rename(tmpName, c.path); err != nil {
		return err
	}
	dir, err := os.Open(filepath.Dir(c.path))
	if err != nil {
		return err
	}
	defer dir.Close()
	return dir.Sync()
}
