package store

import (
	"bufio"
	"encoding/json"
	"errors"
	"io"
	"os"
	"path/filepath"
	"sort"
)

// segment is an immutable, sorted-on-disk list of merged records.
// Records are stored as length-prefixed JSON lines, ordered by
// (chrom, pos, ref, alt-slot), so queries can binary-scan by key.
type segment struct {
	path  string
	count int
	min   string
	max   string
}

const segmentPrefix = "seg-"

// writeSegment atomically writes records to a new segment file.
func writeSegment(dir, name string, records []recordEntry) (*segment, error) {
	sort.SliceStable(records, func(i, j int) bool {
		return records[i].Key < records[j].Key
	})
	tmp, err := os.CreateTemp(dir, name+".tmp-*")
	if err != nil {
		return nil, err
	}
	tmpName := tmp.Name()
	writer := bufio.NewWriterSize(tmp, 1<<16)
	enc := json.NewEncoder(writer)
	for _, entry := range records {
		if err := enc.Encode(entry); err != nil {
			tmp.Close()
			os.Remove(tmpName)
			return nil, err
		}
	}
	if err := writer.Flush(); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return nil, err
	}
	if err := tmp.Sync(); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return nil, err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmpName)
		return nil, err
	}
	final := filepath.Join(dir, name)
	if err := os.Rename(tmpName, final); err != nil {
		return nil, err
	}
	if err := syncDir(dir); err != nil {
		return nil, err
	}
	seg := &segment{path: final, count: len(records)}
	if len(records) > 0 {
		seg.min = records[0].Key
		seg.max = records[len(records)-1].Key
	}
	return seg, nil
}

// scan streams all entries of the segment in key order.
func (s *segment) scan(fn func(recordEntry) error) error {
	return s.scanRange("", "", fn)
}

// lookup finds a single key in the sorted segment.
func (s *segment) lookup(want string) (recordEntry, bool, error) {
	var hit recordEntry
	found := false
	err := s.scanRange(want, successor(want), func(entry recordEntry) error {
		if entry.Key == want {
			hit = entry
			found = true
			return errStopScan
		}
		return nil
	})
	if errors.Is(err, errStopScan) {
		err = nil
	}
	return hit, found, err
}

func successor(key string) string {
	b := []byte(key)
	return string(append(b, 0))
}

// scanRange streams entries in the half-open key range [lower, upper).
// Empty bounds are unbounded. Entries are key-sorted, so scanning stops as
// soon as upper is passed.
func (s *segment) scanRange(lower, upper string, fn func(recordEntry) error) error {
	f, err := os.Open(s.path)
	if err != nil {
		return err
	}
	defer f.Close()

	dec := json.NewDecoder(bufio.NewReaderSize(f, 1<<16))
	for {
		var entry recordEntry
		if err := dec.Decode(&entry); err != nil {
			if errors.Is(err, io.EOF) {
				return nil
			}
			return err
		}
		if lower != "" && entry.Key < lower {
			continue
		}
		if upper != "" && entry.Key >= upper {
			return nil
		}
		if err := fn(entry); err != nil {
			return err
		}
	}
}

func syncDir(dir string) error {
	d, err := os.Open(dir)
	if err != nil {
		return err
	}
	err = d.Sync()
	closeErr := d.Close()
	if err != nil {
		return err
	}
	return closeErr
}
