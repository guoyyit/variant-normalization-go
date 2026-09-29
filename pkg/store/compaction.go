package store

import (
	"os"
	"sort"
	"sync/atomic"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

var segmentSeq atomic.Int64

// scheduleCompact signals the background compactor when too many segments
// accumulated.
func (s *Store) scheduleCompact() {
	select {
	case s.compactCh <- struct{}{}:
	default:
	}
}

// runCompactor merges segments in the background while the store keeps
// serving reads from the old segment set until the new one is published.
func (s *Store) runCompactor(mergeFn func(existing, incoming variant.Record) variant.Record) {
	for {
		select {
		case <-s.stopCh:
			return
		case <-s.compactCh:
		}
		s.mu.RLock()
		need := len(s.segments) >= flushFraction
		s.mu.RUnlock()
		if !need {
			continue
		}
		s.mu.Lock()
		_ = s.compactLocked()
		s.mu.Unlock()
	}
}

// compactLocked merges every existing segment plus the memtable into one new
// segment, then removes the old files. Same-key entries appear in segment age
// order so mergeFn folds them in arrival order.
func (s *Store) compactLocked() error {
	old := append([]*segment{}, s.segments...)

	// A key may have stale copies in several segments: newer copies fully
	// replace older ones (segments hold complete states, not deltas). Walk
	// oldest -> youngest overwriting, then apply the memtable last.
	values := map[string]variant.Record{}
	for _, seg := range old {
		err := seg.scan(func(entry recordEntry) error {
			values[entry.Key] = entry.Record
			return nil
		})
		if err != nil {
			return err
		}
	}

	for key, rec := range s.memtable {
		values[key] = rec
	}

	keys := make([]string, 0, len(values))
	for key := range values {
		keys = append(keys, key)
	}
	sort.Strings(keys)

	entries := make([]recordEntry, 0, len(keys))
	for _, key := range keys {
		// Normalize through the fold identity branch so compacted output is
		// equivalent to freshly inserted state.
		rec := s.mergeFn(variant.Record{}, values[key])
		entries = append(entries, recordEntry{Key: key, Record: rec})
	}

	seg, err := writeSegment(s.dir, s.segmentName(), entries)
	if err != nil {
		return err
	}

	// Publish atomically under the lock, then delete obsolete files.
	s.segments = []*segment{seg}
	s.index = newKeyIndex()
	s.index.promoteMemtable(keys, seg)
	s.memtable = map[string]variant.Record{}
	s.memCount = 0
	if err := s.wal.Truncate(0); err != nil {
		return err
	}
	if _, err := s.wal.Seek(0, 0); err != nil {
		return err
	}

	for _, oldSeg := range old {
		_ = os.Remove(oldSeg.path)
	}
	return syncDir(s.dir)
}
