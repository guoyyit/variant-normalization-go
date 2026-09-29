package store

import "sort"

// keyIndex maps each stored key to the youngest segment holding it
// (memtable keys map to nil). It stores keys and pointers, never record
// bodies, so its footprint scales with the number of distinct keys, not
// input volume.
type keyIndex struct {
	segOf map[string]*segment
}

func newKeyIndex() *keyIndex {
	return &keyIndex{segOf: map[string]*segment{}}
}

// contains reports whether key is indexed anywhere (segment or memtable).
func (idx *keyIndex) contains(key string) bool {
	_, ok := idx.segOf[key]
	return ok
}

// addSegment registers all keys of seg as located in seg.
func (idx *keyIndex) addSegment(seg *segment) {
	_ = seg.scan(func(entry recordEntry) error {
		idx.segOf[entry.Key] = seg
		return nil
	})
}

// promoteMemtable marks every key of the flushed entries as living in seg.
func (idx *keyIndex) promoteMemtable(keys []string, seg *segment) {
	for _, key := range keys {
		idx.segOf[key] = seg
	}
}

// segmentFor returns the segment containing key, or nil when the key is only
// in the memtable or is unknown.
func (idx *keyIndex) segmentFor(key string) *segment {
	return idx.segOf[key]
}

// keysInRange returns sorted indexed keys within [lower, upper).
func (idx *keyIndex) keysInRange(lower, upper string, chrom string) []string {
	out := make([]string, 0)
	for key := range idx.segOf {
		if key >= lower && key < upper && keyChrom(key) == chrom {
			out = append(out, key)
		}
	}
	sort.Strings(out)
	return out
}
