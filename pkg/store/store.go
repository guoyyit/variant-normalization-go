package store

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// recordEntry is one stored, fully merged record with its sortable key.
type recordEntry struct {
	Key    string         `json:"k"`
	Record variant.Record `json:"r"`
}

const (
	defaultMemLimit = 256
	flushFraction   = 4
	posWidth        = 20
	walName         = "wal.log"
	walSyncInterval = 40 * time.Millisecond
)

// Store is a durable interval-queryable record store. Memory usage is bounded
// by memLimit regardless of total data: the memtable is flushed to immutable
// on-disk segments and compacted in the background.
type Store struct {
	dir       string
	memLimit  int
	mu        sync.RWMutex
	memtable  map[string]variant.Record
	memCount  int
	segments  []*segment
	wal       *os.File
	walEnc    *json.Encoder
	flushing  bool
	closed    bool
	compactCh chan struct{}
	stopCh    chan struct{}
	doneCh    chan struct{}
	mergeFn   func(existing, incoming variant.Record) variant.Record
	index     *keyIndex
	dirty     bool
	waiters   []chan struct{}
}

// Option customizes a Store at Open time.
type Option func(*Store)

// WithMemLimit bounds the number of records held in the memtable.
func WithMemLimit(n int) Option {
	return func(s *Store) {
		if n > 0 {
			s.memLimit = n
		}
	}
}

// Open opens (creating if needed) a store in dir and replays its WAL.
func Open(dir string, mergeFn func(existing, incoming variant.Record) variant.Record, opts ...Option) (*Store, error) {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}
	s := &Store{
		dir:       dir,
		memLimit:  defaultMemLimit,
		memtable:  map[string]variant.Record{},
		compactCh: make(chan struct{}, 1),
		stopCh:    make(chan struct{}),
		doneCh:    make(chan struct{}),
		mergeFn:   mergeFn,
		index:     newKeyIndex(),
	}
	for _, opt := range opts {
		opt(s)
	}

	if err := s.loadSegments(); err != nil {
		return nil, err
	}
	if err := s.openWAL(); err != nil {
		return nil, err
	}
	if err := s.replayWAL(); err != nil {
		return nil, err
	}
	go func() {
		defer close(s.doneCh)
		s.runCompactor(mergeFn)
	}()
	go s.runWALSyncer()
	go s.runWALSyncer()
	return s, nil
}

func (s *Store) loadSegments() error {
	files, err := os.ReadDir(s.dir)
	if err != nil {
		return err
	}
	var segs []*segment
	for _, f := range files {
		if !strings.HasPrefix(f.Name(), segmentPrefix) {
			continue
		}
		path := filepath.Join(s.dir, f.Name())
		count, min, max, err := inspectSegment(path)
		if err != nil {
			return err
		}
		segs = append(segs, &segment{path: path, count: count, min: min, max: max})
	}
	sort.SliceStable(segs, func(i, j int) bool {
		return filepath.Base(segs[i].path) < filepath.Base(segs[j].path)
	})
	s.segments = segs
	s.index = newKeyIndex()
	for _, seg := range segs {
		s.index.addSegment(seg)
	}
	for _, seg := range segs {
		name := filepath.Base(seg.path)
		num := strings.TrimPrefix(name, segmentPrefix)
		if n, err := strconv.ParseInt(num, 10, 64); err == nil && n > segmentSeq.Load() {
			segmentSeq.Store(n)
		}
	}
	return nil
}

func inspectSegment(path string) (int, string, string, error) {
	seg := &segment{path: path}
	count := 0
	min, max := "", ""
	err := seg.scan(func(entry recordEntry) error {
		if count == 0 {
			min = entry.Key
		}
		max = entry.Key
		count++
		return nil
	})
	return count, min, max, err
}

func (s *Store) walPath() string { return filepath.Join(s.dir, walName) }

func (s *Store) openWAL() error {
	f, err := os.OpenFile(s.walPath(), os.O_CREATE|os.O_RDWR, 0o644)
	if err != nil {
		return err
	}
	s.wal = f
	s.walEnc = json.NewEncoder(f)
	return nil
}

// replayWAL folds records left in the WAL from a previous process into the
// memtable. Entries are already upsert-merged on disk, so replay is
// idempotent.
func (s *Store) replayWAL() error {
	info, err := s.wal.Stat()
	if err != nil {
		return err
	}
	if info.Size() == 0 {
		return nil
	}
	if _, err := s.wal.Seek(0, 0); err != nil {
		return err
	}
	dec := json.NewDecoder(s.wal)
	for {
		var entry recordEntry
		if err := dec.Decode(&entry); err != nil {
			if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
				break
			}
			return err
		}
		s.memtable[entry.Key] = entry.Record
		s.memCount++
	}
	if _, err := s.wal.Seek(0, 0); err != nil {
		return err
	}
	if err := s.wal.Truncate(0); err != nil {
		return err
	}
	if _, err := s.wal.Seek(0, 0); err != nil {
		return err
	}
	s.walEnc = json.NewEncoder(s.wal)
	return nil
}

// EncodeKey returns the lexicographic sort key for chrom/pos/ref.
func EncodeKey(chrom string, pos int, ref string) string {
	return chrom + "\x00" + fmt.Sprintf("%0*d", posWidth, pos) + "\x00" + ref
}

func keyChrom(key string) string {
	if i := strings.IndexByte(key, '\x00'); i >= 0 {
		return key[:i]
	}
	return key
}

// Upsert atomically folds rec into the stored group at its key.
func (s *Store) Upsert(rec variant.Record) error {
	s.mu.Lock()
	if s.closed {
		s.mu.Unlock()
		return ErrClosed
	}
	key := EncodeKey(rec.Chrom, rec.Pos, rec.Ref)
	folded, err := s.foldAt(key, rec)
	if err != nil {
		s.mu.Unlock()
		return err
	}
	rec = folded

	entry := recordEntry{Key: key, Record: rec}
	if err := s.walEnc.Encode(entry); err != nil {
		s.mu.Unlock()
		return err
	}
	s.markDirtyLocked()

	if _, existed := s.memtable[key]; !existed {
		s.memCount++
	}
	s.memtable[key] = rec

	if s.memCount >= s.memLimit && !s.flushing {
		if err := s.flushLocked(); err != nil {
			s.mu.Unlock()
			return err
		}
	}
	s.mu.Unlock()
	return nil
}

// Sync blocks until every accepted Upsert is durable on disk. The importer
// calls it before advancing its checkpoint, so a crash never loses a record
// whose progress was already persisted.
func (s *Store) Sync() error {
	s.mu.RLock()
	closed := s.closed
	s.mu.RUnlock()
	if closed {
		return ErrClosed
	}
	return s.drainWAL()
}

// markDirtyLocked marks un-synced WAL bytes.
func (s *Store) markDirtyLocked() {
	s.dirty = true
}

// runWALSyncer group-commits pending WAL writes on a fixed cadence, so a
// burst of imports pays one fsync per interval instead of one per record.
func (s *Store) runWALSyncer() {
	ticker := time.NewTicker(walSyncInterval)
	defer ticker.Stop()
	for {
		select {
		case <-s.stopCh:
			return
		case <-ticker.C:
			s.syncWAL()
		}
	}
}

func (s *Store) syncWAL() {
	s.mu.Lock()
	if !s.dirty {
		s.mu.Unlock()
		return
	}
	s.dirty = false
	waiters := s.waiters
	s.waiters = nil
	err := s.wal.Sync()
	s.mu.Unlock()

	for _, wait := range waiters {
		close(wait)
	}
	_ = err
}

// drainWAL forces pending bytes durable before flush/shutdown.
func (s *Store) drainWAL() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.drainWALLocked()
}

// drainWALLocked performs the fsync while assuming the write lock is held.
// The lock is briefly released during the fsync syscall after detaching
// waiters, so concurrent readers keep running.
func (s *Store) drainWALLocked() error {
	if !s.dirty {
		return nil
	}
	s.dirty = false
	waiters := s.waiters
	s.waiters = nil
	s.mu.Unlock()
	defer s.mu.Lock()

	if err := s.wal.Sync(); err != nil {
		for _, wait := range waiters {
			close(wait)
		}
		return err
	}
	for _, wait := range waiters {
		close(wait)
	}
	return nil
}

// foldAt applies the configured fold against the stored state, or against an
// empty record on first insert, so initial annotation recompute also runs
// atomically under the write lock.
func (s *Store) foldAt(key string, incoming variant.Record) (variant.Record, error) {
	if existing, ok := s.memtable[key]; ok {
		return s.mergeFn(existing, incoming), nil
	}
	loaded, found, err := s.readLocked(key)
	if err != nil {
		return variant.Record{}, err
	}
	if found {
		return s.mergeFn(loaded, incoming), nil
	}
	return s.mergeFn(variant.Record{}, incoming), nil
}

// readLocked fetches the newest copy of key from memtable or newest-to-oldest
// segments via the key index. Invoked with the write lock held.
func (s *Store) readLocked(key string) (variant.Record, bool, error) {
	seg := s.index.segmentFor(key)
	if seg == nil {
		return variant.Record{}, false, nil
	}
	entry, found, err := seg.lookup(key)
	if err != nil || !found {
		return variant.Record{}, false, err
	}
	return entry.Record, true, nil
}

var errStopScan = fmt.Errorf("stop scan")

// flushLocked moves the memtable into a new immutable segment. The WAL is
// rotated only after the segment is durably on disk and fsynced.
func (s *Store) flushLocked() error {
	if s.memCount == 0 {
		return nil
	}
	if err := s.drainWALLocked(); err != nil {
		return err
	}
	entries := make([]recordEntry, 0, len(s.memtable))
	for key, rec := range s.memtable {
		entries = append(entries, recordEntry{Key: key, Record: rec})
	}

	// A flush snapshot keeps readers on the old state until publication.
	s.flushing = true
	seg, err := writeSegment(s.dir, s.segmentName(), entries)
	if err != nil {
		s.flushing = false
		return err
	}

	if err := s.wal.Truncate(0); err != nil {
		s.flushing = false
		return err
	}
	if _, err := s.wal.Seek(0, 0); err != nil {
		s.flushing = false
		return err
	}
	if err := s.wal.Sync(); err != nil {
		s.flushing = false
		return err
	}

	s.segments = append(s.segments, seg)
	flushedKeys := make([]string, 0, len(entries))
	for _, entry := range entries {
		flushedKeys = append(flushedKeys, entry.Key)
	}
	s.index.promoteMemtable(flushedKeys, seg)
	s.memtable = map[string]variant.Record{}
	s.memCount = 0
	s.flushing = false

	s.scheduleCompact()
	return nil
}

func (s *Store) segmentName() string {
	return fmt.Sprintf("%s%020d", segmentPrefix, segmentSeq.Add(1))
}

// Get returns the merged group at chrom/pos/ref.
func (s *Store) Get(chrom string, pos int, ref string) (variant.Record, bool, error) {
	key := EncodeKey(chrom, pos, ref)
	s.mu.RLock()
	defer s.mu.RUnlock()
	if rec, ok := s.memtable[key]; ok {
		return rec, true, nil
	}
	return s.readLocked(key)
}

// PageQuery selects records by chromosome and 1-based closed interval.
type PageQuery struct {
	Chrom     string
	Start     int
	End       int
	PageSize  int
	PageToken string
}

// Page is one stable query page plus continuation token.
type Page struct {
	Records   []variant.Record
	NextToken string
}

// Query returns a stable page of normalized records ordered by
// (chrom, pos, ref, alt). Pagination uses an opaque key token.
func (s *Store) Query(q PageQuery) (Page, error) {
	if q.Start < 1 {
		q.Start = 1
	}
	if q.PageSize <= 0 {
		q.PageSize = 100
	}
	fromKey := EncodeKey(q.Chrom, q.Start, "")
	// Exclusive upper bound: smallest key strictly greater than every key
	// whose position is <= End. EncodeKey appends "\x00" + ref, so the
	// successor of position End is prefix(End) with its separator bumped
	// past "\x00".
	toKey := q.Chrom + "\x00" + fmt.Sprintf("%0*d", posWidth, q.End) + "\x01"

	token := q.PageToken
	if token != "" {
		raw, err := decodeToken(token)
		if err != nil {
			return Page{}, err
		}
		fromKey = raw
	}

	s.mu.RLock()
	defer s.mu.RUnlock()

	merged := map[string]variant.Record{}
	var order []string

	for key, rec := range s.memtable {
		if key < fromKey || key >= toKey || keyChrom(key) != q.Chrom {
			continue
		}
		order = append(order, key)
		merged[key] = rec
	}

	indexedKeys := s.index.keysInRange(fromKey, toKey, q.Chrom)
	// Keys are already sorted and disjoint from memtable shadowing.
	for _, key := range indexedKeys {
		if _, shadowed := merged[key]; shadowed {
			continue
		}
		seg := s.index.segmentFor(key)
		entry, found, err := seg.lookup(key)
		if err != nil {
			return Page{}, err
		}
		if found {
			order = append(order, key)
			merged[key] = entry.Record
		}
	}

	sort.Strings(order)

	page := Page{Records: make([]variant.Record, 0, q.PageSize)}
	for i, key := range order {
		if i >= q.PageSize {
			page.NextToken = encodeToken(key)
			break
		}
		page.Records = append(page.Records, merged[key])
	}
	return page, nil
}

// Flush forces memtable persistence (used by tests and clean shutdown).
func (s *Store) Flush() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.flushLocked()
}

// Compact merges all segments into one deterministically merged segment.
func (s *Store) Compact() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	if err := s.flushLocked(); err != nil {
		return err
	}
	return s.compactLocked()
}

// Close flushes and releases all resources.
func (s *Store) Close() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.closed {
		return nil
	}
	s.closed = true
	close(s.stopCh)
	if err := s.flushLocked(); err != nil {
		s.wal.Close()
		return err
	}
	if err := s.wal.Close(); err != nil {
		return err
	}
	<-s.doneCh
	return nil
}

// Count returns the number of distinct keys currently stored.
func (s *Store) Count() (int, error) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	seen := map[string]bool{}
	for key := range s.memtable {
		seen[key] = true
	}
	for key := range s.index.segOf {
		seen[key] = true
	}
	return len(seen), nil
}

// ErrClosed is returned after Close.
var ErrClosed = errClosed

var errClosed = fmt.Errorf("store is closed")
