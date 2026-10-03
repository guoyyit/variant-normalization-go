package store

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// manifestName is the generation pointer file stored in the database root.
// Old databases have no manifest; they are read directly from the root and
// opening them never rewrites anything.
const manifestName = "current.json"

// genPrefix names generation directories, e.g. gen-00000000000000000001.
const genPrefix = "gen-"

// manifest is the durable pointer to the active generation directory.
type manifest struct {
	Generation int    `json:"generation"`
	Dir        string `json:"dir"`
}

// Router fronts one active Store. It supports building a new generation in
// the background while queries and imports keep seeing the complete old
// state, then atomically switches readers and writers to the finished
// generation by publishing a new manifest and swapping the in-process
// pointer under a read/write lock.
type Router struct {
	root   string
	mu     sync.RWMutex
	gen    int
	active atomic.Pointer[Store]
}

// OpenRouter opens an existing database without changing its layout.
//
//   - With a manifest, the active store is the generation directory named in
//     it.
//   - Without one (a database written before rebuild support), the active
//     store is the root directory itself, and nothing is written or moved;
//     the first rebuild performs the first switch.
func OpenRouter(root string, mergeFn func(existing, incoming variant.Record) variant.Record, opts ...Option) (*Router, error) {
	r := &Router{root: root}
	if err := os.MkdirAll(root, 0o755); err != nil {
		return nil, err
	}
	m, err := readManifest(root)
	if err != nil {
		return nil, err
	}
	if m == nil {
		st, err := Open(root, mergeFn, opts...)
		if err != nil {
			return nil, err
		}
		r.gen = 0
		r.active.Store(st)
		return r, nil
	}
	st, err := Open(filepath.Join(root, m.Dir), mergeFn, opts...)
	if err != nil {
		return nil, err
	}
	r.gen = m.Generation
	r.active.Store(st)
	return r, nil
}

func readManifest(root string) (*manifest, error) {
	raw, err := os.ReadFile(filepath.Join(root, manifestName))
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	var m manifest
	if err := json.Unmarshal(raw, &m); err != nil {
		return nil, fmt.Errorf("corrupt %s: %w", manifestName, err)
	}
	if m.Dir == "" || filepath.IsAbs(m.Dir) || strings.Contains(m.Dir, "..") {
		return nil, fmt.Errorf("invalid generation dir %q in %s", m.Dir, manifestName)
	}
	if _, err := os.Stat(filepath.Join(root, m.Dir)); err != nil {
		return nil, fmt.Errorf("active generation missing: %w", err)
	}
	return &m, nil
}

func writeManifestAtomic(root string, m manifest) error {
	tmp, err := os.CreateTemp(root, ".current-*.tmp")
	if err != nil {
		return err
	}
	name := tmp.Name()
	enc := json.NewEncoder(tmp)
	if err := enc.Encode(m); err != nil {
		tmp.Close()
		os.Remove(name)
		return err
	}
	if err := tmp.Sync(); err != nil {
		tmp.Close()
		os.Remove(name)
		return err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(name)
		return err
	}
	if err := os.Rename(name, filepath.Join(root, manifestName)); err != nil {
		return err
	}
	d, err := os.Open(root)
	if err != nil {
		return err
	}
	defer d.Close()
	return d.Sync()
}

// Close releases the active store.
func (r *Router) Close() error {
	r.mu.Lock()
	defer r.mu.Unlock()
	return r.active.Load().Close()
}

// Active returns the store currently serving reads and writes. The returned
// pointer is a consistent snapshot of the generation pointer and stays valid
// after a rebuild switches the router to a newer generation.
func (r *Router) Active() *Store {
	return r.active.Load()
}

// ScanActive streams every active-generation record in key order.
func (r *Router) ScanActive(fn func(variant.Record) error) error {
	return r.Active().Scan(fn)
}

// BeginRebuild opens a fresh shadow generation directory. Imports made
// during a rebuild continue to hit the old active store, so queries cannot
// observe a half-written generation.
func (r *Router) BeginRebuild(mergeFn func(existing, incoming variant.Record) variant.Record, opts ...Option) (*Store, string, error) {
	r.mu.RLock()
	next := r.gen + 1
	r.mu.RUnlock()
	dir := filepath.Join(r.root, fmt.Sprintf("%s%020d", genPrefix, next))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, "", err
	}
	st, err := Open(dir, mergeFn, opts...)
	if err != nil {
		return nil, "", err
	}
	return st, dir, nil
}

// CommitRebuild durably publishes shadow as the new active generation, then
// swaps live traffic onto it. The old generation files are retained on disk
// as a backup and are simply no longer referenced.
func (r *Router) CommitRebuild(shadow *Store, dir string) error {
	if err := shadow.Flush(); err != nil {
		return err
	}
	genDir, err := filepath.Rel(r.root, dir)
	if err != nil {
		return err
	}
	next := r.gen + 1
	if err := writeManifestAtomic(r.root, manifest{Generation: next, Dir: genDir}); err != nil {
		return err
	}
	old := r.active.Load()
	r.mu.Lock()
	r.active.Store(shadow)
	r.gen = next
	r.mu.Unlock()
	if old != shadow {
		return old.Close()
	}
	return nil
}

// DiscardRebuild closes and removes an unfinished shadow generation. The
// manifest and the active generation are untouched.
func (r *Router) DiscardRebuild(shadow *Store, dir string) error {
	if err := shadow.Close(); err != nil {
		return err
	}
	return os.RemoveAll(dir)
}

// Generation returns the active generation number (0 for a legacy root).
func (r *Router) Generation() int {
	r.mu.RLock()
	defer r.mu.RUnlock()
	return r.gen
}

// RootDir returns the database root directory.
func (r *Router) RootDir() string { return r.root }
