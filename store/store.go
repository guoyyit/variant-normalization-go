// Package store 提供磁盘持久化、快照隔离的规范记录存储。
//
// 每个已提交版本由若干按 Key（chrom,pos,ref,alt）排序的“段”文件组成；
// manifest.json 原子指向当前版本。查询通过对所有段做 k 路归并获得一致视图，
// 写入在新段上进行并在提交时一次性切换清单，因此并发查询绝不会看到
// 未完成的合并中间状态。
package store

import (
	"context"
	"encoding/gob"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"
)

const manifestName = "manifest.json"

// ErrNotFound 表示数据目录尚未初始化。
var ErrNotFound = errors.New("store not initialized")

// segmentDesc 描述一个已提交的有序段。
type segmentDesc struct {
	Name string `json:"name"`
}

// manifest 是当前已提交版本的清单。
type manifest struct {
	Generation int64         `json:"generation"`
	Segments   []segmentDesc `json:"segments"`
}

// Store 为磁盘持久化、快照隔离的规范记录存储。
type Store struct {
	dir string
	mu  sync.RWMutex
	cur manifest
}

// Open 打开（或创建）数据目录中的存储。
func Open(dir string) (*Store, error) {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}
	s := &Store{dir: dir}
	if err := s.loadManifest(); err != nil && !errors.Is(err, os.ErrNotExist) {
		return nil, err
	}
	return s, nil
}

func (s *Store) manifestPath() string { return filepath.Join(s.dir, manifestName) }

func (s *Store) loadManifest() error {
	data, err := os.ReadFile(s.manifestPath())
	if err != nil {
		return err
	}
	var m manifest
	if err := json.Unmarshal(data, &m); err != nil {
		return fmt.Errorf("decode manifest: %w", err)
	}
	sort.Slice(m.Segments, func(i, j int) bool { return m.Segments[i].Name < m.Segments[j].Name })
	s.cur = m
	return nil
}

// Generation 返回当前已提交版本号。
func (s *Store) Generation() int64 {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.cur.Generation
}

// Snapshot 返回当前已提交版本的只读快照。
func (s *Store) Snapshot(ctx context.Context) (Snapshot, error) {
	s.mu.RLock()
	m := s.cur
	s.mu.RUnlock()
	paths := make([]string, len(m.Segments))
	for i, seg := range m.Segments {
		paths[i] = filepath.Join(s.dir, seg.Name)
	}
	return newSnapshot(ctx, m.Generation, paths)
}

// currentSegments 在锁内返回当前段文件名副本（供提交者读取旧版本）。
func (s *Store) currentSegments() []string {
	s.mu.RLock()
	defer s.mu.RUnlock()
	out := make([]string, len(s.cur.Segments))
	for i, seg := range s.cur.Segments {
		out[i] = filepath.Join(s.dir, seg.Name)
	}
	return out
}

// currentGeneration 在锁内返回当前版本号。
func (s *Store) currentGeneration() int64 {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.cur.Generation
}

// commit 以新段列表原子替换当前版本；新清单 fsync 后重命名生效。
func (s *Store) commit(segNames []string) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	next := manifest{Generation: s.cur.Generation + 1}
	for _, name := range segNames {
		next.Segments = append(next.Segments, segmentDesc{Name: filepath.Base(name)})
	}
	sort.Slice(next.Segments, func(i, j int) bool { return next.Segments[i].Name < next.Segments[j].Name })

	tmp, err := os.CreateTemp(s.dir, ".manifest-*")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	enc := json.NewEncoder(tmp)
	enc.SetIndent("", "  ")
	if err := enc.Encode(next); err != nil {
		_ = tmp.Close()
		_ = os.Remove(tmpName)
		return err
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		_ = os.Remove(tmpName)
		return err
	}
	if err := tmp.Close(); err != nil {
		_ = os.Remove(tmpName)
		return err
	}
	if err := os.Rename(tmpName, s.manifestPath()); err != nil {
		return err
	}
	// 清单已切换：清理不再属于任何版本的旧段文件。
	keep := map[string]bool{}
	for _, seg := range next.Segments {
		keep[seg.Name] = true
	}
	entries, _ := os.ReadDir(s.dir)
	for _, e := range entries {
		name := e.Name()
		if !keep[name] && isSegmentFile(name) {
			_ = os.Remove(filepath.Join(s.dir, name))
		}
	}
	s.cur = next
	return nil
}

func isSegmentFile(name string) bool {
	return (strings.HasPrefix(name, "committed-") ||
		strings.HasPrefix(name, "incoming-") ||
		strings.HasPrefix(name, "seg-")) &&
		strings.HasSuffix(name, ".gob")
}

// segmentPath 返回目录内一个段文件的绝对路径。
func (s *Store) segmentPath(name string) string {
	return filepath.Join(s.dir, name)
}

// Close 关闭存储（当前实现无独占资源，保留接口以便扩展）。
func (s *Store) Close() error { return nil }

// ---- 段文件编码（gob 自带帧分隔，逐条 Encode/Decode） ----

func newDecoder(r io.Reader) *gob.Decoder { return gob.NewDecoder(r) }
