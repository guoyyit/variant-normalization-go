package store

import (
	"encoding/gob"
	"os"
	"path/filepath"

	"github.com/guoyyit/variant-normalization-go/variant"
)

// WriteSortedSegment 将一组已按 Key 排序的 Canonical 写入一个段文件。
// 文件名由调用方指定（位于存储目录），写完执行 fsync 以保证崩溃安全。
func (s *Store) WriteSortedSegment(name string, recs []*variant.Canonical) (string, error) {
	path := s.segmentPath(name)
	tmp, err := os.CreateTemp(s.dir, ".seg-*")
	if err != nil {
		return "", err
	}
	tmpName := tmp.Name()
	enc := gob.NewEncoder(tmp)
	for _, c := range recs {
		if err := enc.Encode(c); err != nil {
			_ = tmp.Close()
			_ = os.Remove(tmpName)
			return "", err
		}
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		_ = os.Remove(tmpName)
		return "", err
	}
	if err := tmp.Close(); err != nil {
		_ = os.Remove(tmpName)
		return "", err
	}
	if err := os.Rename(tmpName, path); err != nil {
		return "", err
	}
	return path, nil
}

// SegmentWriter 流式写一个有序段（内存不随输入总量增长，由外部分块调用）。
type SegmentWriter struct {
	f   *os.File
	enc *gob.Encoder
	n   int
}

// NewSegmentWriter 在目录中创建临时段，Close 后原子改名为 name。
func (s *Store) NewSegmentWriter(name string) (*SegmentWriter, error) {
	tmp, err := os.CreateTemp(s.dir, ".seg-*")
	if err != nil {
		return nil, err
	}
	return &SegmentWriter{f: tmp, enc: gob.NewEncoder(tmp)}, nil
}

func (w *SegmentWriter) Write(c *variant.Canonical) error {
	if err := w.enc.Encode(c); err != nil {
		return err
	}
	w.n++
	return nil
}

// Len 返回已写入条数。
func (w *SegmentWriter) Len() int { return w.n }

// Finish fsync、关闭并把临时段原子改名为目标段。
func (w *SegmentWriter) Finish(s *Store, name string) (string, error) {
	if err := w.f.Sync(); err != nil {
		_ = w.f.Close()
		return "", err
	}
	if err := w.f.Close(); err != nil {
		return "", err
	}
	path := s.segmentPath(name)
	return path, os.Rename(w.f.Name(), path)
}

// Commit 用新段集合原子替换当前已提交版本。
func (s *Store) Commit(segNames []string) error {
	names := make([]string, len(segNames))
	for i, n := range segNames {
		names[i] = filepath.Base(n)
	}
	return s.commit(names)
}

// CurrentSegmentPaths 返回当前已提交版本的全部段路径。
func (s *Store) CurrentSegmentPaths() []string { return s.currentSegments() }

// Dir 返回存储目录。
func (s *Store) Dir() string { return s.dir }

// SegmentPath 返回目录内某个段文件的绝对路径。
func (s *Store) SegmentPath(name string) string { return s.segmentPath(name) }
