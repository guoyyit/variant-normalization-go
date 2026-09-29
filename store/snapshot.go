package store

import (
	"bufio"
	"container/heap"
	"context"
	"encoding/gob"
	"io"
	"os"

	"github.com/guoyyit/variant-normalization-go/variant"
)

// Query 为染色体位置区间查询条件。
type Query struct {
	Chrom    string
	StartPos int
	EndPos   int
	// AfterKey 为上一页最后一条记录的 Key，用于稳定翻页。
	AfterKey string
	Limit    int
}

// Snapshot 是某个已提交版本的只读视图。
type Snapshot interface {
	// Query 按染色体与位置区间查询，结果稳定排序并支持翻页。
	Query(q Query) ([]*variant.Canonical, string, error)
	Generation() int64
	Close() error
}

// snapshot 对一组有序段做 k 路归并。
type snapshot struct {
	gen   int64
	files []*os.File
}

func newSnapshot(ctx context.Context, gen int64, paths []string) (*snapshot, error) {
	sn := &snapshot{gen: gen}
	for _, p := range paths {
		f, err := os.Open(p)
		if err != nil {
			_ = sn.Close()
			return nil, err
		}
		sn.files = append(sn.files, f)
	}
	return sn, nil
}

func (s *snapshot) Generation() int64 { return s.gen }

func (s *snapshot) Close() error {
	var first error
	for _, f := range s.files {
		if err := f.Close(); err != nil && first == nil {
			first = err
		}
	}
	s.files = nil
	return first
}

// segmentIter 顺序读取一个有序段。
type segmentIter struct {
	br  *bufio.Reader
	dec *gob.Decoder
	cur *variant.Canonical
	err error
	idx int
}

func newSegmentIter(r io.Reader, idx int) *segmentIter {
	it := &segmentIter{idx: idx}
	it.br = bufio.NewReaderSize(r, 64*1024)
	it.dec = newDecoder(it.br)
	it.advance()
	return it
}

func (it *segmentIter) advance() {
	var c variant.Canonical
	if err := it.dec.Decode(&c); err != nil {
		it.cur = nil
		if err != io.EOF {
			it.err = err
		}
		return
	}
	it.cur = &c
}

// mergeHeap 按 Key 归并多个段；同 Key 时段下标小者优先（确定且稳定）。
type mergeHeap []*segmentIter

func (h mergeHeap) Len() int { return len(h) }
func (h mergeHeap) Less(i, j int) bool {
	ri, rj := &h[i].cur.Record, &h[j].cur.Record
	// 全局排序：染色体、坐标、REF、ALT（完整 Key 已含定宽坐标）。
	if ki, kj := ri.Key(), rj.Key(); ki != kj {
		return ki < kj
	}
	return h[i].idx < h[j].idx
}
func (h mergeHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] }
func (h *mergeHeap) Push(x any)   { *h = append(*h, x.(*segmentIter)) }
func (h *mergeHeap) Pop() any {
	old := *h
	n := len(old)
	it := old[n-1]
	*h = old[:n-1]
	return it
}

// Query 执行区间查询并返回本页结果与下一页游标（最后一页为空串）。
func (s *snapshot) Query(q Query) ([]*variant.Canonical, string, error) {
	h := &mergeHeap{}
	for i, f := range s.files {
		// 每次查询使用独立 SectionReader，避免共享文件偏移并保证从段首开始。
		sr := io.NewSectionReader(f, 0, 1<<63-1)
		it := newSegmentIter(sr, i)
		if it.err != nil {
			return nil, "", it.err
		}
		if it.cur != nil {
			heap.Push(h, it)
		}
	}

	var out []*variant.Canonical
	var lastKey string
	emit := func(c *variant.Canonical) bool {
		if q.Limit > 0 && len(out) >= q.Limit {
			return false
		}
		out = append(out, c)
		lastKey = c.Record.Key()
		return true
	}

	pop := func() error {
		it := heap.Pop(h).(*segmentIter)
		it.advance()
		if it.err != nil {
			return it.err
		}
		if it.cur != nil {
			heap.Push(h, it)
		}
		return nil
	}

	for h.Len() > 0 {
		c := (*h)[0].cur
		key := c.Record.Key()

		// 1) 先推进翻页游标：跳过所有 key <= AfterKey 的记录，不做区间判断。
		if q.AfterKey != "" && key <= q.AfterKey {
			if err := pop(); err != nil {
				return nil, "", err
			}
			continue
		}

		// 2) 全局有序，越过目标染色体即可结束。
		if c.Record.Chrom > q.Chrom {
			break
		}
		// 3) 区间过滤（POS 严格按位置，跨染色体不影响）。
		inRange := c.Record.Chrom == q.Chrom &&
			(q.StartPos == 0 || c.Record.Pos >= q.StartPos) &&
			(q.EndPos == 0 || c.Record.Pos <= q.EndPos)
		if inRange && !emit(c) {
			// 达到页大小：不消费堆顶，下一页用游标从此记录继续。
			return out, lastKey, nil
		}
		if err := pop(); err != nil {
			return nil, "", err
		}
	}
	// 堆耗尽：本页之后已无任何记录，返回空游标表示结束。
	return out, "", nil
}
