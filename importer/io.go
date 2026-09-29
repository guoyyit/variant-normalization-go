package importer

import (
	"bufio"
	"encoding/gob"
	"io"
	"os"

	"github.com/guoyyit/variant-normalization-go/variant"
)

// entry 是一个已打开段的顺序迭代器。
type entry struct {
	f   *os.File
	dec *gob.Decoder
	cur *variant.Canonical
	err error
}

func openEntry(path string) (*entry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	e := &entry{f: f, dec: gob.NewDecoder(bufio.NewReaderSize(f, 64*1024))}
	e.next()
	return e, nil
}

func (e *entry) next() {
	var c variant.Canonical
	if err := e.dec.Decode(&c); err != nil {
		e.cur = nil
		if err != io.EOF {
			e.err = err
		}
		return
	}
	e.cur = &c
}

func (e *entry) close() { _ = e.f.Close() }

// entries 管理多路段的最小堆归并。
type entries struct {
	list []*entry
}

func openEntries(paths []string) (*entries, error) {
	es := &entries{}
	for _, p := range paths {
		e, err := openEntry(p)
		if os.IsNotExist(err) {
			continue
		}
		if err != nil {
			es.close()
			return nil, err
		}
		if e.cur != nil {
			es.push(e)
		} else {
			e.close()
			if e.err != nil {
				es.close()
				return nil, e.err
			}
		}
	}
	return es, nil
}

func (es *entries) push(e *entry) {
	// 简单插入堆（段数通常很少；用 container/heap 也可，这里保持依赖最小）。
	es.list = append(es.list, e)
	es.up(len(es.list) - 1)
}

func (es *entries) pop() *entry {
	n := len(es.list)
	top := es.list[0]
	es.list[0] = es.list[n-1]
	es.list = es.list[:n-1]
	if len(es.list) > 0 {
		es.down(0)
	}
	return top
}

func (es *entries) less(i, j int) bool {
	return es.list[i].cur.Record.Key() < es.list[j].cur.Record.Key()
}

func (es *entries) swap(i, j int) { es.list[i], es.list[j] = es.list[j], es.list[i] }

func (es *entries) up(i int) {
	for i > 0 {
		p := (i - 1) / 2
		if !es.less(i, p) {
			return
		}
		es.swap(i, p)
		i = p
	}
}

func (es *entries) down(i int) {
	n := len(es.list)
	for {
		l, r, best := 2*i+1, 2*i+2, i
		if l < n && es.less(l, best) {
			best = l
		}
		if r < n && es.less(r, best) {
			best = r
		}
		if best == i {
			return
		}
		es.swap(i, best)
		i = best
	}
}

// channel 按全局 Key 顺序产出所有段中的记录。
func (es *entries) channel() func(func(*variant.Canonical, error) bool) {
	return func(yield func(*variant.Canonical, error) bool) {
		for len(es.list) > 0 {
			e := es.list[0]
			c := e.cur
			if !yield(c, nil) {
				return
			}
			es.pop()
			e.next()
			if e.err != nil {
				yield(nil, e.err)
				return
			}
			if e.cur != nil {
				es.push(e)
			} else {
				e.close()
			}
		}
	}
}

func (es *entries) close() {
	for _, e := range es.list {
		e.close()
	}
	es.list = nil
}
