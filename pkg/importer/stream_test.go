package importer

import (
	"context"
	"io"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
)

// generatedReader emits header then n records without buffering them all in
// a single string: only one line is materialized at a time.
type generatedReader struct {
	n       int
	pos     int
	header  bool
	pending []byte
}

func (g *generatedReader) Read(p []byte) (int, error) {
	if len(g.pending) > 0 {
		n := copy(p, g.pending)
		g.pending = g.pending[n:]
		return n, nil
	}
	if !g.header {
		g.pending = []byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
		g.header = true
		return g.Read(p)
	}
	g.pos++
	if g.pos > g.n {
		return 0, io.EOF
	}
	g.pending = []byte("chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1\n")
	return g.Read(p)
}

func TestImporter_StreamMemoryBounded(t *testing.T) {
	st, err := store.Open(filepath.Join(t.TempDir(), "db"),
		Fold(config.Default()), store.WithMemLimit(64))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	im, err := New(st, refFor(1, strings.Repeat("N", 64)+"CAAAACGT"),
		config.Default(), filepath.Join(t.TempDir(), "cp.json"))
	if err != nil {
		t.Fatal(err)
	}

	var ms1, ms2 runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&ms1)

	rep, err := im.Import(context.Background(), "bulk", "bulk.vcf", &generatedReader{n: 20000})
	if err != nil {
		t.Fatal(err)
	}
	// Same-source idempotence means only the first record contributes; the
	// point here is streaming 20k lines with bounded resident memory.
	if rep.Lines < 20000 {
		t.Fatalf("processed lines=%d", rep.Lines)
	}

	runtime.GC()
	runtime.ReadMemStats(&ms2)
	growth := int64(ms2.HeapAlloc) - int64(ms1.HeapAlloc)
	// Input is ~20k * ~40B ~= 800 KB; assert heap growth is far below input
	// size linearity (segments are on disk, memtable capped at 64 groups).
	if growth > 4<<20 {
		t.Fatalf("heap grew %d bytes, memory not bounded by input", growth)
	}

	count, err := st.Count()
	if err != nil || count != 1 {
		t.Fatalf("count=%d err=%v", count, err)
	}
}

// numberedReader emits n records at distinct positions (distinct keys).
type numberedReader struct {
	n       int
	pos     int
	header  bool
	pending []byte
}

func (g *numberedReader) Read(p []byte) (int, error) {
	if len(g.pending) > 0 {
		n := copy(p, g.pending)
		g.pending = g.pending[n:]
		return n, nil
	}
	if !g.header {
		g.pending = []byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
		g.header = true
		return g.Read(p)
	}
	g.pos++
	if g.pos > g.n {
		return 0, io.EOF
	}
	g.pending = []byte("chr1\t" + strconv.Itoa(g.pos) + "\t.\tA\tT\t.\t.\t.\n")
	return g.Read(p)
}

// seekableNumber wraps numberedReader to satisfy io.ReadSeeker for Resume.
type seekableNumber struct {
	data []byte
	at   int64
}

func (s *seekableNumber) Read(p []byte) (int, error) {
	if s.at >= int64(len(s.data)) {
		return 0, io.EOF
	}
	n := copy(p, s.data[s.at:])
	s.at += int64(n)
	return n, nil
}

func (s *seekableNumber) Seek(offset int64, whence int) (int64, error) {
	switch whence {
	case io.SeekStart:
		s.at = offset
	case io.SeekCurrent:
		s.at += offset
	}
	return s.at, nil
}

func TestImporter_CrashResumeNoDuplicates(t *testing.T) {
	const total = 53
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "db"), Fold(config.Default()),
		store.WithMemLimit(8))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	cpPath := filepath.Join(dir, "cp.json")
	im, err := New(st, refFor(1, strings.Repeat("N", 256)),
		config.Default(), cpPath)
	if err != nil {
		t.Fatal(err)
	}
	im.WithSyncBatch(10)

	// First run: cancel after 25 records emulate a crash mid-batch.
	ctx, cancel := context.WithCancel(context.Background())
	stopped := false
	rr := &cancelReader{n: total, ctx: ctx, cancel: cancel, stopAt: 25, stopped: &stopped}
	_, _ = im.Import(context.Background(), "src", "in.vcf", rr)

	// Resume from checkpoint with a fresh seekable copy of the full input.
	data := buildVCF(total)
	seeker := &seekableNumber{data: data}
	im2, err := New(st, refFor(1, strings.Repeat("N", 256)),
		config.Default(), cpPath)
	if err != nil {
		t.Fatal(err)
	}
	im2.WithSyncBatch(10)
	rep, err := im2.Resume(context.Background(), "src", "in.vcf", seeker)
	if err != nil {
		t.Fatal(err)
	}
	if rep.Imported == 0 {
		t.Fatal("resume imported nothing")
	}

	count, err := st.Count()
	if err != nil {
		t.Fatal(err)
	}
	if count != total {
		t.Fatalf("after resume count=%d want %d (dup or loss)", count, total)
	}
}

type cancelReader struct {
	n       int
	pos     int
	header  bool
	pending []byte
	ctx     context.Context
	cancel  context.CancelFunc
	stopAt  int
	stopped *bool
}

func (c *cancelReader) Read(p []byte) (int, error) {
	if len(c.pending) > 0 {
		n := copy(p, c.pending)
		c.pending = c.pending[n:]
		return n, nil
	}
	if !c.header {
		c.pending = []byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
		c.header = true
		return c.Read(p)
	}
	if c.pos >= c.stopAt && !*c.stopped {
		*c.stopped = true
		c.cancel()
		return 0, io.EOF
	}
	c.pos++
	if c.pos > c.n {
		return 0, io.EOF
	}
	c.pending = []byte("chr1\t" + strconv.Itoa(c.pos) + "\t.\tA\tT\t.\t.\t.\n")
	return c.Read(p)
}

func buildVCF(n int) []byte {
	var b strings.Builder
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i := 1; i <= n; i++ {
		b.WriteString("chr1\t" + strconv.Itoa(i) + "\t.\tA\tT\t.\t.\t.\n")
	}
	return []byte(b.String())
}
