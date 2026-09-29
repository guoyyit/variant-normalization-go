package store

import (
	"fmt"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func testMerge(existing, incoming variant.Record) variant.Record {
	if len(existing.Sources) == 0 {
		return incoming.Clone()
	}
	out := existing.Clone()
	out.Sources = append(out.Sources, incoming.Sources...)
	out.Alt = append(out.Alt, incoming.Alt...)
	return out
}

func newTestStore(t *testing.T, memLimit int) *Store {
	t.Helper()
	dir := t.TempDir()
	st, err := Open(dir, testMerge, WithMemLimit(memLimit))
	if err != nil {
		t.Fatal(err)
	}
	return st
}

func mkRec(chrom string, pos int, ref, alt, source string) variant.Record {
	return variant.Record{
		Chrom: chrom, Pos: pos, Ref: ref, Alt: []string{alt},
		Sources: []variant.SourceRef{{Source: source, File: source + ".vcf"}},
	}
}

func TestStore_UpsertAndQuery(t *testing.T) {
	st := newTestStore(t, 4)
	for i := 1; i <= 10; i++ {
		if err := st.Upsert(mkRec("chr1", i, "A", "T", "s1")); err != nil {
			t.Fatal(err)
		}
	}

	page, err := st.Query(PageQuery{Chrom: "chr1", Start: 3, End: 7, PageSize: 2})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 2 {
		t.Fatalf("page size = %d", len(page.Records))
	}
	if page.Records[0].Pos != 3 || page.Records[1].Pos != 4 {
		t.Fatal("page order wrong")
	}

	got := []int{}
	token := page.NextToken
	for token != "" {
		page, err = st.Query(PageQuery{Chrom: "chr1", Start: 3, End: 7, PageSize: 2, PageToken: token})
		if err != nil {
			t.Fatal(err)
		}
		for _, rec := range page.Records {
			got = append(got, rec.Pos)
		}
		token = page.NextToken
	}
	if fmt.Sprint(got) != "[5 6 7]" {
		t.Fatalf("pagination = %v", got)
	}
}

func TestStore_IntervalUpperBound(t *testing.T) {
	st := newTestStore(t, 2)
	if err := st.Upsert(mkRec("chr1", 5, "ZZ", "T", "s")); err != nil {
		t.Fatal(err)
	}
	if err := st.Upsert(mkRec("chr1", 6, "A", "T", "s")); err != nil {
		t.Fatal(err)
	}
	page, err := st.Query(PageQuery{Chrom: "chr1", Start: 5, End: 5, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 5 || page.Records[0].Ref != "ZZ" {
		t.Fatalf("boundary query wrong: %+v", page.Records)
	}
}

func TestStore_FlushShadowingAndMerge(t *testing.T) {
	dir := t.TempDir()
	st, err := Open(dir, testMerge, WithMemLimit(2))
	if err != nil {
		t.Fatal(err)
	}
	// First flush: key v1 persisted.
	if err := st.Upsert(mkRec("chr1", 1, "A", "T", "first")); err != nil {
		t.Fatal(err)
	}
	if err := st.Upsert(mkRec("chr1", 2, "A", "T", "first")); err != nil {
		t.Fatal(err)
	}
	// Second flush in another segment updates key1.
	if err := st.Upsert(mkRec("chr1", 1, "A", "G", "second")); err != nil {
		t.Fatal(err)
	}
	if err := st.Upsert(mkRec("chr1", 3, "A", "T", "second")); err != nil {
		t.Fatal(err)
	}

	rec, ok, err := st.Get("chr1", 1, "A")
	if err != nil || !ok {
		t.Fatal(err)
	}
	// Merge of v1 + update: sources from both, alts T and G.
	if len(rec.Sources) != 2 || len(rec.Alt) != 2 {
		t.Fatalf("record = %+v", rec)
	}

	if err := st.Compact(); err != nil {
		t.Fatal(err)
	}
	rec, ok, err = st.Get("chr1", 1, "A")
	if err != nil || !ok {
		t.Fatal(err)
	}
	if len(rec.Sources) != 2 || len(rec.Alt) != 2 {
		t.Fatalf("after compact record = %+v", rec)
	}

	count, err := st.Count()
	if err != nil || count != 3 {
		t.Fatalf("count=%d err=%v", count, err)
	}
	if err := st.Close(); err != nil {
		t.Fatal(err)
	}

	// Reopen: state must survive.
	st2, err := Open(dir, testMerge)
	if err != nil {
		t.Fatal(err)
	}
	defer st2.Close()
	count, err = st2.Count()
	if err != nil || count != 3 {
		t.Fatalf("after reopen count=%d err=%v", count, err)
	}
	rec, ok, _ = st2.Get("chr1", 1, "A")
	if !ok || len(rec.Sources) != 2 {
		t.Fatalf("after reopen record = %+v", rec)
	}
}

func TestStore_WALReplay(t *testing.T) {
	dir := t.TempDir()
	st, err := Open(dir, testMerge, WithMemLimit(1000))
	if err != nil {
		t.Fatal(err)
	}
	if err := st.Upsert(mkRec("chr2", 5, "C", "G", "w")); err != nil {
		t.Fatal(err)
	}
	// Close without flushing: WAL holds the record. But Close flushes, so
	// simulate a crash by dropping the handle directly.
	st.mu.Lock()
	_ = st.wal.Close()
	st.mu.Unlock()

	// After a crash-like shutdown, reopen replays WAL into memtable.
	st2, err := Open(dir, testMerge)
	if err != nil {
		t.Fatal(err)
	}
	defer st2.Close()
	rec, ok, err := st2.Get("chr2", 5, "C")
	if err != nil || !ok || rec.Alt[0] != "G" {
		t.Fatalf("wal replay failed: %+v %v", rec, err)
	}
}
