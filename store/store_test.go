package store

import (
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/guoyyit/variant-normalization-go/variant"
)

func rec(chrom string, pos int, ref string, alt string) *variant.Canonical {
	return &variant.Canonical{
		Record:  variant.Record{Chrom: chrom, Pos: pos, Ref: ref, Alt: []string{alt}},
		Sources: []variant.Provenance{{SourceID: "s", Line: int64(pos)}},
	}
}

func TestStore_QueryAndPaging(t *testing.T) {
	dir := t.TempDir()
	s, err := Open(dir)
	if err != nil {
		t.Fatal(err)
	}
	// 两个段，故意乱序，验证归并后全局有序。
	seg1 := []*variant.Canonical{
		rec("chr1", 100, "A", "T"),
		rec("chr1", 300, "G", "A"),
		rec("chr2", 10, "C", "G"),
	}
	seg2 := []*variant.Canonical{
		rec("chr1", 200, "C", "T"),
		rec("chr1", 300, "G", "C"), // 同 pos:ref 不同 alt，应保留为独立记录
	}
	p1, err := s.WriteSortedSegment("seg-1.gob", seg1)
	if err != nil {
		t.Fatal(err)
	}
	p2, err := s.WriteSortedSegment("seg-2.gob", seg2)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Commit([]string{p1, p2}); err != nil {
		t.Fatal(err)
	}

	sn, err := s.Snapshot(context.Background())
	if err != nil {
		t.Fatal(err)
	}
	defer sn.Close()

	got, cursor, err := sn.Query(Query{Chrom: "chr1", StartPos: 100, EndPos: 300, Limit: 2})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 2 || got[0].Record.Pos != 100 || got[1].Record.Pos != 200 {
		t.Fatalf("page1 = %+v", got)
	}
	if cursor == "" {
		t.Fatal("expected cursor")
	}
	got2, cursor2, err := sn.Query(Query{Chrom: "chr1", StartPos: 100, EndPos: 300, Limit: 2, AfterKey: cursor})
	if err != nil {
		t.Fatal(err)
	}
	if len(got2) != 2 || got2[0].Record.Pos != 300 || got2[1].Record.Pos != 300 {
		t.Fatalf("page2 = %+v", got2)
	}
	if got2[0].Record.Alt[0] != "A" || got2[1].Record.Alt[0] != "C" {
		t.Fatalf("same-pos alts not sorted: %s %s", got2[0].Record.Alt, got2[1].Record.Alt)
	}
	if cursor2 != "" {
		t.Fatal("last page cursor must be empty")
	}
}

func TestStore_Persistence(t *testing.T) {
	dir := t.TempDir()
	s, err := Open(dir)
	if err != nil {
		t.Fatal(err)
	}
	p, err := s.WriteSortedSegment("seg.gob", []*variant.Canonical{rec("chr1", 5, "A", "T")})
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Commit([]string{p}); err != nil {
		t.Fatal(err)
	}
	_ = s.Close()

	// 重新打开，清单与数据仍在。
	s2, err := Open(dir)
	if err != nil {
		t.Fatal(err)
	}
	if s2.Generation() != 1 {
		t.Fatalf("generation = %d", s2.Generation())
	}
	sn, _ := s2.Snapshot(context.Background())
	defer sn.Close()
	got, _, err := sn.Query(Query{Chrom: "chr1"})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 1 || got[0].Record.Pos != 5 || got[0].Sources[0].SourceID != "s" {
		t.Fatalf("persisted = %+v", got)
	}
}

func TestStore_SnapshotIsolation(t *testing.T) {
	dir := t.TempDir()
	s, _ := Open(dir)
	p, _ := s.WriteSortedSegment("seg1.gob", []*variant.Canonical{rec("chr1", 1, "A", "T")})
	_ = s.Commit([]string{p})

	sn, _ := s.Snapshot(context.Background())
	defer sn.Close()

	// 提交新版本：旧快照不受影响，仍看到旧版本。
	p2, _ := s.WriteSortedSegment("seg2.gob", []*variant.Canonical{
		rec("chr1", 1, "A", "T"),
		rec("chr1", 2, "C", "G"),
	})
	_ = s.Commit([]string{p2})

	got, _, _ := sn.Query(Query{Chrom: "chr1"})
	if len(got) != 1 {
		t.Fatalf("old snapshot saw new data: %d", len(got))
	}
	sn2, _ := s.Snapshot(context.Background())
	defer sn2.Close()
	got2, _, _ := sn2.Query(Query{Chrom: "chr1"})
	if len(got2) != 2 {
		t.Fatalf("new snapshot = %d", len(got2))
	}
}

func TestStore_OldSnapshotAfterGC(t *testing.T) {
	dir := t.TempDir()
	s, _ := Open(dir)
	p1, _ := s.WriteSortedSegment("committed-1.gob", []*variant.Canonical{rec("chr1", 1, "A", "T")})
	_ = s.Commit([]string{p1})

	old, _ := s.Snapshot(context.Background())

	// 新版本提交会删除旧段文件；旧快照在 Unix 上仍可通过已打开句柄读取。
	p2, _ := s.WriteSortedSegment("committed-2.gob", []*variant.Canonical{
		rec("chr1", 1, "A", "T"),
		rec("chr1", 2, "C", "G"),
	})
	_ = s.Commit([]string{p2})

	got, _, err := old.Query(Query{Chrom: "chr1"})
	if err != nil {
		t.Fatalf("old snapshot unreadable after GC: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("old snapshot = %d records", len(got))
	}
	_ = old.Close()

	if _, err := os.Stat(filepath.Join(dir, "committed-1.gob")); !os.IsNotExist(err) {
		t.Fatalf("old segment should be removed")
	}
}

func TestStore_EmptyBeforeCommit(t *testing.T) {
	dir := t.TempDir()
	s, _ := Open(filepath.Join(dir, "x"))
	sn, err := s.Snapshot(context.Background())
	if err != nil {
		t.Fatal(err)
	}
	got, cursor, err := sn.Query(Query{Chrom: "chr1"})
	if err != nil || len(got) != 0 || cursor != "" {
		t.Fatalf("empty snapshot: %v %d %q", err, len(got), cursor)
	}
}
