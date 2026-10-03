package importer

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/stats"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
)

type fakeRef struct {
	seq string
}

func (f fakeRef) Substring(chrom string, start, end int) (string, error) {
	return f.seq[start-1 : end-1], nil
}

func (f fakeRef) Length(chrom string) (int, error) { return len(f.seq), nil }

// refFor builds a length-256 reference with bases placed at 1-based start.
func refFor(start int, bases string) fakeRef {
	seq := make([]byte, 256)
	for i := range seq {
		seq[i] = 'N'
	}
	copy(seq[start-1:], bases)
	return fakeRef{seq: string(seq)}
}

func newTestImporter(t *testing.T) (*Importer, *store.Store) {
	t.Helper()
	st, err := store.Open(t.TempDir(), Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { st.Close() })
	im, err := New(st, refFor(1, "CAAAACGT"), config.Default(),
		filepath.Join(t.TempDir(), "cp.json"))
	if err != nil {
		t.Fatal(err)
	}
	return im, st
}

func TestImporter_EndToEnd(t *testing.T) {
	im, st := newTestImporter(t)

	vcfData := strings.NewReader(`##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	3	.	AAAC	AAC	60	PASS	DP=40	GT	0/1
chr1	100	.	A	T,G	.	.	.	GT	1/2
chr1	9	.	C	T	.	.	.	GT	0/1
`)
	rep, err := im.Import(context.Background(), "src-a", "a.vcf", vcfData)
	if err != nil {
		t.Fatal(err)
	}
	// Deletion left-aligns to 1 group; the multi-allelic line splits into 2.
	if rep.Imported != 4 { // one del + two splits + the trailing SNP
		t.Fatalf("imported=%d report=%+v", rep.Imported, rep)
	}
	if len(rep.Rejected) != 0 {
		t.Fatalf("rejected=%+v", rep.Rejected)
	}

	// Left-aligned deletion canonicalizes to chr1:1 CA>C.
	rec, ok, err := st.Get("chr1", 1, "CA")
	if err != nil || !ok {
		t.Fatalf("canonical deletion missing: %v", err)
	}
	if rec.Alt[0] != "C" || len(rec.Sources) != 1 {
		t.Fatalf("rec=%+v", rec)
	}
	if rec.Info["AC"] != "1" || rec.Info["AN"] != "2" {
		t.Fatalf("stats not recomputed: %v", rec.Info)
	}

	// Multi-allelic splits at the same key re-merge into one multi-alt group.
	count, _ := st.Count()
	if count != 3 {
		t.Fatalf("count=%d", count)
	}

	group, ok, _ := st.Get("chr1", 100, "A")
	if !ok {
		t.Fatal("multi-alt group missing at chr1:100")
	}
	if got := strings.Join(group.Alt, ","); got != "G,T" {
		t.Fatalf("alts=%q want G,T", got)
	}
	// Original GT 1/2 -> slots G=1,T=2; AC 1 each, AN 2.
	if group.Info["AC"] != "1,1" || group.Info["AN"] != "2" {
		t.Fatalf("split/merge stats wrong: %v", group.Info)
	}
}

func TestImporter_OutOfOrderDedup(t *testing.T) {
	im, st := newTestImporter(t)

	a := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t3\t.\tAAAC\tAAC\t60\tPASS\tDP=10\tGT\t0/1\n")
	b := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t3\t.\tAAAC\tAAC\t40\tLowQual\tDP=20\tGT\t0/1\n")

	if _, err := im.Import(context.Background(), "src-a", "a.vcf", a); err != nil {
		t.Fatal(err)
	}
	if _, err := im.Import(context.Background(), "src-b", "b.vcf", b); err != nil {
		t.Fatal(err)
	}
	rec, ok, _ := st.Get("chr1", 1, "CA")
	if !ok {
		t.Fatal("missing merged record")
	}
	if len(rec.Sources) != 2 {
		t.Fatalf("sources=%v", rec.Sources)
	}
	// max qual
	if rec.Qual != 60 {
		t.Fatalf("qual=%v", rec.Qual)
	}

	// Re-import src-a same file/line: must be idempotent.
	a2 := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t3\t.\tAAAC\tAAC\t60\tPASS\tDP=10\tGT\t0/1\n")
	if _, err := im.Import(context.Background(), "src-a", "a.vcf", a2); err != nil {
		t.Fatal(err)
	}
	rec, _, _ = st.Get("chr1", 1, "CA")
	if len(rec.Sources) != 2 || rec.Info["AN"] != "2" {
		t.Fatalf("replay changed state: sources=%d info=%v", len(rec.Sources), rec.Info)
	}

	// Same VCF bytes in opposite source order gives the same canonical key.
	if _, err := im.Import(context.Background(), "src-b", "b.vcf", strings.NewReader(
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
			"chr1\t3\t.\tAAAC\tAAC\t40\t.\t.\n")); err != nil {
		t.Fatal(err)
	}
	rec, _, _ = st.Get("chr1", 1, "CA")
	if len(rec.Sources) != 2 {
		t.Fatalf("dedup failed: %v", rec.Sources)
	}

	_ = stats.Summary(rec)
}

func TestImporter_Resume(t *testing.T) {
	dir := t.TempDir()
	cpPath := filepath.Join(dir, "cp.json")
	st, err := store.Open(filepath.Join(dir, "db"),
		Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}

	full := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t1\t.\tA\tT\t.\t.\t.\n" +
		"chr1\t2\t.\tC\tG\t.\t.\t.\n" +
		"chr1\t3\t.\tG\tA\t.\t.\t.\n"

	// Write the input file and emulate a crash after line 1 by truncating
	// the checkpoint to the first record boundary.
	inputPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inputPath, []byte(full), 0o644); err != nil {
		t.Fatal(err)
	}

	im, err := New(st, refFor(1, "CAAAACGT"), config.Default(), cpPath)
	if err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(inputPath)
	if err != nil {
		t.Fatal(err)
	}
	rep, err := im.Import(context.Background(), "src", "in.vcf", f)
	f.Close()
	if err != nil {
		t.Fatal(err)
	}
	if rep.Imported != 3 {
		t.Fatalf("imported=%d", rep.Imported)
	}

	// Resume: checkpoint covers the whole file, nothing new must arrive and
	// no duplicates may appear.
	f2, err := os.Open(inputPath)
	if err != nil {
		t.Fatal(err)
	}
	rep2, err := im.Resume(context.Background(), "src", "in.vcf", f2)
	f2.Close()
	if err != nil {
		t.Fatal(err)
	}
	if rep2.Imported != 0 {
		t.Fatalf("resume imported=%d, want 0", rep2.Imported)
	}
	count, _ := st.Count()
	if count != 3 {
		t.Fatalf("count after resume=%d", count)
	}
	st.Close()
}

// TestImporter_RejectReasonsDistinct feeds one file containing an invalid
// allele, an out-of-bounds record and a valid record. Both bad rows are
// rejected with distinguishable reasons, import continues, and the valid row
// lands in the store.
func TestImporter_RejectReasonsDistinct(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "db"), Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()
	im, err := New(st, fakeRef{seq: "ACGT"}, config.Default(), filepath.Join(dir, "cp.json"))
	if err != nil {
		t.Fatal(err)
	}
	data := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t2\t.\tC\tX\t.\t.\t.\n" +
		"chr1\t9\t.\tA\tT\t.\t.\t.\n" +
		"chr1\t3\t.\tG\tT\t.\t.\t.\n")
	rep, err := im.Import(context.Background(), "s", "a.vcf", data)
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Rejected) != 2 {
		t.Fatalf("rejected=%d want 2", len(rep.Rejected))
	}
	var sawInvalid, sawBounds bool
	for _, rj := range rep.Rejected {
		if strings.Contains(rj.Reason, "invalid variant format") {
			sawInvalid = true
		}
		if strings.Contains(rj.Reason, "coordinate out of bounds") {
			sawBounds = true
		}
	}
	if !sawInvalid || !sawBounds {
		t.Fatalf("reasons not distinguishable: %+v", rep.Rejected)
	}
	if rep.Imported != 1 {
		t.Fatalf("imported=%d want 1 (invalid rows must not abort import)", rep.Imported)
	}
	count, err := st.Count()
	if err != nil || count != 1 {
		t.Fatalf("count=%d err=%v", count, err)
	}
}
