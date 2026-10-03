package importer

import (
	"context"
	"path/filepath"
	"strings"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// TestE2E_CrossSourceMergeDifferentForms verifies that the same variant
// reported from different sources in different textual forms merges into a
// single record with recomputed annotations.
func TestE2E_CrossSourceMergeDifferentForms(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "db"), Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	ref := refFor(1, "CAAAACT")
	im, err := New(st, ref, config.Default(), filepath.Join(dir, "cp.json"))
	if err != nil {
		t.Fatal(err)
	}

	// Same deletion of one A in the run, reported in different textual forms.
	srcA := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t3\t.\tAAAC\tAAC\t60\tPASS\t.\tGT:DP\t0/1:40\n")
	srcB := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n" +
		"chr1\t2\t.\tAA\tA\t50\tPASS\t.\tGT:DP\t1/1:30\n")
	srcC := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS3\n" +
		"chr1\t4\t.\tAA\tA\t.\t.\t.\tGT:DP\t0/1:20\n")

	for _, src := range []struct {
		name, file string
		r          *strings.Reader
	}{
		{"lab-a", "a.vcf", srcA},
		{"lab-b", "b.vcf", srcB},
		{"lab-c", "c.vcf", srcC},
	} {
		if _, err := im.Import(context.Background(), src.name, src.file, src.r); err != nil {
			t.Fatal(err)
		}
	}

	// All forms converge to chr1:1 CA>C.
	rec, ok, err := st.Get("chr1", 1, "CA")
	if err != nil || !ok {
		t.Fatalf("canonical record missing: %v", err)
	}
	if rec.Alt[0] != "C" {
		t.Fatalf("alt=%q want C", rec.Alt[0])
	}
	if len(rec.Sources) != 3 {
		t.Fatalf("sources=%d want 3", len(rec.Sources))
	}
	// S1: 0/1 -> AC=1; S2: 1/1 -> AC=2; S3: 0/1 -> AC=1. Total AC=4, AN=6.
	if rec.Info["AC"] != "4" || rec.Info["AN"] != "6" {
		t.Fatalf("AC=%q AN=%q want 4/6", rec.Info["AC"], rec.Info["AN"])
	}
	if rec.Info["DP"] != "90" {
		t.Fatalf("DP=%q want 90", rec.Info["DP"])
	}
}

// TestE2E_IntervalQueryCorrectness verifies that interval queries return the
// correct records based on canonical positions, with no missed or false
// detections due to normalization.
func TestE2E_IntervalQueryCorrectness(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "db"), Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	ref := refFor(1, "CAAAACT")
	im, err := New(st, ref, config.Default(), filepath.Join(dir, "cp.json"))
	if err != nil {
		t.Fatal(err)
	}

	// Report variants at positions that normalization will rewrite.
	// - insertion at 7 T>CT -> canonical 5 A>AC
	// - deletion at 3 AAAC>AAC -> canonical 1 CA>C
	// - SNP at 6 C>G -> canonical 6 C>G
	// - SNP at 7 T>G -> canonical 7 T>G
	vcf := strings.NewReader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t7\t.\tT\tCT\t.\t.\t.\n" +
		"chr1\t3\t.\tAAAC\tAAC\t.\t.\t.\n" +
		"chr1\t6\t.\tC\tG\t.\t.\t.\n" +
		"chr1\t7\t.\tT\tG\t.\t.\t.\n")
	if _, err := im.Import(context.Background(), "src", "in.vcf", vcf); err != nil {
		t.Fatal(err)
	}

	// Query [1,1]: should return the deletion (canonical pos 1).
	page, err := st.Query(store.PageQuery{Chrom: "chr1", Start: 1, End: 1, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 1 {
		t.Fatalf("[1,1] got %d records, want 1 at pos 1", len(page.Records))
	}

	// Query [5,5]: should return the insertion (canonical pos 5).
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 5, End: 5, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 5 {
		t.Fatalf("[5,5] got %d records, want 1 at pos 5", len(page.Records))
	}

	// Query [6,6]: should return the SNP at pos 6.
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 6, End: 6, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 6 {
		t.Fatalf("[6,6] got %d records, want 1 at pos 6", len(page.Records))
	}

	// Query [7,7]: should return the SNP at pos 7.
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 7, End: 7, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 7 {
		t.Fatalf("[7,7] got %d records, want 1 at pos 7", len(page.Records))
	}

	// Query [2,4]: should return nothing (no canonical positions here).
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 2, End: 4, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 0 {
		t.Fatalf("[2,4] got %d records, want 0", len(page.Records))
	}
}

// TestE2E_QuerySortingAndPagination verifies stable sorting and pagination
// across many records.
func TestE2E_QuerySortingAndPagination(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "db"), Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	ref := refFor(1, strings.Repeat("ACGT", 64))
	im, err := New(st, ref, config.Default(), filepath.Join(dir, "cp.json"))
	if err != nil {
		t.Fatal(err)
	}

	// Import SNPs at every position.
	var vcf strings.Builder
	vcf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i := 1; i <= 100; i++ {
		base := string(ref.seq[i-1])
		alt := "A"
		if base == "A" {
			alt = "C"
		}
		vcf.WriteString("chr1\t" + itoa(i) + "\t.\t" + base + "\t" + alt + "\t.\t.\t.\n")
	}
	if _, err := im.Import(context.Background(), "src", "in.vcf", strings.NewReader(vcf.String())); err != nil {
		t.Fatal(err)
	}

	// Paginate with page size 7 and collect all records.
	var all []variant.Record
	token := ""
	for {
		page, err := st.Query(store.PageQuery{
			Chrom: "chr1", Start: 1, End: 100, PageSize: 7, PageToken: token,
		})
		if err != nil {
			t.Fatal(err)
		}
		all = append(all, page.Records...)
		if page.NextToken == "" {
			break
		}
		token = page.NextToken
	}
	if len(all) != 100 {
		t.Fatalf("paginated got %d records, want 100", len(all))
	}

	// Verify stable sorting by (chrom, pos, ref).
	for i := 1; i < len(all); i++ {
		prev := all[i-1]
		cur := all[i]
		if prev.Pos > cur.Pos || (prev.Pos == cur.Pos && prev.Ref > cur.Ref) {
			t.Fatalf("sorting violated at index %d: %d:%s > %d:%s",
				i, prev.Pos, prev.Ref, cur.Pos, cur.Ref)
		}
	}

	// Verify no duplicates.
	seen := map[string]bool{}
	for _, rec := range all {
		key := rec.Chrom + ":" + itoa(rec.Pos) + ":" + rec.Ref
		if seen[key] {
			t.Fatalf("duplicate record: %s", key)
		}
		seen[key] = true
	}
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}
