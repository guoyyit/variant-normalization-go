package importer

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
)

// writeVCF writes a VCF file with the given records to a temp file.
func writeVCF(t *testing.T, dir, name string, records []string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	var b strings.Builder
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	for _, rec := range records {
		b.WriteString(rec + "\n")
	}
	if err := os.WriteFile(path, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestRebuild_ProducesSameAsFreshImport verifies that rebuilding a store
// produces records identical to a fresh import of the same raw data.
func TestRebuild_ProducesSameAsFreshImport(t *testing.T) {
	dir := t.TempDir()
	dbDir := filepath.Join(dir, "db")
	cpPath := filepath.Join(dir, "cp.json")

	records := []string{
		"chr1\t3\t.\tAAAC\tAAC\t60\tPASS\t.\tGT\t0/1",
		"chr1\t7\t.\tT\tCT\t.\t.\t.\tGT\t1/1",
		"chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1",
		"chr1\t6\t.\tC\tG\t.\t.\t.\tGT\t1/1",
	}
	vcfPath := writeVCF(t, dir, "in.vcf", records)

	// Import into the legacy store.
	st, err := store.Open(dbDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	ref := refFor(1, "CAAAACT")
	im, err := New(st, ref, config.Default(), cpPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := im.Import(context.Background(), "src", "in.vcf",
		mustOpen(t, vcfPath)); err != nil {
		t.Fatal(err)
	}
	st.Close()

	// Rebuild the legacy store.
	files := []SourceFile{{Source: "src", File: "in.vcf", Path: vcfPath}}
	if _, err := Rebuild(context.Background(), dbDir, ref, config.Default(),
		cpPath, files); err != nil {
		t.Fatal(err)
	}

	// Fresh import into a new store for comparison.
	freshDir := filepath.Join(dir, "fresh")
	freshSt, err := store.Open(freshDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	freshIM, err := New(freshSt, ref, config.Default(),
		filepath.Join(dir, "fresh-cp.json"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := freshIM.Import(context.Background(), "src", "in.vcf",
		mustOpen(t, vcfPath)); err != nil {
		t.Fatal(err)
	}
	freshSt.Close()

	// Compare every record.
	rebuilt, err := store.Open(dbDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer rebuilt.Close()
	fresh, err := store.Open(freshDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer fresh.Close()

	rebuiltCount, _ := rebuilt.Count()
	freshCount, _ := fresh.Count()
	if rebuiltCount != freshCount {
		t.Fatalf("count mismatch: rebuilt=%d fresh=%d", rebuiltCount, freshCount)
	}

	// Verify each rebuilt record matches the fresh import.
	keys := [][3]any{
		{"chr1", 1, "CA"},
		{"chr1", 5, "A"},
		{"chr1", 6, "C"},
	}
	for _, k := range keys {
		chrom := k[0].(string)
		pos := k[1].(int)
		ref := k[2].(string)
		r1, ok1, err := rebuilt.Get(chrom, pos, ref)
		if err != nil || !ok1 {
			t.Fatalf("rebuilt missing %s:%d %s: %v", chrom, pos, ref, err)
		}
		r2, ok2, err := fresh.Get(chrom, pos, ref)
		if err != nil || !ok2 {
			t.Fatalf("fresh missing %s:%d %s: %v", chrom, pos, ref, err)
		}
		if r1.Pos != r2.Pos || r1.Ref != r2.Ref || r1.Alt[0] != r2.Alt[0] {
			t.Fatalf("record mismatch at %s:%d %s: rebuilt=%d:%s>%s fresh=%d:%s>%s",
				chrom, pos, ref, r1.Pos, r1.Ref, r1.Alt[0], r2.Pos, r2.Ref, r2.Alt[0])
		}
		if r1.Info["AC"] != r2.Info["AC"] || r1.Info["AN"] != r2.Info["AN"] {
			t.Fatalf("stats mismatch at %s:%d %s: rebuilt AC=%s AN=%s fresh AC=%s AN=%s",
				chrom, pos, ref, r1.Info["AC"], r1.Info["AN"], r2.Info["AC"], r2.Info["AN"])
		}
	}
}

// TestRebuild_Resume verifies that a crashed rebuild resumes from the
// checkpoint without duplicates or omissions.
func TestRebuild_Resume(t *testing.T) {
	dir := t.TempDir()
	dbDir := filepath.Join(dir, "db")
	cpPath := filepath.Join(dir, "cp.json")

	// Generate a large VCF with many records.
	var b strings.Builder
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	for i := 1; i <= 10000; i++ {
		b.WriteString("chr1\t" + itoa(i) + "\t.\tA\tT\t.\t.\t.\tGT\t0/1\n")
	}
	vcfPath := writeVCF(t, dir, "in.vcf", nil)
	if err := os.WriteFile(vcfPath, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}

	// Create the legacy store (empty).
	st, err := store.Open(dbDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	st.Close()

	ref := longRefFor(16384)

	// Start the rebuild with a short timeout to simulate a crash.
	ctx, cancel := context.WithTimeout(context.Background(), 1*time.Millisecond)
	defer cancel()
	files := []SourceFile{{Source: "src", File: "in.vcf", Path: vcfPath}}
	_, _ = Rebuild(ctx, dbDir, ref, config.Default(), cpPath, files)

	// Resume the rebuild.
	rep, err := RebuildResume(context.Background(), dbDir, ref, config.Default(),
		cpPath, files)
	if err != nil {
		t.Fatal(err)
	}
	if rep.Imported == 0 {
		t.Fatal("resume imported nothing")
	}

	// Verify the rebuilt store has all 500 records.
	rebuilt, err := store.Open(dbDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	defer rebuilt.Close()
	count, err := rebuilt.Count()
	if err != nil {
		t.Fatal(err)
	}
	if count != 10000 {
		t.Fatalf("count=%d want 10000 (dup or loss)", count)
	}
}

// TestRebuild_ConcurrentQueries verifies that queries during a rebuild see
// complete states (either old or new), never partial.
func TestRebuild_ConcurrentQueries(t *testing.T) {
	dir := t.TempDir()
	dbDir := filepath.Join(dir, "db")
	cpPath := filepath.Join(dir, "cp.json")

	// Generate a large VCF.
	var b strings.Builder
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	for i := 1; i <= 1000; i++ {
		b.WriteString("chr1\t" + itoa(i) + "\t.\tA\tT\t.\t.\t.\tGT\t0/1\n")
	}
	vcfPath := writeVCF(t, dir, "in.vcf", nil)
	if err := os.WriteFile(vcfPath, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}

	// Create the legacy store with some data.
	st, err := store.Open(dbDir, Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	ref := longRefFor(1024)
	im, err := New(st, ref, config.Default(), cpPath)
	if err != nil {
		t.Fatal(err)
	}
	legacyRecords := []string{
		"chr1\t1\t.\tA\tT\t.\t.\t.\tGT\t0/1",
		"chr1\t2\t.\tA\tT\t.\t.\t.\tGT\t0/1",
	}
	legacyPath := writeVCF(t, dir, "legacy.vcf", legacyRecords)
	if _, err := im.Import(context.Background(), "legacy", "legacy.vcf",
		mustOpen(t, legacyPath)); err != nil {
		t.Fatal(err)
	}
	st.Close()

	// Start the rebuild in the background.
	done := make(chan error, 1)
	go func() {
		files := []SourceFile{{Source: "src", File: "in.vcf", Path: vcfPath}}
		_, err := Rebuild(context.Background(), dbDir, ref, config.Default(),
			cpPath, files)
		done <- err
	}()

	// Run queries concurrently and verify they see complete states.
	queriesOK := 0
	for i := 0; i < 20; i++ {
		q, err := store.Open(dbDir, Fold(config.Default()))
		if err != nil {
			t.Fatal(err)
		}
		page, err := q.Query(store.PageQuery{Chrom: "chr1", Start: 1, End: 1000, PageSize: 1000})
		q.Close()
		if err != nil {
			t.Fatal(err)
		}
		// A complete state has either 2 (legacy) or 1000 (rebuilt) records.
		if len(page.Records) != 2 && len(page.Records) != 1000 {
			t.Fatalf("partial state observed: %d records", len(page.Records))
		}
		queriesOK++
	}

	// Wait for the rebuild to finish.
	if err := <-done; err != nil {
		t.Fatal(err)
	}
	if queriesOK != 20 {
		t.Fatalf("only %d/20 queries saw complete states", queriesOK)
	}
}

func mustOpen(t *testing.T, path string) *os.File {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	return f
}

// longRefFor creates a reference of the given length with all "A" bases.
func longRefFor(length int) fakeRef {
	return fakeRef{seq: strings.Repeat("A", length)}
}
