package store

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func identityFold(_, in variant.Record) variant.Record { return in }

func TestRouter_LegacyThenRebuildGenerations(t *testing.T) {
	root := filepath.Join(t.TempDir(), "db")

	// Write a legacy database directly in the root (no manifest).
	old, err := Open(root, identityFold)
	if err != nil {
		t.Fatal(err)
	}
	if err := old.Upsert(variant.Record{Chrom: "c", Pos: 1, Ref: "A", Alt: []string{"T"}}); err != nil {
		t.Fatal(err)
	}
	if err := old.Close(); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(root, manifestName)); !os.IsNotExist(err) {
		t.Fatalf("legacy db must have no manifest, stat err=%v", err)
	}

	r, err := OpenRouter(root, identityFold)
	if err != nil {
		t.Fatal(err)
	}
	if r.Generation() != 0 || r.Active() == nil {
		t.Fatalf("legacy open gen=%d", r.Generation())
	}
	// Querying must not create a manifest.
	if _, err := r.Active().Query(PageQuery{Chrom: "c", Start: 1, End: 1}); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(root, manifestName)); !os.IsNotExist(err) {
		t.Fatal("opening/querying legacy db silently created a manifest")
	}

	shadow, dir, err := r.BeginRebuild(identityFold)
	if err != nil {
		t.Fatal(err)
	}
	if err := shadow.Upsert(variant.Record{Chrom: "c", Pos: 2, Ref: "A", Alt: []string{"G"}}); err != nil {
		t.Fatal(err)
	}
	// Before commit the active generation is unchanged.
	if n, _ := r.Active().Count(); n != 1 {
		t.Fatalf("active changed during rebuild: count=%d", n)
	}
	if err := r.CommitRebuild(shadow, dir); err != nil {
		t.Fatal(err)
	}
	if r.Generation() != 1 {
		t.Fatalf("gen=%d want 1", r.Generation())
	}
	if n, _ := r.Active().Count(); n != 1 {
		t.Fatalf("new generation count=%d want 1", n)
	}
	page, err := r.Active().Query(PageQuery{Chrom: "c", Start: 1, End: 10, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 2 {
		t.Fatalf("rebuilt content wrong: %+v", page.Records)
	}
	if err := r.Close(); err != nil {
		t.Fatal(err)
	}

	// Reopen follows the manifest.
	r2, err := OpenRouter(root, identityFold)
	if err != nil {
		t.Fatal(err)
	}
	defer r2.Close()
	if r2.Generation() != 1 {
		t.Fatalf("reopen gen=%d want 1", r2.Generation())
	}
}
