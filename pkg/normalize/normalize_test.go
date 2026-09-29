package normalize

import (
	"errors"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

type memRef struct {
	seq    string
	length int
}

func (m memRef) Substring(chrom string, start, end int) (string, error) {
	return m.seq[start-1 : end-1], nil
}

func (m memRef) Length(chrom string) (int, error) { return m.length, nil }

// seq: positions 1..8 = C A A A A A C T
var testRef = memRef{seq: "CAAAAACT", length: 8}

func TestNormalize_LeftAlignDeletion(t *testing.T) {
	ref := memRef{seq: "CAAAACT", length: 7} // C A A A A C T
	// Deletion of one A in the homopolymer run canonicalizes to the
	// leftmost representation: REF=CA ALT=C at pos 1.
	rec := variant.Record{Chrom: "chr1", Pos: 3, Ref: "AAAC", Alt: []string{"AAC"}}
	out, err := Normalize(rec, ref)
	if err != nil {
		t.Fatal(err)
	}
	if out.Pos != 1 || out.Ref != "CA" || out.Alt[0] != "C" {
		t.Fatalf("got pos=%d ref=%s alt=%s", out.Pos, out.Ref, out.Alt[0])
	}
}

func TestNormalize_TrimSNP(t *testing.T) {
	ref := memRef{seq: "ACGTACGT", length: 8}
	rec := variant.Record{Chrom: "1", Pos: 2, Ref: "CG", Alt: []string{"TG"}}
	out, err := Normalize(rec, ref)
	if err != nil {
		t.Fatal(err)
	}
	// anchor base must remain
	if out.Pos != 2 || out.Ref != "C" || out.Alt[0] != "T" {
		t.Fatalf("got pos=%d ref=%s alt=%s", out.Pos, out.Ref, out.Alt[0])
	}
}

func TestNormalize_Stable(t *testing.T) {
	ref := memRef{seq: "CAAAACT", length: 7}
	first, err := Normalize(variant.Record{Chrom: "c", Pos: 3, Ref: "AAAC", Alt: []string{"AAC"}}, ref)
	if err != nil {
		t.Fatal(err)
	}
	second, err := Normalize(first, ref)
	if err != nil {
		t.Fatal(err)
	}
	if first.Pos != second.Pos || first.Ref != second.Ref || first.Alt[0] != second.Alt[0] {
		t.Fatal("normalization is not idempotent")
	}
}

func TestNormalize_InvalidFormat(t *testing.T) {
	_, err := Normalize(variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"A"}}, testRef)
	if !errors.Is(err, ErrInvalidFormat) {
		t.Fatalf("want ErrInvalidFormat, got %v", err)
	}
	_, err = Normalize(variant.Record{Chrom: "1", Pos: 1, Ref: "X", Alt: []string{"A"}}, testRef)
	if !errors.Is(err, ErrInvalidFormat) {
		t.Fatalf("want ErrInvalidFormat, got %v", err)
	}
	_, err = Normalize(variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"A", "T"}}, testRef)
	if !errors.Is(err, ErrInvalidFormat) {
		t.Fatalf("want ErrInvalidFormat, got %v", err)
	}
}

func TestNormalize_OutOfBounds(t *testing.T) {
	_, err := Normalize(variant.Record{Chrom: "1", Pos: 8, Ref: "CT", Alt: []string{"C"}}, testRef)
	if !errors.Is(err, ErrOutOfBounds) {
		t.Fatalf("want ErrOutOfBounds, got %v", err)
	}
}
