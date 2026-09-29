package normalize

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func norm(t *testing.T, seq string, pos int, ref, alt string) variant.Record {
	t.Helper()
	r := memRef{seq: seq, length: len(seq)}
	out, err := Normalize(variant.Record{Chrom: "1", Pos: pos, Ref: ref, Alt: []string{alt}}, r)
	if err != nil {
		t.Fatalf("normalize %d:%s>%s: %v", pos, ref, alt, err)
	}
	return out
}

func TestCanonicalCases(t *testing.T) {
	// Reference:  C A A A A C T  (positions 1..7)
	const seq = "CAAAACT"

	cases := []struct {
		name     string
		pos      int
		ref, alt string
		wPos     int
		wRef     string
		wAlt     string
	}{
		// Deletion of one A from the A-run: leftmost is pos 1 CA>C.
		{"del_run", 3, "AAAC", "AAC", 1, "CA", "C"},
		// Insertion of A after the run: canonical anchor at leftmost.
		// Inserting A inside the A-run is the same variant regardless of
		// reported position; canonical anchor is pos 1: C>CA.
		{"ins_run", 5, "AC", "AAC", 1, "C", "CA"},
		// Simple SNP unchanged.
		{"snp", 7, "T", "G", 7, "T", "G"},
		// SNP with trailing matching base trimmed.
		{"snp_trim", 6, "CT", "GT", 6, "C", "G"},
		// Inserting C... T>CT at the tail still walks left through the
		// identical anchor to its canonical position.
		{"ins_tail", 7, "T", "CT", 1, "C", "CA"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := norm(t, seq, tc.pos, tc.ref, tc.alt)
			if got.Pos != tc.wPos || got.Ref != tc.wRef || got.Alt[0] != tc.wAlt {
				t.Fatalf("got %d:%s>%s want %d:%s>%s",
					got.Pos, got.Ref, got.Alt[0], tc.wPos, tc.wRef, tc.wAlt)
			}
		})
	}
}

func TestCanonical_CrossRepresentations(t *testing.T) {
	const seq = "CAAAACT"
	// Two different textual spellings of the same deletion must converge.
	a := norm(t, seq, 3, "AAAC", "AAC")
	b := norm(t, seq, 2, "AAC", "AC")
	if a.Pos != b.Pos || a.Ref != b.Ref || a.Alt[0] != b.Alt[0] {
		t.Fatalf("representations diverge: %d:%s>%s vs %d:%s>%s",
			a.Pos, a.Ref, a.Alt[0], b.Pos, b.Ref, b.Alt[0])
	}
}

func TestCanonical_SymbolicPassthrough(t *testing.T) {
	r := memRef{seq: "ACGT", length: 4}
	out, err := Normalize(variant.Record{Chrom: "1", Pos: 2, Ref: "C", Alt: []string{"<DEL>"}}, r)
	if err != nil {
		t.Fatal(err)
	}
	if out.Pos != 2 || out.Ref != "C" || out.Alt[0] != "<DEL>" {
		t.Fatalf("symbolic altered: %d:%s>%s", out.Pos, out.Ref, out.Alt[0])
	}
}
