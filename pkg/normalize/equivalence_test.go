package normalize

import (
	"strings"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// applyHaplotype reconstructs the haplotype sequence by applying the variant
// to the reference. This is independent of the normalization implementation:
// it only splices the ALT allele into the reference at the reported window.
func applyHaplotype(seq string, rec variant.Record) string {
	var b strings.Builder
	b.WriteString(seq[:rec.Pos-1])
	b.WriteString(rec.Alt[0])
	b.WriteString(seq[rec.Pos-1+len(rec.Ref):])
	return b.String()
}

// stringsOver returns all strings of length n over the alphabet.
func stringsOver(alphabet string, n int) []string {
	if n == 0 {
		return []string{""}
	}
	short := stringsOver(alphabet, n-1)
	out := make([]string, 0, len(short)*len(alphabet))
	for _, s := range short {
		for _, c := range alphabet {
			out = append(out, s+string(c))
		}
	}
	return out
}

// enumerateSpellings generates every legal spelling for the reference:
// every position, every REF length (matching the reference), and every ALT
// string over the alphabet.
func enumerateSpellings(seq, alphabet string, maxRefLen, maxAltLen int) []variant.Record {
	var recs []variant.Record
	for pos := 1; pos <= len(seq); pos++ {
		for refLen := 1; refLen <= maxRefLen && pos+refLen-1 <= len(seq); refLen++ {
			ref := seq[pos-1 : pos-1+refLen]
			for altLen := 1; altLen <= maxAltLen; altLen++ {
				for _, alt := range stringsOver(alphabet, altLen) {
					if alt == ref {
						continue
					}
					recs = append(recs, variant.Record{
						Chrom: "1", Pos: pos, Ref: ref, Alt: []string{alt},
					})
				}
			}
		}
	}
	return recs
}

// TestEquivalenceEnumeration is the independent equivalence acceptance test.
// It enumerates all length-constrained allele pairs over a finite alphabet at
// every position, groups them by the haplotype they produce, and verifies for
// each spelling that:
//
//  1. Semantic equivalence: input and output produce the same haplotype.
//  2. Idempotence: normalizing the output again changes nothing.
//  3. Convergence: all spellings of the same haplotype converge to the same
//     canonical (CHROM, POS, REF, ALT).
//
// Expected values are never taken from the implementation; they come solely
// from haplotype reconstruction against the reference.
func TestEquivalenceEnumeration(t *testing.T) {
	refs := []struct {
		seq, alphabet string
	}{
		{"ACACACACACAC", "AC"},
		{"ACGTACGT", "ACGT"},
		{"ACGTNACGT", "ACGT"},
		{"CAAAACT", "ACGT"},
	}
	for _, tc := range refs {
		t.Run(tc.seq, func(t *testing.T) {
			r := memRef{seq: tc.seq, length: len(tc.seq)}
			spellings := enumerateSpellings(tc.seq, tc.alphabet, 3, 3)

			groups := map[string][]variant.Record{}
			for _, rec := range spellings {
				groups[applyHaplotype(tc.seq, rec)] =
					append(groups[applyHaplotype(tc.seq, rec)], rec)
			}

			for hap, group := range groups {
				var canonical variant.Record
				for i, rec := range group {
					out, err := Normalize(rec, r)
					if err != nil {
						t.Fatalf("normalize %d:%s>%s: %v",
							rec.Pos, rec.Ref, rec.Alt[0], err)
					}

					// 1. Semantic equivalence via haplotype reconstruction.
					hIn := applyHaplotype(tc.seq, rec)
					hOut := applyHaplotype(tc.seq, out)
					if hIn != hOut {
						t.Fatalf("haplotype mismatch for %d:%s>%s: in=%q out=%q",
							rec.Pos, rec.Ref, rec.Alt[0], hIn, hOut)
					}

					// 2. Idempotence.
					out2, err := Normalize(out, r)
					if err != nil {
						t.Fatalf("re-normalize %d:%s>%s: %v",
							out.Pos, out.Ref, out.Alt[0], err)
					}
					if out.Pos != out2.Pos || out.Ref != out2.Ref ||
						out.Alt[0] != out2.Alt[0] {
						t.Fatalf("not idempotent: %d:%s>%s -> %d:%s>%s",
							out.Pos, out.Ref, out.Alt[0],
							out2.Pos, out2.Ref, out2.Alt[0])
					}

					// 3. Convergence across all spellings of this haplotype.
					if i == 0 {
						canonical = out
					} else if out.Pos != canonical.Pos || out.Ref != canonical.Ref ||
						out.Alt[0] != canonical.Alt[0] {
						t.Fatalf("convergence failed for haplotype %q: "+
							"%d:%s>%s vs %d:%s>%s",
							hap, canonical.Pos, canonical.Ref, canonical.Alt[0],
							out.Pos, out.Ref, out.Alt[0])
					}
				}
			}
		})
	}
}
