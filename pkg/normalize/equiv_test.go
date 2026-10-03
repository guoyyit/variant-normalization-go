package normalize

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// These acceptance tests are independent of the Normalize implementation:
// ground truth is derived from first principles only.
//
//   - reconstructHaplotype applies a (POS, REF, ALT) record to the full
//     reference string. Normalization is semantically correct iff the input
//     record and its normalized output yield byte-identical haplotypes.
//   - equivalentSpellings independently re-anchors the record's window over
//     shared flanking context. Every returned spelling is verified by
//     reconstruction, not taken from Normalize. The canonical key is the
//     independently smallest (POS, REF, ALT) in the equivalence class.

func reconstructHaplotype(seq, ref, alt string, pos int) (string, bool) {
	end := pos - 1 + len(ref)
	if pos < 1 || end > len(seq) || seq[pos-1:end] != ref {
		return "", false
	}
	return seq[:pos-1] + alt + seq[end:], true
}

type spelling struct {
	pos      int
	ref, alt string
}

// equivalentSpellings enumerates every valid VCF spelling of the variant
// represented by (pos, ref, alt). A spelling is a choice of common prefix
// length L and suffix length R shared verbatim between the reference and the
// independently reconstructed target haplotype; the middle slices are the
// REF/ALT windows. Every candidate is independently verified by
// reconstruction, and only fully trimmed (canonical-window) spellings are
// kept: a shared end base is removable only while both windows would stay
// non-empty. Enumerating all L/R alignments also covers the multiple
// alignments of indels inside repeat runs, not just one fixed alignment.
func equivalentSpellings(seq, ref, alt string, pos int) []spelling {
	target, ok := reconstructHaplotype(seq, ref, alt, pos)
	if !ok {
		return nil
	}

	seen := map[spelling]bool{}
	maxLR := len(seq)
	if len(target) < maxLR {
		maxLR = len(target)
	}
	for l := 0; l <= maxLR; l++ {
		if seq[:l] != target[:l] {
			break
		}
		for r := 0; l+r <= len(seq) && l+r <= len(target); r++ {
			if seq[len(seq)-r:] != target[len(target)-r:] {
				break
			}
			rw := seq[l : len(seq)-r]
			aw := target[l : len(target)-r]
			if len(rw) < 1 || len(aw) < 1 || rw == aw {
				continue
			}
			if len(rw) > 1 && len(aw) > 1 && rw[0] == aw[0] {
				continue
			}
			if len(rw) > 1 && len(aw) > 1 && rw[len(rw)-1] == aw[len(aw)-1] {
				continue
			}
			p := l + 1
			got, ok := reconstructHaplotype(seq, rw, aw, p)
			if !ok || got != target {
				continue
			}
			seen[spelling{p, rw, aw}] = true
		}
	}

	out := make([]spelling, 0, len(seen))
	for s := range seen {
		out = append(out, s)
	}
	return out
}

func allAllelesOver(alpha []byte, maxLen int) []string {
	var out []string
	var rec func(prefix []byte, n int)
	rec = func(prefix []byte, n int) {
		if n == 0 {
			out = append(out, string(prefix))
			return
		}
		for _, b := range alpha {
			rec(append(prefix, b), n-1)
		}
	}
	for n := 1; n <= maxLen; n++ {
		rec(nil, n)
	}
	return out
}

// TestNormalize_HaplotypeEquivalence_Enumerated enumerates a bounded space
// of binary reference strings and allele pairs. For every legal input it
// independently reconstructs both haplotypes, checks idempotence, enumerates
// every equivalent spelling independently, and requires convergence to the
// independently chosen minimal key.
func TestNormalize_HaplotypeEquivalence_Enumerated(t *testing.T) {
	const refLen = 9
	const maxAllele = 3

	alleles := allAllelesOver([]byte{'A', 'C'}, maxAllele)
	count := 0
	for mask := 0; mask < 1<<refLen; mask++ {
		seq := make([]byte, refLen)
		for i := range seq {
			if mask&(1<<i) == 0 {
				seq[i] = 'A'
			} else {
				seq[i] = 'C'
			}
		}
		reference := string(seq)
		genomeRef := memRef{seq: reference, length: refLen}

		for _, r := range alleles {
			for pos := 1; pos+len(r)-1 <= refLen; pos++ {
				if reference[pos-1:pos-1+len(r)] != r {
					continue
				}
				for _, a := range alleles {
					if a == r {
						continue
					}
					count++
					rec := variant.Record{Chrom: "c", Pos: pos, Ref: r, Alt: []string{a}}
					out, err := Normalize(rec, genomeRef)
					if err != nil {
						t.Errorf("legal input %q %d:%s>%s rejected: %v", reference, pos, r, a, err)
						continue
					}

					// 1. Semantic equivalence via independent reconstruction.
					hIn, okIn := reconstructHaplotype(reference, r, a, pos)
					hOut, okOut := reconstructHaplotype(reference, out.Ref, out.Alt[0], out.Pos)
					if !okIn || !okOut || hIn != hOut {
						t.Fatalf("haplotype changed: ref=%q in=%d:%s>%s out=%d:%s>%s hIn=%q hOut=%q",
							reference, pos, r, a, out.Pos, out.Ref, out.Alt[0], hIn, hOut)
					}

					// 2. Idempotence.
					out2, err := Normalize(out, genomeRef)
					if err != nil {
						t.Fatalf("re-normalize error: %v", err)
					}
					if out2.Pos != out.Pos || out2.Ref != out.Ref || out2.Alt[0] != out.Alt[0] {
						t.Fatalf("not idempotent: %d:%s>%s then %d:%s>%s",
							out.Pos, out.Ref, out.Alt[0], out2.Pos, out2.Ref, out2.Alt[0])
					}

					// 3. Convergence over the independently enumerated
					// equivalence class, and equality with the independently
					// chosen minimal spelling.
					reps := equivalentSpellings(reference, r, a, pos)
					if len(reps) == 0 {
						t.Fatalf("no representations for ref=%q %d:%s>%s", reference, pos, r, a)
					}
					min := reps[0]
					for _, sp := range reps[1:] {
						if sp.pos < min.pos ||
							(sp.pos == min.pos && (sp.ref < min.ref ||
								(sp.ref == min.ref && sp.alt < min.alt))) {
							min = sp
						}
						got, err := Normalize(
							variant.Record{Chrom: "c", Pos: sp.pos, Ref: sp.ref, Alt: []string{sp.alt}},
							genomeRef)
						if err != nil {
							t.Fatalf("equivalent spelling %d:%s>%s rejected: %v", sp.pos, sp.ref, sp.alt, err)
						}
						if got.Pos != out.Pos || got.Ref != out.Ref || got.Alt[0] != out.Alt[0] {
							t.Fatalf("equivalent spellings diverge for ref=%q: %d:%s>%s -> %d:%s>%s vs %d:%s>%s",
								reference, sp.pos, sp.ref, sp.alt,
								got.Pos, got.Ref, got.Alt[0],
								out.Pos, out.Ref, out.Alt[0])
						}
					}
					if out.Pos != min.pos || out.Ref != min.ref || out.Alt[0] != min.alt {
						t.Fatalf("canonical key not minimal for ref=%q %d:%s>%s: got %d:%s>%s want %d:%s>%s",
							reference, pos, r, a,
							out.Pos, out.Ref, out.Alt[0], min.pos, min.ref, min.alt)
					}
				}
			}
		}
	}
	t.Logf("enumerated %d legal (reference,pos,refAllele,altAllele) combinations", count)
}
