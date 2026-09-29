package split

import (
	"strings"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// Split breaks a multi-allelic record into single-alt records,
// remapping each sample's genotype accordingly.
func Split(rec variant.Record) []variant.Record {
	if len(rec.Alt) <= 1 {
		return []variant.Record{rec.Clone()}
	}

	out := make([]variant.Record, 0, len(rec.Alt))
	for altIdx := range rec.Alt {
		split := rec.Clone()
		split.Alt = []string{rec.Alt[altIdx]}

		for i := range split.Samples {
			split.Samples[i].GT = remapGT(rec.Samples[i].GT, altIdx)
			split.Samples[i].Fields = remapFields(rec.Samples[i].Fields, altIdx, len(rec.Alt))
		}

		// Allele-counting INFO fields are recomputed after merging.
		delete(split.Info, "AC")
		delete(split.Info, "AN")
		delete(split.Info, "AF")

		out = append(out, split)
	}
	return out
}

// remapGT maps genotype indices: 0 (REF) stays 0, the chosen alt i+1 becomes
// 1, every other allele becomes -1 (missing). Missing (-1) stays missing.
func remapGT(gt []int, altIdx int) []int {
	out := make([]int, len(gt))
	for j, allele := range gt {
		switch {
		case allele == -1:
			out[j] = -1
		case allele == 0:
			out[j] = 0
		case allele == altIdx+1:
			out[j] = 1
		default:
			out[j] = -1
		}
	}
	return out
}

// remapFields performs best-effort per-allele FORMAT remapping.
//
// Without the VCF header the Number type of each field is unknown; the
// convention used here is:
//   - AD-like fields (one value per allele incl. REF) keep the REF and
//     the selected alt values ("R" cardinality);
//   - other comma-separated fields keep the selected alt value
//     ("A" cardinality).
func remapFields(fields map[string]string, altIdx, altCount int) map[string]string {
	if len(fields) == 0 {
		return nil
	}
	out := make(map[string]string, len(fields))
	for key, value := range fields {
		if !strings.Contains(value, ",") {
			out[key] = value
			continue
		}
		parts := strings.Split(value, ",")
		switch {
		case strings.EqualFold(key, "AD") && len(parts) == altCount+1:
			out[key] = parts[0] + "," + parts[altIdx+1]
		case altIdx+1 < len(parts):
			out[key] = parts[altIdx+1]
		default:
			out[key] = value
		}
	}
	return out
}

// CountGT returns the number of times alt index 1 occurs in a remapped
// genotype, ignoring missing calls.
func CountGT(gt []int) int {
	count := 0
	for _, allele := range gt {
		if allele == 1 {
			count++
		}
	}
	return count
}

// CalledChromosomes returns the number of called (non-missing) chromosomes.
func CalledChromosomes(gt []int) int {
	count := 0
	for _, allele := range gt {
		if allele >= 0 {
			count++
		}
	}
	return count
}
