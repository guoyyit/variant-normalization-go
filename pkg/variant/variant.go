package variant

import (
	"fmt"
	"strings"
)

// SourceRef identifies the origin of a reported record.
type SourceRef struct {
	Source string
	File   string
	Line   int64
}

// SampleCall is one sample's genotype and per-sample fields.
// GT alleles are 0-based indices into a record's Alt list;
// -1 denotes a missing haplotype.
type SampleCall struct {
	Sample string
	Phased bool
	GT     []int
	Fields map[string]string
}

// Record is a VCF-style variant record using 1-based closed coordinates.
type Record struct {
	Chrom    string
	Pos      int
	ID       []string
	Ref      string
	Alt      []string
	HasQual  bool
	Qual     float64
	Filter   []string
	Info     map[string]string
	Samples  []SampleCall
	Sources  []SourceRef
	ConfigID string
}

// Key is the identity under which records are merged after normalization:
// chromosome, 1-based start position and reference allele.
func (r Record) Key() string {
	return fmt.Sprintf("%s\x00%d\x00%s", r.Chrom, r.Pos, r.Ref)
}

// Clone returns a deep copy of the record.
func (r Record) Clone() Record {
	out := r
	out.ID = cloneStrings(r.ID)
	out.Alt = cloneStrings(r.Alt)
	out.Filter = cloneStrings(r.Filter)
	out.Sources = append([]SourceRef(nil), r.Sources...)
	out.Info = cloneInfo(r.Info)
	out.Samples = make([]SampleCall, len(r.Samples))
	for i, sm := range r.Samples {
		cp := sm
		cp.GT = append([]int(nil), sm.GT...)
		cp.Fields = cloneInfo(sm.Fields)
		out.Samples[i] = cp
	}
	return out
}

// IsSymbolic reports whether alt is a symbolic allele such as <DEL>.
func IsSymbolic(allele string) bool {
	return strings.HasPrefix(allele, "<") && strings.HasSuffix(allele, ">")
}

// IsBreakend reports whether alt is a breakend allele.
func IsBreakend(allele string) bool {
	return strings.ContainsAny(allele, "[]")
}

// IsNoAlt reports the missing-alt placeholder "*".
func IsNoAlt(allele string) bool { return allele == "*" }

func cloneStrings(in []string) []string {
	if in == nil {
		return nil
	}
	return append([]string(nil), in...)
}

func cloneInfo(in map[string]string) map[string]string {
	if in == nil {
		return nil
	}
	out := make(map[string]string, len(in))
	for k, v := range in {
		out[k] = v
	}
	return out
}
