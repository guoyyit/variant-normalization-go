package merge

import (
	"sort"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// Merge combines two records sharing chrom/pos/ref into one deterministic
// record. The result is independent of argument order: ALT slots are sorted,
// scalar conflicts are resolved by source priority, and every contributing
// source is retained.
func Merge(a, b variant.Record, cfg config.Config) variant.Record {
	if len(a.Alt) == 0 {
		return b.Clone()
	}
	if len(b.Alt) == 0 {
		return a.Clone()
	}

	out := a.Clone()

	altOrder := altUnion(a.Alt, b.Alt)
	out.Alt = make([]string, len(altOrder))
	for allele, slot := range altOrder {
		out.Alt[slot-1] = allele
	}

	out.ID = mergeIDs(a.ID, b.ID)
	out.Sources = mergeSources(a.Sources, b.Sources)
	out.Samples = mergeSamples(a, b, altOrder, cfg)
	out.Filter = mergeFilter(a.Filter, b.Filter, cfg)
	out.Info = mergeInfo(a.Info, b.Info)
	out.HasQual, out.Qual = mergeQual(a, b, cfg)
	out.ConfigID = cfg.ID
	return out
}

// altUnion builds the canonical 1-based slot map for the union of both ALT
// lists. Slots are assigned lexicographically so the merge is symmetric.
func altUnion(a, b []string) map[string]int {
	seen := map[string]bool{}
	var union []string
	for _, allele := range append(append([]string{}, a...), b...) {
		if !seen[allele] {
			seen[allele] = true
			union = append(union, allele)
		}
	}
	sort.Strings(union)
	idx := make(map[string]int, len(union))
	for slot, allele := range union {
		idx[allele] = slot + 1
	}
	return idx
}

func mergeIDs(a, b []string) []string {
	seen := map[string]bool{}
	var out []string
	for _, id := range append(append([]string{}, a...), b...) {
		if id == "" || id == "." || seen[id] {
			continue
		}
		seen[id] = true
		out = append(out, id)
	}
	sort.Strings(out)
	return out
}

// mergeSources keeps every distinct (source,file) provenance pair,
// deterministically ordered.
func mergeSources(a, b []variant.SourceRef) []variant.SourceRef {
	seen := map[[2]string]bool{}
	out := append([]variant.SourceRef{}, a...)
	for _, src := range a {
		seen[[2]string{src.Source, src.File}] = true
	}
	for _, src := range b {
		key := [2]string{src.Source, src.File}
		if seen[key] {
			continue
		}
		seen[key] = true
		out = append(out, src)
	}
	sort.SliceStable(out, func(i, j int) bool {
		if out[i].Source != out[j].Source {
			return out[i].Source < out[j].Source
		}
		return out[i].File < out[j].File
	})
	return out
}

// recordPriority derives a representative source priority for scalar
// arbitration: the highest priority among contributing sources.
func recordPriority(rec variant.Record, cfg config.Config) int {
	if len(rec.Sources) == 0 {
		return cfg.DefaultPriority
	}
	best := cfg.Priority(rec.Sources[0].Source)
	for _, src := range rec.Sources[1:] {
		if p := cfg.Priority(src.Source); p > best {
			best = p
		}
	}
	return best
}

// mergeSamples combines genotype calls keyed by sample name. The
// lower-priority record is ingested first, so conflicting called alleles are
// overridden by the higher-priority call; ties favor the second ingest, but
// callers feed equal-priority records in a fixed order.
func mergeSamples(a, b variant.Record, altOrder map[string]int, cfg config.Config) []variant.SampleCall {
	byName := map[string]variant.SampleCall{}
	var order []string

	ingest := func(rec variant.Record) {
		for _, sm := range rec.Samples {
			call := sm
			if len(rec.Alt) == 1 {
				call.GT = remapGTTo(sm.GT, altOrder[rec.Alt[0]])
			}
			call.Fields = cloneFields(sm.Fields)
			existing, ok := byName[sm.Sample]
			if !ok {
				byName[sm.Sample] = call
				order = append(order, sm.Sample)
				continue
			}
			byName[sm.Sample] = reconcileSample(existing, call)
		}
	}

	if recordPriority(a, cfg) <= recordPriority(b, cfg) {
		ingest(a)
		ingest(b)
	} else {
		ingest(b)
		ingest(a)
	}

	sort.Strings(order)
	out := make([]variant.SampleCall, 0, len(order))
	for _, name := range order {
		out = append(out, byName[name])
	}
	return out
}

// reconcileSample merges two calls of the same sample. A called allele
// replaces a missing one; conflicting non-missing alleles keep the incoming
// value (the higher-priority record is always ingested second).
func reconcileSample(cur, incoming variant.SampleCall) variant.SampleCall {
	n := len(cur.GT)
	if len(incoming.GT) > n {
		n = len(incoming.GT)
	}
	merged := make([]int, n)
	for i := range merged {
		merged[i] = -1
	}
	copy(merged, cur.GT)
	for i := 0; i < len(incoming.GT); i++ {
		inv := incoming.GT[i]
		curV := -1
		if i < len(cur.GT) {
			curV = cur.GT[i]
		}
		if inv != -1 && (curV == -1 || curV != inv) {
			merged[i] = inv
		}
	}
	cur.GT = merged
	if len(cur.Fields) == 0 {
		cur.Fields = map[string]string{}
	}
	for key, value := range incoming.Fields {
		if _, exists := cur.Fields[key]; !exists {
			cur.Fields[key] = value
		}
	}
	if incoming.Phased {
		cur.Phased = true
	}
	return cur
}

// remapGTTo renumbers local single-alt genotype indices (0/1) into the output
// multi-alt space; local alt 1 maps to slot.
func remapGTTo(gt []int, slot int) []int {
	out := make([]int, len(gt))
	for i, allele := range gt {
		if allele == 1 {
			out[i] = slot
		} else {
			out[i] = allele
		}
	}
	return out
}

func cloneFields(in map[string]string) map[string]string {
	if len(in) == 0 {
		return nil
	}
	out := make(map[string]string, len(in))
	for k, v := range in {
		out[k] = v
	}
	return out
}

// mergeQual aggregates quality according to cfg.QualAgg. "weighted" uses
// INFO/DP from each side as weight when available, otherwise max.
func mergeQual(a, b variant.Record, cfg config.Config) (bool, float64) {
	if cfg.QualAgg == "min" && a.HasQual && b.HasQual {
		if a.Qual < b.Qual {
			return true, a.Qual
		}
		return true, b.Qual
	}
	if cfg.QualAgg == "weighted" && a.HasQual && b.HasQual {
		dpa := infoInt(a.Info, "DP")
		dpb := infoInt(b.Info, "DP")
		if dpa > 0 || dpb > 0 {
			return true, (a.Qual*dpa + b.Qual*dpb) / (dpa + dpb)
		}
	}
	switch {
	case a.HasQual && b.HasQual:
		if b.Qual > a.Qual {
			return true, b.Qual
		}
		return true, a.Qual
	case a.HasQual:
		return true, a.Qual
	default:
		return b.HasQual, b.Qual
	}
}

func infoInt(info map[string]string, key string) float64 {
	if info == nil {
		return 0
	}
	v, ok := info[key]
	if !ok || v == "" {
		return 0
	}
	n := 0
	for _, r := range v {
		if r < '0' || r > '9' {
			return 0
		}
		n = n*10 + int(r-'0')
	}
	return float64(n)
}

// mergeFilter:
//   - mode "all" (default): PASS only when every record PASSes; concrete
//     filter codes are unioned, deduped and sorted;
//   - mode "any": PASS if any record PASSes.
func mergeFilter(a, b []string, cfg config.Config) []string {
	passA := isPass(a)
	passB := isPass(b)
	if cfg.FilterMode == "any" && (passA || passB) {
		return []string{"PASS"}
	}
	seen := map[string]bool{}
	var out []string
	collect := func(filters []string) {
		for _, f := range filters {
			if f == "PASS" || f == "." || f == "" || seen[f] {
				continue
			}
			seen[f] = true
			out = append(out, f)
		}
	}
	collect(a)
	collect(b)
	sort.Strings(out)
	if len(out) == 0 {
		return []string{"PASS"}
	}
	return out
}

func isPass(filters []string) bool {
	if len(filters) == 0 {
		return true
	}
	for _, f := range filters {
		if f != "PASS" {
			return false
		}
	}
	return true
}

// mergeInfo keeps allele-agnostic INFO keys with identical values and keys
// provided by only one side. AC/AN/AF are recomputed later; conflicting
// scalar values are dropped rather than silently chosen.
func mergeInfo(a, b map[string]string) map[string]string {
	reserved := map[string]bool{"AC": true, "AN": true, "AF": true}
	keys := map[string]bool{}
	for k := range a {
		keys[k] = true
	}
	for k := range b {
		keys[k] = true
	}
	out := map[string]string{}
	for k := range keys {
		if reserved[k] {
			continue
		}
		va, oka := a[k]
		vb, okb := b[k]
		switch {
		case oka && okb && va == vb:
			out[k] = va
		case oka && !okb:
			out[k] = va
		case okb && !oka:
			out[k] = vb
		}
	}
	if len(out) == 0 {
		return nil
	}
	return out
}

// SameKey reports whether two records merge into the same canonical group.
func SameKey(a, b variant.Record) bool {
	return a.Chrom == b.Chrom && a.Pos == b.Pos && a.Ref == b.Ref
}
