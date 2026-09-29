package merge

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/config"
	"github.com/guoyyit/variant-normalization-go/variant"
)

func canon(source string, line int64, chrom string, pos int, ref string, alt []string,
	qual float64, hasQual bool, filter []string, samples []variant.Sample, provID []string) *variant.Canonical {
	return &variant.Canonical{
		Record: variant.Record{
			Chrom: chrom, Pos: pos, Ref: ref, Alt: alt,
			Qual: qual, HasQual: hasQual, Filter: filter,
			Format:  []string{"GT"},
			Samples: samples,
			Info:    map[string]string{},
			Source:  variant.Source{ID: source, Line: line},
		},
		Sources: []variant.Provenance{{SourceID: source, Line: line, ID: provID}},
	}
}

func gtSample(name, gt string, alleles ...int) variant.Sample {
	return variant.Sample{
		Name:   name,
		GT:     variant.Genotype{Alleles: alleles},
		Fields: map[string]string{"GT": gt},
	}
}

func TestMerge_DedupSameAlt(t *testing.T) {
	mg := NewMerger(config.Default())
	a := canon("srcA", 1, "c", 10, "A", []string{"T"}, 30, true, []string{"PASS"},
		[]variant.Sample{gtSample("s1", "0/1", 0, 1)}, []string{"rs1"})
	b := canon("srcB", 7, "c", 10, "A", []string{"T"}, 50, true, []string{"VQSR"},
		[]variant.Sample{gtSample("s2", "1/1", 1, 1)}, []string{"rs2"})
	out := mg.Merge([]*variant.Canonical{a, b})
	if len(out.Record.Alt) != 1 || out.Record.Alt[0] != "T" {
		t.Fatalf("alt = %v", out.Record.Alt)
	}
	if len(out.Sources) != 2 {
		t.Fatalf("sources = %d", len(out.Sources))
	}
	if !out.Record.HasQual || out.Record.Qual != 50 {
		t.Fatalf("qual = %v", out.Record.Qual)
	}
	if out.Record.Info["AC"] != "3" || out.Record.Info["AN"] != "4" || out.Record.Info["NS"] != "2" {
		t.Fatalf("stats = %+v", out.Record.Info)
	}
	if out.Record.Info["AF"] != "0.750000" {
		t.Fatalf("af = %s", out.Record.Info["AF"])
	}
	if len(out.Record.Filter) != 1 || out.Record.Filter[0] != "VQSR" {
		t.Fatalf("filter = %v", out.Record.Filter)
	}
}

func TestMerge_DifferentAlts(t *testing.T) {
	mg := NewMerger(config.Default())
	a := canon("a", 1, "c", 10, "A", []string{"T"}, 0, false, []string{"PASS"},
		[]variant.Sample{gtSample("s1", "0/1", 0, 1)}, nil)
	b := canon("b", 2, "c", 10, "A", []string{"C"}, 0, false, []string{"PASS"},
		[]variant.Sample{gtSample("s2", "1/1", 1, 1)}, nil)
	out := mg.Merge([]*variant.Canonical{a, b})
	if len(out.Record.Alt) != 2 || out.Record.Alt[0] != "C" || out.Record.Alt[1] != "T" {
		t.Fatalf("alts not sorted: %v", out.Record.Alt)
	}
	// s2 的 ALT C 在合并后为索引1，GT 1/1 -> 1/1；s1 的 T 为索引2，GT 0/1 -> 0/2
	byName := map[string]variant.Sample{}
	for _, s := range out.Record.Samples {
		byName[s.Name] = s
	}
	if g := byName["s1"].Fields["GT"]; g != "0/2" {
		t.Fatalf("s1 GT = %s", g)
	}
	if g := byName["s2"].Fields["GT"]; g != "1/1" {
		t.Fatalf("s2 GT = %s", g)
	}
	if out.Record.Info["AC"] != "2,1" {
		t.Fatalf("AC = %s", out.Record.Info["AC"])
	}
}

func TestMerge_OrderIndependent(t *testing.T) {
	mg := NewMerger(config.Default())
	mk := func() []*variant.Canonical {
		return []*variant.Canonical{
			canon("z", 3, "c", 10, "A", []string{"T"}, 9, true, nil,
				[]variant.Sample{gtSample("z", "1/1", 1, 1)}, nil),
			canon("a", 1, "c", 10, "A", []string{"C"}, 1, true, nil,
				[]variant.Sample{gtSample("a", "0/1", 0, 1)}, nil),
			canon("m", 2, "c", 10, "A", []string{"G"}, 5, true, nil,
				[]variant.Sample{gtSample("m", "0/0", 0, 0)}, nil),
		}
	}
	o1 := mg.Merge(mk())
	o2 := mg.Merge(reverse(mk()))
	if o1.Record.Key() != o2.Record.Key() {
		t.Fatalf("keys differ: %s vs %s", o1.Record.Key(), o2.Record.Key())
	}
	if !equalSources(o1.Sources, o2.Sources) {
		t.Fatalf("sources order differs: %v vs %v", o1.Sources, o2.Sources)
	}
	// Provenance 必须按来源排序
	if o1.Sources[0].SourceID != "a" || o1.Sources[2].SourceID != "z" {
		t.Fatalf("provenance not sorted: %v", o1.Sources)
	}
}

func TestMerge_PriorityBase(t *testing.T) {
	cfg := &config.Config{
		SourcePriority:  map[string]int{"gold": 0},
		DefaultPriority: 100,
		ConfigVersion:   "v2",
	}
	mg := NewMerger(cfg)
	low := canon("shallow", 1, "c", 10, "A", []string{"T"}, 99, true, nil, nil, nil)
	gold := canon("gold", 1, "c", 10, "A", []string{"T"}, 10, true, nil, nil, []string{"goldID"})
	out := mg.Merge([]*variant.Canonical{low, gold})
	if len(out.Record.ID) != 1 || out.Record.ID[0] != "goldID" {
		t.Fatalf("base id = %v", out.Record.ID)
	}
	if out.ConfigVersion != "v2" {
		t.Fatalf("config version = %s", out.ConfigVersion)
	}
}

func TestMerge_MinQual(t *testing.T) {
	cfg := &config.Config{DefaultPriority: 100, MinQual: 30, ConfigVersion: "x"}
	mg := NewMerger(cfg)
	a := canon("s", 1, "c", 10, "A", []string{"T"}, 10, true, nil, nil, nil)
	out := mg.Merge([]*variant.Canonical{a})
	if out.Record.Filter[0] != "LowQual" {
		t.Fatalf("filter = %v", out.Record.Filter)
	}
}

func reverse(in []*variant.Canonical) []*variant.Canonical {
	out := append([]*variant.Canonical(nil), in...)
	for i, j := 0, len(out)-1; i < j; i, j = i+1, j-1 {
		out[i], out[j] = out[j], out[i]
	}
	return out
}

func equalSources(a, b []variant.Provenance) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].SourceID != b[i].SourceID || a[i].Line != b[i].Line {
			return false
		}
	}
	return true
}
