package merge

import (
	"reflect"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func src(name string) []variant.SourceRef {
	return []variant.SourceRef{{Source: name, File: name + ".vcf", Line: 1}}
}

func TestMerge_OrderIndependent(t *testing.T) {
	cfg := config.Default()
	a := variant.Record{Chrom: "1", Pos: 10, Ref: "A", Alt: []string{"T"},
		Filter: []string{"PASS"}, Sources: src("a")}
	b := variant.Record{Chrom: "1", Pos: 10, Ref: "A", Alt: []string{"G"},
		Filter: []string{"PASS"}, Sources: src("b")}

	m1 := Merge(a, b, cfg)
	m2 := Merge(b, a, cfg)
	if !reflect.DeepEqual(m1.Alt, m2.Alt) || !reflect.DeepEqual(m1.Alt, []string{"G", "T"}) {
		t.Fatalf("alt order not symmetric: %v vs %v", m1.Alt, m2.Alt)
	}
	if len(m1.Sources) != 2 || m1.Sources[0].Source != "a" || m1.Sources[1].Source != "b" {
		t.Fatalf("provenance wrong: %+v", m1.Sources)
	}
	if !reflect.DeepEqual(m1.Filter, []string{"PASS"}) {
		t.Fatalf("filter = %v", m1.Filter)
	}

	// Re-merge must be idempotent.
	m3 := Merge(m1, a, cfg)
	if len(m3.Sources) != 2 {
		t.Fatal("duplicate sources not deduped")
	}
}

func TestMerge_FilterModes(t *testing.T) {
	cfg := config.Default() // all PASS -> PASS; else codes
	a := variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"},
		Filter: []string{"PASS"}, Sources: src("a")}
	b := variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"},
		Filter: []string{"LowQual"}, Sources: src("b")}
	m := Merge(a, b, cfg)
	if !reflect.DeepEqual(m.Filter, []string{"LowQual"}) {
		t.Fatalf("strict mode filter = %v", m.Filter)
	}

	cfg.FilterMode = "any"
	m = Merge(a, b, cfg)
	if !reflect.DeepEqual(m.Filter, []string{"PASS"}) {
		t.Fatalf("lenient mode filter = %v", m.Filter)
	}
}

func TestMerge_SampleGT(t *testing.T) {
	cfg := config.Default()
	a := variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"},
		Sources: src("a"),
		Samples: []variant.SampleCall{{Sample: "s1", GT: []int{0, 1}}}}
	b := variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"G"},
		Sources: src("b"),
		Samples: []variant.SampleCall{{Sample: "s1", GT: []int{1, 1}}}}
	m := Merge(a, b, cfg)
	// sorted alts: G=1, T=2. b-side GT {1,1} -> {1,1}; a-side {0,1} -> {0,2};
	// equal priority: b ingested second wins conflicts -> {1,1}
	want := []int{1, 1}
	if !reflect.DeepEqual(m.Samples[0].GT, want) {
		t.Fatalf("GT = %v, want %v", m.Samples[0].GT, want)
	}

	// Priority makes the high-priority source win regardless of side.
	cfg.SourcePriority = map[string]int{"a": 10}
	m = Merge(a, b, cfg)
	want = []int{0, 2}
	if !reflect.DeepEqual(m.Samples[0].GT, want) {
		t.Fatalf("priority GT = %v, want %v", m.Samples[0].GT, want)
	}
	m = Merge(b, a, cfg)
	if !reflect.DeepEqual(m.Samples[0].GT, want) {
		t.Fatalf("priority GT not symmetric: %v", m.Samples[0].GT)
	}
}

func TestMerge_QualAgg(t *testing.T) {
	cfg := config.Default()
	a := variant.Record{Alt: []string{"T"}, HasQual: true, Qual: 30,
		Info: map[string]string{"DP": "10"}, Sources: src("a")}
	b := variant.Record{Alt: []string{"T"}, HasQual: true, Qual: 10,
		Info: map[string]string{"DP": "30"}, Sources: src("b")}
	m := Merge(a, b, cfg)
	if m.Qual != 30 {
		t.Fatalf("max qual = %v", m.Qual)
	}
	cfg.QualAgg = "weighted"
	m = Merge(a, b, cfg)
	// (30*10 + 10*30)/40 = 15
	if m.Qual != 15 {
		t.Fatalf("weighted qual = %v", m.Qual)
	}
}
