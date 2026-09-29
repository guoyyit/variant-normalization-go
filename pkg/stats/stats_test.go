package stats

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func TestRecompute_AC_AN_AF(t *testing.T) {
	cfg := config.Default()
	rec := variant.Record{
		Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T", "G"},
		Samples: []variant.SampleCall{
			{Sample: "s1", GT: []int{1, 2}},
			{Sample: "s2", GT: []int{0, 1}, Fields: map[string]string{"DP": "10"}},
			{Sample: "s3", GT: []int{-1, 0}, Fields: map[string]string{"DP": "20"}},
		},
	}
	Recompute(&rec, cfg)
	if rec.Info["AC"] != "2,1" {
		t.Fatalf("AC = %q", rec.Info["AC"])
	}
	if rec.Info["AN"] != "5" {
		t.Fatalf("AN = %q", rec.Info["AN"])
	}
	// 2/5 = 0.4, 1/5 = 0.2
	if rec.Info["AF"] != "0.400000,0.200000" {
		t.Fatalf("AF = %q", rec.Info["AF"])
	}
	if rec.Info["DP"] != "30" {
		t.Fatalf("DP = %q", rec.Info["DP"])
	}
	if len(rec.Filter) != 1 || rec.Filter[0] != "PASS" {
		t.Fatalf("filter = %v", rec.Filter)
	}
}

func TestRecompute_Thresholds(t *testing.T) {
	cfg := config.Default()
	cfg.MinQual = 50
	cfg.MinAF = 0.5
	rec := variant.Record{
		Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"},
		HasQual: true, Qual: 20,
		Samples: []variant.SampleCall{{Sample: "s1", GT: []int{0, 1}}},
	}
	Recompute(&rec, cfg)
	// AF=0.5 passes min-af boundary; qual 20 < 50 fails.
	found := map[string]bool{}
	for _, f := range rec.Filter {
		found[f] = true
	}
	if !found["MinQual"] {
		t.Fatalf("missing MinQual in %v", rec.Filter)
	}

	rec2 := variant.Record{
		Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"},
		HasQual: true, Qual: 99,
		Samples: []variant.SampleCall{
			{Sample: "s1", GT: []int{1, 1}},
			{Sample: "s2", GT: []int{0, 0}},
		},
	}
	Recompute(&rec2, cfg)
	found = map[string]bool{}
	for _, f := range rec2.Filter {
		found[f] = true
	}
	if found["MinAF"] {
		t.Fatalf("unexpected MinAF in %v", rec2.Filter)
	}
	if !found["PASS"] {
		t.Fatalf("expected PASS in %v", rec2.Filter)
	}
}
