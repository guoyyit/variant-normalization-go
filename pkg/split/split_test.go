package split

import (
	"reflect"
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func TestSplit_RemapGT(t *testing.T) {
	rec := variant.Record{
		Chrom: "1", Pos: 100, Ref: "A", Alt: []string{"T", "G"},
		Samples: []variant.SampleCall{
			{Sample: "s1", Phased: true, GT: []int{1, 2}},
			{Sample: "s2", GT: []int{0, 1}, Fields: map[string]string{"AD": "10,5,2", "DP": "17"}},
			{Sample: "s3", GT: []int{-1, 0}},
		},
	}
	parts := Split(rec)
	if len(parts) != 2 {
		t.Fatalf("want 2 splits, got %d", len(parts))
	}
	if !reflect.DeepEqual(parts[0].Samples[0].GT, []int{1, -1}) {
		t.Fatalf("s1 split0 GT = %v", parts[0].Samples[0].GT)
	}
	if !reflect.DeepEqual(parts[1].Samples[0].GT, []int{-1, 1}) {
		t.Fatalf("s1 split1 GT = %v", parts[1].Samples[0].GT)
	}
	if !reflect.DeepEqual(parts[0].Samples[2].GT, []int{-1, 0}) {
		t.Fatalf("s3 split0 GT = %v", parts[0].Samples[2].GT)
	}
	if got := parts[0].Samples[1].Fields["AD"]; got != "10,5" {
		t.Fatalf("AD remap split0 = %q", got)
	}
	if got := parts[1].Samples[1].Fields["AD"]; got != "10,2" {
		t.Fatalf("AD remap split1 = %q", got)
	}
	if got := parts[0].Samples[1].Fields["DP"]; got != "17" {
		t.Fatalf("DP should be preserved, got %q", got)
	}
	if parts[0].Samples[0].Phased != true {
		t.Fatal("phasing must be preserved")
	}

	// Information-preservation: every non-missing input allele maps to
	// exactly one split carrying it.
	totalAlt := 0
	for _, p := range parts {
		for _, sm := range p.Samples {
			totalAlt += CountGT(sm.GT)
		}
	}
	if totalAlt != 3 { // s1:1/2 and s2:0/1 -> three alt alleles
		t.Fatalf("alt allele count changed: %d", totalAlt)
	}
}

func TestSplit_SingleAltPassthrough(t *testing.T) {
	rec := variant.Record{Chrom: "1", Pos: 1, Ref: "A", Alt: []string{"T"}}
	if len(Split(rec)) != 1 {
		t.Fatal("single-alt record should yield one split")
	}
}
