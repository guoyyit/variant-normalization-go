package split

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/variant"
	"github.com/guoyyit/variant-normalization-go/vcfio"
)

func multiRecord() *variant.Record {
	return &variant.Record{
		Chrom:  "chr1",
		Pos:    100,
		Ref:    "A",
		Alt:    []string{"T", "C"},
		Filter: []string{"PASS"},
		Format: []string{"GT", "AD", "PL", "DP"},
		Samples: []variant.Sample{
			{Name: "S1", Fields: map[string]string{
				"GT": "0/1", "AD": "10,5,2", "PL": "0,30,300,40,400,4000", "DP": "17",
			}, GT: variant.Genotype{Alleles: []int{0, 1}}},
			{Name: "S2", Fields: map[string]string{
				"GT": "1|2", "AD": "10,5,2", "PL": "0,30,300,40,400,4000", "DP": "17",
			}, GT: variant.Genotype{Alleles: []int{1, 2}, Phased: true}},
		},
	}
}

func header() *vcfio.Header {
	return &vcfio.Header{FormatNumber: map[string]string{
		"AD": "R", "PL": "G", "DP": "1",
	}}
}

func TestSplit_Count(t *testing.T) {
	got, err := Split(multiRecord(), header())
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 2 {
		t.Fatalf("got %d children", len(got))
	}
	if got[0].Alt[0] != "T" || got[1].Alt[0] != "C" {
		t.Fatalf("alts = %s / %s", got[0].Alt, got[1].Alt)
	}
}

func TestSplit_Genotype(t *testing.T) {
	got, _ := Split(multiRecord(), header())
	// 第一条 ALT=T：S1 0/1 保持 0/1；S2 1|2 -> 1|. （2 不属于该拆分）
	if g := got[0].Samples[0].Fields["GT"]; g != "0/1" {
		t.Fatalf("S1 child0 GT = %s", g)
	}
	if g := got[0].Samples[1].Fields["GT"]; g != "1|." {
		t.Fatalf("S2 child0 GT = %s", g)
	}
	// 第二条 ALT=C：S1 0/1 -> 0/.；S2 1|2 -> .|1
	if g := got[1].Samples[0].Fields["GT"]; g != "0/." {
		t.Fatalf("S1 child1 GT = %s", g)
	}
	if g := got[1].Samples[1].Fields["GT"]; g != ".|1" {
		t.Fatalf("S2 child1 GT = %s", g)
	}
	// 语义一致：拆分后基因型结构里的等位数与相位不变
	if !got[1].Samples[1].GT.Phased {
		t.Fatalf("phasing lost")
	}
}

func TestSplit_Fields(t *testing.T) {
	got, _ := Split(multiRecord(), header())
	if ad := got[1].Samples[0].Fields["AD"]; ad != "10,2" {
		t.Fatalf("AD R-remap = %s", ad)
	}
	// PL(G)：双等位取 00,01,11；ALT=C(下标2) 在三等位G序中 0/2->下标3=40, 2/2->下标5=4000
	if pl := got[1].Samples[0].Fields["PL"]; pl != "0,40,4000" {
		t.Fatalf("PL G-remap = %s", pl)
	}
	if dp := got[0].Samples[0].Fields["DP"]; dp != "17" {
		t.Fatalf("DP scalar = %s", dp)
	}
}

func TestSplit_SingleAllele(t *testing.T) {
	rec := &variant.Record{Ref: "A", Alt: []string{"T"}}
	got, err := Split(rec, nil)
	if err != nil || len(got) != 1 || got[0] != rec {
		t.Fatalf("single allele should pass through: %v %v", got, err)
	}
}

func TestSplit_InfoPreserved(t *testing.T) {
	rec := multiRecord()
	rec.Info = map[string]string{"DP": "100"}
	got, _ := Split(rec, header())
	if got[0].Info["DP"] != "100" {
		t.Fatalf("INFO not preserved")
	}
}
