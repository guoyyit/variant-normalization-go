package vcfio

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"testing"

	"github.com/guoyyit/variant-normalization-go/variant"
)

const sampleVCF = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	rs1	A	T,C	50	PASS	AF=0.25,0.1;DP=100	GT:AD:PL	0/1:10,5,2:0,30,300	1|2:10,2,8:0,30,300
chr1	200	.	G	A	.	LowQual	.	GT	./1
`

func TestReader_Plain(t *testing.T) {
	r, err := NewReader(bytes.NewReader([]byte(sampleVCF)))
	if err != nil {
		t.Fatal(err)
	}
	if got := r.Header().FormatNumber["AD"]; got != "R" {
		t.Fatalf("AD Number = %q, want R", got)
	}
	if len(r.Samples()) != 2 || r.Samples()[0] != "S1" {
		t.Fatalf("samples = %v", r.Samples())
	}

	rec, err := r.Read()
	if err != nil {
		t.Fatal(err)
	}
	if rec.Chrom != "chr1" || rec.Pos != 100 || len(rec.Alt) != 2 {
		t.Fatalf("bad record: %+v", rec)
	}
	if rec.ID[0] != "rs1" || rec.Info["AF"] != "0.25,0.1" {
		t.Fatalf("id/info parse failed: %+v", rec)
	}
	if !rec.HasQual || rec.Qual != 50 {
		t.Fatalf("qual = %v has=%v", rec.Qual, rec.HasQual)
	}
	if gt := rec.Samples[1].GT; !gt.Phased || len(gt.Alleles) != 2 || gt.Alleles[1] != 2 {
		t.Fatalf("S2 GT = %+v", gt)
	}

	rec2, err := r.Read()
	if err != nil {
		t.Fatal(err)
	}
	if rec2.HasQual || rec2.Filter[0] != "LowQual" {
		t.Fatalf("qual/filter: %+v", rec2)
	}
	if gt := rec2.Samples[0].GT; gt.Alleles[0] != -1 || gt.Alleles[1] != 1 {
		t.Fatalf("missing GT = %+v", gt)
	}
	if _, err := r.Read(); !errors.Is(err, io.EOF) {
		t.Fatalf("want EOF, got %v", err)
	}
}

func TestReader_Gzip(t *testing.T) {
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	_, _ = gw.Write([]byte(sampleVCF))
	_ = gw.Close()

	r, err := NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	rec, err := r.Read()
	if err != nil || rec.Pos != 100 {
		t.Fatalf("gzip read: rec=%+v err=%v", rec, err)
	}
}

func TestReader_Invalid(t *testing.T) {
	r, _ := NewReader(bytes.NewReader([]byte("chr1\tx\t.\tA\tT\t.\t.\t.\n")))
	_, err := r.Read()
	var ve *variant.ValidationError
	if !errors.As(err, &ve) || ve.Kind != variant.KindInvalidFormat {
		t.Fatalf("want invalid_format, got %v", err)
	}
}
