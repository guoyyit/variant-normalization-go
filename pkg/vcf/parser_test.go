package vcf

import (
	"io"
	"reflect"
	"strings"
	"testing"
)

const sampleVCF = `##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	rs1	A	T,G	50	PASS	DP=100;NS=2	GT:DP	0|1:30	1/2:20
1	200	.	C	<DEL>	.	LowQual	.	GT	./.
`

func TestParser(t *testing.T) {
	p := NewParser(strings.NewReader(sampleVCF))

	rec, err := p.Next()
	if err != nil {
		t.Fatal(err)
	}
	if rec.Chrom != "1" || rec.Pos != 100 || !reflect.DeepEqual(rec.Alt, []string{"T", "G"}) {
		t.Fatalf("fixed fields wrong: %+v", rec)
	}
	if !rec.HasQual || rec.Qual != 50 {
		t.Fatal("qual wrong")
	}
	if rec.Info["DP"] != "100" || len(rec.ID) != 1 || rec.ID[0] != "rs1" {
		t.Fatal("info/id wrong")
	}
	if len(rec.Samples) != 2 {
		t.Fatalf("samples = %d", len(rec.Samples))
	}
	if !rec.Samples[0].Phased || !reflect.DeepEqual(rec.Samples[0].GT, []int{0, 1}) {
		t.Fatalf("s1 = %+v", rec.Samples[0])
	}
	if rec.Samples[0].Fields["DP"] != "30" {
		t.Fatal("format field missing")
	}
	if rec.Samples[1].Phased || !reflect.DeepEqual(rec.Samples[1].GT, []int{1, 2}) {
		t.Fatalf("s2 = %+v", rec.Samples[1])
	}

	rec2, err := p.Next()
	if err != nil {
		t.Fatal(err)
	}
	if rec2.Alt[0] != "<DEL>" || rec2.HasQual {
		t.Fatalf("symbolic line wrong: %+v", rec2)
	}
	if !reflect.DeepEqual(rec2.Samples[0].GT, []int{-1, -1}) {
		t.Fatalf("missing GT = %v", rec2.Samples[0].GT)
	}

	if _, err := p.Next(); err != io.EOF {
		t.Fatalf("want EOF, got %v", err)
	}
}

func TestParser_RejectBadLine(t *testing.T) {
	p := NewParser(strings.NewReader("1\tx\t.\tA\tT\t.\t.\t.\n"))
	if _, err := p.Next(); err == nil {
		t.Fatal("expected parse error")
	}
}
