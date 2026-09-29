package normalize

import (
	"errors"
	"testing"

	"github.com/guoyyit/variant-normalization-go/reference"
	"github.com/guoyyit/variant-normalization-go/variant"
)

// repeatRef: 位置 3..8 为 ACACAC 重复区。
func repeatRef() *reference.Memory {
	return reference.NewMemory(map[string]string{"chr1": "GGACACACATT"})
}

func TestNormalize_SNPTrim(t *testing.T) {
	n := New(nil)
	// TCA->TGA 应右修剪为 C->G（POS 2）
	got, err := n.Normalize(&variant.Record{Chrom: "c", Pos: 1, Ref: "TCA", Alt: []string{"TGA"}})
	if err != nil {
		t.Fatal(err)
	}
	if got.Pos != 2 || got.Ref != "C" || got.Alt[0] != "G" {
		t.Fatalf("got pos=%d ref=%s alt=%v", got.Pos, got.Ref, got.Alt)
	}
}

func TestNormalize_DeletionLeftAlign(t *testing.T) {
	// AC 重复区中的 4bp 删除：ACACA(pos3..7)->A。
	// 无论从 pos3 还是从右侧等价位置 pos5 描述，都应归一到同一规范表示。
	n := New(repeatRef())
	got, err := n.Normalize(&variant.Record{Chrom: "chr1", Pos: 3, Ref: "ACACA", Alt: []string{"A"}})
	if err != nil {
		t.Fatal(err)
	}
	right, err := n.Normalize(&variant.Record{Chrom: "chr1", Pos: 5, Ref: "ACACA", Alt: []string{"A"}})
	if err != nil {
		t.Fatal(err)
	}
	if got.Pos != right.Pos || got.Ref != right.Ref || got.Alt[0] != right.Alt[0] {
		t.Fatalf("deletion not left-aligned: %d:%s/%s vs %d:%s/%s",
			got.Pos, got.Ref, got.Alt[0], right.Pos, right.Ref, right.Alt[0])
	}
	// 规范结果必须最左：pos2 GACAC/G（左邻 G 与 alt 首碱基 A 不同）
	if got.Pos != 2 || got.Ref != "GACAC" || got.Alt[0] != "G" {
		t.Fatalf("got pos=%d ref=%s alt=%v", got.Pos, got.Ref, got.Alt)
	}
}

func TestNormalize_InsertionLeftAlign(t *testing.T) {
	// AC 重复区中插入 AC：C->CAC，从 pos4 与 pos6 两种等价写法应一致。
	n := New(repeatRef())
	got, err := n.Normalize(&variant.Record{Chrom: "chr1", Pos: 4, Ref: "C", Alt: []string{"CAC"}})
	if err != nil {
		t.Fatal(err)
	}
	right, err := n.Normalize(&variant.Record{Chrom: "chr1", Pos: 6, Ref: "C", Alt: []string{"CAC"}})
	if err != nil {
		t.Fatal(err)
	}
	if got.Pos != right.Pos || got.Ref != right.Ref || got.Alt[0] != right.Alt[0] {
		t.Fatalf("insertion not aligned: %d:%s/%s vs %d:%s/%s",
			got.Pos, got.Ref, got.Alt[0], right.Pos, right.Ref, right.Alt[0])
	}
	if got.Pos != 2 || got.Ref != "G" || got.Alt[0] != "GAC" {
		t.Fatalf("got pos=%d ref=%s alt=%v", got.Pos, got.Ref, got.Alt)
	}
}

func TestNormalize_Deterministic(t *testing.T) {
	// 重复规范化保持幂等。
	n := New(repeatRef())
	a, err := n.Normalize(&variant.Record{Chrom: "chr1", Pos: 5, Ref: "ACACA", Alt: []string{"A"}})
	if err != nil {
		t.Fatal(err)
	}
	b, err := n.Normalize(a)
	if err != nil {
		t.Fatal(err)
	}
	if a.Pos != b.Pos || a.Ref != b.Ref || a.Alt[0] != b.Alt[0] {
		t.Fatalf("normalization not idempotent: %v vs %v", a, b)
	}
}

func TestNormalize_InvalidFormat(t *testing.T) {
	n := New(nil)
	cases := []*variant.Record{
		{Chrom: "c", Pos: 1, Ref: "AC", Alt: []string{"AZ"}},
		{Chrom: "c", Pos: 0, Ref: "A", Alt: []string{"T"}},
		{Chrom: "c", Pos: 1, Ref: "A", Alt: []string{"<DEL>"}},
		{Chrom: "c", Pos: 1, Ref: "", Alt: []string{"T"}},
		{Chrom: "c", Pos: 1, Ref: "A", Alt: nil},
	}
	for i, rec := range cases {
		_, err := n.Normalize(rec)
		var ve *variant.ValidationError
		if !errors.As(err, &ve) || ve.Kind != variant.KindInvalidFormat {
			t.Fatalf("case %d: want invalid_format, got %v", i, err)
		}
	}
}

func TestNormalize_OutOfRange(t *testing.T) {
	n := New(repeatRef())
	_, err := n.Normalize(&variant.Record{Chrom: "chr9", Pos: 1, Ref: "A", Alt: []string{"T"}})
	var ve *variant.ValidationError
	if !errors.As(err, &ve) || ve.Kind != variant.KindOutOfRange {
		t.Fatalf("unknown chrom: got %v", err)
	}
	_, err = n.Normalize(&variant.Record{Chrom: "chr1", Pos: 100, Ref: "AA", Alt: []string{"T"}})
	if !errors.As(err, &ve) || ve.Kind != variant.KindOutOfRange {
		t.Fatalf("past contig end: got %v", err)
	}
	// REF 与参考不一致属于格式非法
	_, err = n.Normalize(&variant.Record{Chrom: "chr1", Pos: 3, Ref: "C", Alt: []string{"T"}})
	if !errors.As(err, &ve) || ve.Kind != variant.KindInvalidFormat {
		t.Fatalf("ref mismatch: got %v", err)
	}

	// 多等位记录应先拆分
	_, err = n.Normalize(&variant.Record{Chrom: "chr1", Pos: 1, Ref: "A", Alt: []string{"T", "C"}})
	if !errors.As(err, &ve) || ve.Kind != variant.KindInvalidFormat {
		t.Fatalf("multiallelic: got %v", err)
	}
}
