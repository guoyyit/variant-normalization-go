package normalize

import (
	"fmt"
	"strings"

	"github.com/guoyyit/variant-normalization-go/reference"
	"github.com/guoyyit/variant-normalization-go/variant"
)

// Normalizer 对 VCF 记录做等位基因修剪与左对齐。
type Normalizer struct {
	ref reference.Genome
}

// New 创建 Normalizer；ref 用于左对齐时读取上游参考碱基。
func New(ref reference.Genome) *Normalizer {
	return &Normalizer{ref: ref}
}

func invalidf(format string, args ...any) error {
	return &variant.ValidationError{Kind: variant.KindInvalidFormat, Msg: fmt.Sprintf(format, args...)}
}

func outOfRangef(format string, args ...any) error {
	return &variant.ValidationError{Kind: variant.KindOutOfRange, Msg: fmt.Sprintf(format, args...)}
}

var bases = map[byte]bool{'A': true, 'C': true, 'G': true, 'T': true, 'N': true}

func validAllele(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		if !bases[s[i]] {
			return false
		}
	}
	return true
}

// Normalize 返回记录的规范表示；非法或越界时返回 *variant.ValidationError。
func (n *Normalizer) Normalize(in *variant.Record) (*variant.Record, error) {
	rec := *in // 拷贝，避免修改输入
	rec.ID = append([]string(nil), in.ID...)
	rec.Alt = append([]string(nil), in.Alt...)
	rec.Filter = append([]string(nil), in.Filter...)
	rec.Format = append([]string(nil), in.Format...)

	if rec.Chrom == "" {
		return nil, invalidf("missing CHROM")
	}
	if rec.Pos < 1 {
		return nil, invalidf("POS must be >= 1, got %d", rec.Pos)
	}
	if len(rec.Ref) == 0 {
		return nil, invalidf("missing REF")
	}
	rec.Ref = strings.ToUpper(rec.Ref)
	if !validAllele(rec.Ref) {
		return nil, invalidf("invalid REF allele %q", in.Ref)
	}
	if len(rec.Alt) == 0 {
		return nil, invalidf("missing ALT")
	}
	for i, a := range rec.Alt {
		if a == "." || a == "" {
			return nil, invalidf("missing ALT[%d]", i)
		}
		rec.Alt[i] = strings.ToUpper(a)
		if isSymbolic(rec.Alt[i]) {
			return nil, invalidf("symbolic allele not supported: %q", a)
		}
		if !validAllele(rec.Alt[i]) {
			return nil, invalidf("invalid ALT allele %q", a)
		}
	}
	if len(rec.Alt) != 1 {
		return nil, invalidf("normalizer expects a single-allelic record; split first")
	}

	if n.ref != nil {
		if n.ref.Length(rec.Chrom) == 0 {
			return nil, outOfRangef("unknown chromosome %q", rec.Chrom)
		}
		if rec.Pos+len(rec.Ref)-1 > n.ref.Length(rec.Chrom) {
			return nil, outOfRangef("%s:%d-%d exceeds contig length %d",
				rec.Chrom, rec.Pos, rec.Pos+len(rec.Ref)-1, n.ref.Length(rec.Chrom))
		}
		refBase, err := n.ref.Prefix(rec.Chrom, rec.Pos+len(rec.Ref)-1, len(rec.Ref))
		if err != nil {
			return nil, outOfRangef("cannot read reference: %v", err)
		}
		if refBase != rec.Ref {
			return nil, invalidf("REF %q does not match reference %q at %s:%d",
				rec.Ref, refBase, rec.Chrom, rec.Pos)
		}
	}

	rec.Ref, rec.Alt[0], rec.Pos = n.leftAlign(rec.Chrom, rec.Pos, rec.Ref, rec.Alt[0])

	return &rec, nil
}

// leftAlign 对单个 (REF, ALT) 执行等位基因修剪与左对齐。
//
// 采用 Tan/Abecasis 等“规范化”的标准迭代过程（vt / bcftools norm 等价）：
//  1. 去共同后缀（右 parsimony，每条等位保留至少 1 个锚定碱基）；
//  2. 当各等位首碱基相同且未到染色体起点时，向上游扩展一个碱基窗口
//     （REF 与 ALT 同时前置参考碱基），再次去共同后缀；首碱基不同即终止；
//  3. 去共同前缀（左 parsimony，POS 相应右移）。
//
// 终止时首碱基不同保证了结果是最左、最短表示，且与输入写法无关。
func (n *Normalizer) leftAlign(chrom string, pos int, ref, alt string) (string, string, int) {
	trimRight := func(r, s string) (string, string) {
		for len(r) > 1 && len(s) > 1 && r[len(r)-1] == s[len(s)-1] {
			r, s = r[:len(r)-1], s[:len(s)-1]
		}
		return r, s
	}

	ref, alt = trimRight(ref, alt)
	for pos > 1 && ref[0] == alt[0] {
		if n.ref == nil {
			break
		}
		up, err := n.ref.Prefix(chrom, pos-1, 1)
		if err != nil {
			break
		}
		up = strings.ToUpper(up)
		pos--
		ref, alt = trimRight(up+ref, up+alt)
	}
	for len(ref) > 1 && len(alt) > 1 && ref[0] == alt[0] {
		ref, alt = ref[1:], alt[1:]
		pos++
	}
	return ref, alt, pos
}

func isSymbolic(a string) bool {
	switch {
	case strings.Contains(a, "<") || strings.Contains(a, ">"):
		return true
	case strings.Contains(a, "[") || strings.Contains(a, "]"):
		return true // 断裂端 / 远程插入
	case a == "*":
		return true
	}
	return false
}
