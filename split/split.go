package split

import (
	"strconv"
	"strings"

	"github.com/guoyyit/variant-normalization-go/variant"
	"github.com/guoyyit/variant-normalization-go/vcfio"
)

// Split 把一条多等位记录拆成多条单等位记录，并完成基因型与 FORMAT 字段重映射。
func Split(rec *variant.Record, header *vcfio.Header) ([]*variant.Record, error) {
	if len(rec.Alt) == 0 {
		return nil, &variant.ValidationError{Kind: variant.KindInvalidFormat, Msg: "cannot split record without ALT"}
	}
	if len(rec.Alt) == 1 {
		return []*variant.Record{rec}, nil
	}

	numberOf := func(field string) string {
		if header != nil {
			return header.FormatNumber[field]
		}
		return ""
	}

	out := make([]*variant.Record, 0, len(rec.Alt))
	for altIdx := 0; altIdx < len(rec.Alt); altIdx++ {
		child := *rec
		child.Alt = []string{rec.Alt[altIdx]}
		child.ID = append([]string(nil), rec.ID...)
		child.Filter = append([]string(nil), rec.Filter...)
		child.Format = append([]string(nil), rec.Format...)

		samples := make([]variant.Sample, len(rec.Samples))
		for i, s := range rec.Samples {
			cs := variant.Sample{Name: s.Name, Fields: map[string]string{}, GT: remapGT(s.GT, altIdx+1)}
			for _, f := range child.Format {
				if f == "GT" {
					cs.Fields["GT"] = formatGT(cs.GT)
					continue
				}
				raw, ok := s.Fields[f]
				if !ok {
					continue
				}
				cs.Fields[f] = remapField(raw, numberOf(f), altIdx, len(rec.Alt))
			}
			samples[i] = cs
		}
		child.Samples = samples
		out = append(out, &child)
	}
	return out, nil
}

// remapGT 把多等位基因型映射为双等位：保留 0 与目标等位(原值->1)，其余变为缺失(-1)。
func remapGT(gt variant.Genotype, target int) variant.Genotype {
	out := variant.Genotype{Phased: gt.Phased, Alleles: make([]int, len(gt.Alleles))}
	for i, a := range gt.Alleles {
		switch {
		case a < 0:
			out.Alleles[i] = -1
		case a == 0:
			out.Alleles[i] = 0
		case a == target:
			out.Alleles[i] = 1
		default:
			out.Alleles[i] = -1
		}
	}
	return out
}

func formatGT(gt variant.Genotype) string {
	if len(gt.Alleles) == 0 {
		return "."
	}
	sep := "/"
	if gt.Phased {
		sep = "|"
	}
	parts := make([]string, len(gt.Alleles))
	for i, a := range gt.Alleles {
		if a < 0 {
			parts[i] = "."
		} else {
			parts[i] = strconv.Itoa(a)
		}
	}
	return strings.Join(parts, sep)
}

// remapField 按 VCF FORMAT 字段的 Number 声明重映射逗号分隔取值。
// altIdx 为 0-based 目标 ALT 下标，nAlt 为拆分前 ALT 数量。
func remapField(raw, number string, altIdx, nAlt int) string {
	if raw == "." || raw == "" {
		return raw
	}
	parts := strings.Split(raw, ",")
	switch number {
	case "R": // ref + 每个 ALT：保留 ref 与目标 ALT
		if len(parts) != nAlt+1 {
			return raw
		}
		return parts[0] + "," + parts[altIdx+1]
	case "A": // 每个 ALT：只保留目标
		if len(parts) != nAlt {
			return raw
		}
		return parts[altIdx]
	case "G": // 基因型组合数 = (n+1)(n+2)/2，双等位为 3：00,01,11
		if len(parts) != (nAlt+1)*(nAlt+2)/2 {
			return raw
		}
		// VCF Number=G 顺序为 00,01,11,02,12,22...，
		// 按第二等位 j 分组：等位对 (i,j) 下标为 j*(j+1)/2 + i。
		t := altIdx + 1
		return strings.Join([]string{
			parts[0],           // 0/0
			parts[t*(t+1)/2],   // 0/target
			parts[t*(t+1)/2+t], // target/target
		}, ",")
	default:
		return raw
	}
}
