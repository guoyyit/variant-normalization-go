package variant

import "strconv"

// ErrorKind 区分记录被拒绝的原因。
type ErrorKind string

const (
	// KindInvalidFormat 记录格式非法（字段缺失、碱基非法、基因型无法解析等）。
	KindInvalidFormat ErrorKind = "invalid_format"
	// KindOutOfRange 坐标越界（POS 超出参考序列范围等）。
	KindOutOfRange ErrorKind = "out_of_range"
)

// ValidationError 描述一条记录被拒绝的可区分原因。
type ValidationError struct {
	Kind ErrorKind
	Msg  string
}

func (e *ValidationError) Error() string {
	return string(e.Kind) + ": " + e.Msg
}

// Genotype 表示一个样本的基因型，Alleles 为 [REF, ALT...] 中的等位基因下标，-1 表示缺失。
type Genotype struct {
	Alleles []int
	Phased  bool
}

// Sample 表示一条记录中的单个样本。
type Sample struct {
	Name   string
	GT     Genotype
	Fields map[string]string
}

// Record 表示一条 VCF 文本记录（可能是多等位的）。
type Record struct {
	Chrom   string
	Pos     int
	ID      []string
	Ref     string
	Alt     []string
	Qual    float64
	HasQual bool
	Filter  []string
	Info    map[string]string
	Format  []string
	Samples []Sample
	Source  Source
}

// Source 标识记录的来源与在来源中的位置。
type Source struct {
	ID   string
	Line int64
}

// Provenance 记录合并后保留的来源追溯信息。
type Provenance struct {
	SourceID string
	Line     int64
	ID       []string
}

// Canonical 是规范化、合并后的入库记录。
type Canonical struct {
	Record  Record
	Sources []Provenance
	// ConfigVersion 记录生成该合并结果所用配置的版本。
	ConfigVersion string
}

// Key 返回记录的全局排序与去重键（同一变异的规范表示）。
func (r *Record) Key() string {
	return r.Chrom + "\x00" + sortPos(r.Pos) + "\x00" + r.Ref + "\x00" + joinAlt(r.Alt)
}

// GroupKey 返回同位置合并分组键：chrom:pos:ref（ALT 不同的记录在此合并）。
func (r *Record) GroupKey() string {
	return r.Chrom + "\x00" + sortPos(r.Pos) + "\x00" + r.Ref
}

// sortPos 把坐标编码为零填充定宽字符串，使字符串字典序与数值序一致。
// 14 位足以覆盖所有现实基因组坐标；负数前缀符号仅用于健壮性。
func sortPos(pos int) string {
	if pos < 0 {
		return "-" + padLeft(strconv.Itoa(-pos), 13)
	}
	return padLeft(strconv.Itoa(pos), 14)
}

func padLeft(s string, width int) string {
	if len(s) >= width {
		return s
	}
	b := make([]byte, width-len(s), width)
	for i := range b {
		b[i] = '0'
	}
	return string(append(b, s...))
}

func joinAlt(alts []string) string {
	out := ""
	for i, a := range alts {
		if i > 0 {
			out += ","
		}
		out += a
	}
	return out
}
