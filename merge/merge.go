// Package merge 实现跨来源规范记录的确定性去重与合并。
//
// 合并分组键为 chrom:pos:ref（规范坐标 + 参考等位）。同组记录合并为一条
// 多 ALT 记录，ALT 按字典序升序排列，因此结果与记录到达顺序无关。
package merge

import (
	"sort"
	"strconv"
	"strings"

	"github.com/guoyyit/variant-normalization-go/config"
	"github.com/guoyyit/variant-normalization-go/variant"
)

// Merger 按既定策略合并同一位置的规范记录。
type Merger struct {
	cfg *config.Config
}

// NewMerger 创建合并器。
func NewMerger(cfg *config.Config) *Merger {
	if cfg == nil {
		cfg = config.Default()
	}
	return &Merger{cfg: cfg}
}

// GroupKey 返回记录的合并分组键：chrom:pos:ref。
func GroupKey(rec *variant.Record) string { return rec.GroupKey() }

// Merge 合并同一分组内的规范记录，结果与输入顺序无关。
//
// 策略：
//   - 相同规范 ALT（chrom:pos:ref:alt）只保留一份；样本按样本名跨来源去重；
//   - 不同 ALT 在同一分组内合并为多 ALT 记录，ALT 字典序升序；
//   - 以来源优先级最高（数值最小，平手取来源 ID 字典序最小）的记录为基底；
//   - ID/FILTER 取并集；QUAL 取最大；INFO 中可派生统计字段一律重算；
//   - 低于配置 MinQual 的结果 FILTER 标记 LowQual；
//   - 所有来源的 Provenance 全量保留并按 (source,line) 排序。
func (m *Merger) Merge(recs []*variant.Canonical) *variant.Canonical {
	if len(recs) == 0 {
		return nil
	}

	// 1. 选出基底记录：优先级最高，平手来源 ID 最小，再平手 ALT 最小。
	base := recs[0]
	for _, c := range recs[1:] {
		if m.less(c, base) {
			base = c
		}
	}

	// 2. 收集去重后的 ALT 集合（排序保证确定性）。
	altSet := map[string]struct{}{}
	for _, c := range recs {
		for _, a := range c.Record.Alt {
			altSet[a] = struct{}{}
		}
	}
	alts := make([]string, 0, len(altSet))
	for a := range altSet {
		alts = append(alts, a)
	}
	sort.Strings(alts)
	altIndex := map[string]int{}
	for i, a := range alts {
		altIndex[a] = i + 1
	}

	out := *base
	rec := out.Record
	rec.Alt = alts

	// 3. Provenance 全量保留、去重、确定性排序。
	var provs []variant.Provenance
	seenProv := func(p variant.Provenance) bool {
		for _, q := range provs {
			if q.SourceID == p.SourceID && q.Line == p.Line && eqStrings(q.ID, p.ID) {
				return true
			}
		}
		return false
	}
	for _, c := range recs {
		for _, p := range c.Sources {
			if seenProv(p) {
				continue
			}
			provs = append(provs, p)
		}
	}
	sort.Slice(provs, func(i, j int) bool {
		if provs[i].SourceID != provs[j].SourceID {
			return provs[i].SourceID < provs[j].SourceID
		}
		return provs[i].Line < provs[j].Line
	})

	// 4. ID 与 FILTER 取并集。
	ids := collect(recs, func(c *variant.Canonical) []string { return c.Record.ID })
	for _, c := range recs {
		for _, p := range c.Sources {
			ids = append(ids, p.ID...)
		}
	}
	rec.ID = unionStrings(ids)
	rec.Filter = unionStrings(collect(recs, func(c *variant.Canonical) []string { return c.Record.Filter }))

	// 5. QUAL 取最大。
	rec.Qual, rec.HasQual = 0, false
	for _, c := range recs {
		if c.Record.HasQual && (!rec.HasQual || c.Record.Qual > rec.Qual) {
			rec.Qual, rec.HasQual = c.Record.Qual, true
		}
	}

	// 6. 样本合并与基因型重映射到统一 ALT 索引。
	rec.Samples, rec.Format = mergeSamples(recs, altIndex)

	// 7. INFO：保留非统计类字段，删除所有可派生历史字段后重算。
	rec.Info = mergeInfo(recs)
	recomputeInfo(&rec)

	// 8. 质量阈值过滤。
	if rec.HasQual && m.cfg.MinQual > 0 && rec.Qual < m.cfg.MinQual {
		rec.Filter = setFilter(rec.Filter, "LowQual")
	} else {
		rec.Filter = unsetFilter(rec.Filter, "LowQual")
	}
	if len(rec.Filter) == 0 {
		rec.Filter = []string{"PASS"}
	}

	out.Record = rec
	out.Sources = provs
	out.ConfigVersion = m.cfg.ConfigVersion
	return &out
}

// less 定义记录的确定性优先序。
func (m *Merger) less(a, b *variant.Canonical) bool {
	sa := sourceOf(a)
	sb := sourceOf(b)
	pa, pb := m.cfg.Priority(sa), m.cfg.Priority(sb)
	if pa != pb {
		return pa < pb
	}
	if sa != sb {
		return sa < sb
	}
	if a.Record.Pos != b.Record.Pos {
		return a.Record.Pos < b.Record.Pos
	}
	return strings.Join(a.Record.Alt, ",") < strings.Join(b.Record.Alt, ",")
}

func eqStrings(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func sourceOf(c *variant.Canonical) string {
	if len(c.Sources) > 0 {
		return c.Sources[0].SourceID
	}
	return c.Record.Source.ID
}

func collect(recs []*variant.Canonical, f func(*variant.Canonical) []string) []string {
	var all []string
	for _, c := range recs {
		all = append(all, f(c)...)
	}
	return all
}

func unionStrings(in []string) []string {
	seen := map[string]struct{}{}
	var out []string
	for _, s := range in {
		if s == "" || s == "." || s == "PASS" {
			continue
		}
		if _, ok := seen[s]; ok {
			continue
		}
		seen[s] = struct{}{}
		out = append(out, s)
	}
	sort.Strings(out)
	return out
}

func setFilter(filters []string, name string) []string {
	for _, f := range filters {
		if f == name {
			return filters
		}
	}
	return append(filters, name)
}

func unsetFilter(filters []string, name string) []string {
	out := filters[:0]
	for _, f := range filters {
		if f != name && f != "PASS" {
			out = append(out, f)
		}
	}
	return out
}

// recomputeInfo 依据合并后记录重算 AC、AF、AN、NS 等统计字段。
func recomputeInfo(rec *variant.Record) {
	for _, k := range []string{"AC", "AF", "AN", "NS"} {
		delete(rec.Info, k)
	}
	ac := make([]int, len(rec.Alt))
	an := 0
	for _, s := range rec.Samples {
		called := false
		for _, a := range s.GT.Alleles {
			if a >= 1 && a <= len(rec.Alt) {
				ac[a-1]++
				an++
				called = true
			} else if a == 0 {
				an++
				called = true
			}
		}
		_ = called
	}
	if len(rec.Samples) > 0 {
		rec.Info["NS"] = itoa(len(rec.Samples))
		acStr := make([]string, len(ac))
		afStr := make([]string, len(ac))
		for i, v := range ac {
			acStr[i] = itoa(v)
			if an > 0 {
				afStr[i] = strconvFormat(float64(v) / float64(an))
			} else {
				afStr[i] = "0"
			}
		}
		rec.Info["AC"] = strings.Join(acStr, ",")
		rec.Info["AN"] = itoa(an)
		rec.Info["AF"] = strings.Join(afStr, ",")
	}
}

// mergeSamples 按样本名跨来源去重并把 GT 映射到合并后的 ALT 索引。
func mergeSamples(recs []*variant.Canonical, altIndex map[string]int) ([]variant.Sample, []string) {
	byName := map[string]variant.Sample{}
	var order []string
	formatSet := map[string]struct{}{}
	var format []string

	for _, c := range recs {
		for _, f := range c.Record.Format {
			if _, ok := formatSet[f]; ok {
				continue
			}
			formatSet[f] = struct{}{}
			format = append(format, f)
		}
		for _, s := range c.Record.Samples {
			if _, exists := byName[s.Name]; exists {
				// 同名样本已由更高优先级来源提供，保留先到（基底优先）。
				continue
			}
			ns := variant.Sample{Name: s.Name, Fields: map[string]string{}}
			// 该来源记录里的 ALT 下标 -> 合并后下标。
			localToGlobal := map[int]int{0: 0}
			for li, la := range c.Record.Alt {
				localToGlobal[li+1] = altIndex[la]
			}
			ng := variant.Genotype{Phased: s.GT.Phased, Alleles: make([]int, len(s.GT.Alleles))}
			for i, a := range s.GT.Alleles {
				if a < 0 {
					ng.Alleles[i] = -1
				} else if g, ok := localToGlobal[a]; ok {
					ng.Alleles[i] = g
				} else {
					ng.Alleles[i] = -1
				}
			}
			ns.GT = ng
			for k, v := range s.Fields {
				if k == "GT" {
					ns.Fields["GT"] = writeGT(ng)
					continue
				}
				ns.Fields[k] = v
			}
			byName[s.Name] = ns
			order = append(order, s.Name)
		}
	}
	sort.Strings(order)
	samples := make([]variant.Sample, 0, len(order))
	for _, name := range order {
		samples = append(samples, byName[name])
	}
	// FORMAT 固定顺序：GT 优先，其余字典序。
	_, hasGT := formatSet["GT"]
	var rest []string
	for _, f := range format {
		if f != "GT" {
			rest = append(rest, f)
		}
	}
	sort.Strings(rest)
	if hasGT {
		format = append([]string{"GT"}, rest...)
	} else {
		format = rest
	}
	return samples, format
}

// mergeInfo 保留各来源 INFO 中的非派生字段；基底来源优先，冲突时丢弃派生字段。
func mergeInfo(recs []*variant.Canonical) map[string]string {
	derived := map[string]struct{}{
		"AC": {}, "AF": {}, "AN": {}, "NS": {},
	}
	out := map[string]string{}
	for _, c := range recs {
		for k, v := range c.Record.Info {
			if _, drop := derived[k]; drop {
				continue
			}
			if _, exists := out[k]; !exists {
				out[k] = v
			}
		}
	}
	return out
}

func writeGT(g variant.Genotype) string {
	if len(g.Alleles) == 0 {
		return "."
	}
	sep := "/"
	if g.Phased {
		sep = "|"
	}
	parts := make([]string, len(g.Alleles))
	for i, a := range g.Alleles {
		if a < 0 {
			parts[i] = "."
		} else {
			parts[i] = itoa(a)
		}
	}
	return strings.Join(parts, sep)
}

func itoa(n int) string {
	return strconv.Itoa(n)
}

func strconvFormat(v float64) string {
	// 频率统一保留 6 位小数，保证结果确定。
	return strconv.FormatFloat(v, 'f', 6, 64)
}
