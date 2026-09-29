package importer

import (
	"sort"

	"github.com/guoyyit/variant-normalization-go/merge"
	"github.com/guoyyit/variant-normalization-go/normalize"
	"github.com/guoyyit/variant-normalization-go/split"
	"github.com/guoyyit/variant-normalization-go/variant"
	"github.com/guoyyit/variant-normalization-go/vcfio"
)

// pipeline 串联 拆分 -> 规范化 -> 包装为带来源的 Canonical。
type pipeline struct {
	splitter func(*variant.Record, *vcfio.Header) ([]*variant.Record, error)
	norm     *normalize.Normalizer
}

func newPipeline(norm *normalize.Normalizer) *pipeline {
	return &pipeline{splitter: split.Split, norm: norm}
}

// process 处理一条原始 VCF 记录，返回规范化后的单等位 Canonical 列表。
func (p *pipeline) process(rec *variant.Record, header *vcfio.Header) ([]*variant.Canonical, error) {
	children, err := p.splitter(rec, header)
	if err != nil {
		return nil, err
	}
	out := make([]*variant.Canonical, 0, len(children))
	for _, child := range children {
		nr, err := p.norm.Normalize(child)
		if err != nil {
			return nil, err
		}
		out = append(out, &variant.Canonical{
			Record: *nr,
			Sources: []variant.Provenance{{
				SourceID: rec.Source.ID,
				Line:     rec.Source.Line,
				ID:       append([]string(nil), rec.ID...),
			}},
		})
	}
	return out, nil
}

// dedupMerge 按 GroupKey 排序并合并一批规范化记录，结果顺序确定。
func dedupMerge(mg *merge.Merger, batch []*variant.Canonical) []*variant.Canonical {
	groups := map[string][]*variant.Canonical{}
	var keys []string
	for _, c := range batch {
		k := merge.GroupKey(&c.Record)
		if _, ok := groups[k]; !ok {
			keys = append(keys, k)
		}
		groups[k] = append(groups[k], c)
	}
	sort.Strings(keys)
	out := make([]*variant.Canonical, 0, len(keys))
	for _, k := range keys {
		out = append(out, mg.Merge(groups[k]))
	}
	return out
}
