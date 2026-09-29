package stats

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// Recompute recomputes INFO/AC,AN,AF, QUAL and Filter from record contents.
func Recompute(rec *variant.Record, cfg config.Config) {
	if rec.Info == nil {
		rec.Info = map[string]string{}
	}

	altCount := len(rec.Alt)
	ac := make([]int, altCount)
	an := 0
	totalDP := 0

	samples := append([]variant.SampleCall(nil), rec.Samples...)
	sort.SliceStable(samples, func(i, j int) bool { return samples[i].Sample < samples[j].Sample })

	for _, sm := range samples {
		for _, allele := range sm.GT {
			if allele < 0 {
				continue
			}
			an++
			if allele >= 1 && allele <= altCount {
				ac[allele-1]++
			}
		}
		if dp := formatInt(sm.Fields, "DP"); dp > 0 {
			totalDP += dp
		}
	}

	acParts := make([]string, altCount)
	afParts := make([]string, altCount)
	for i, count := range ac {
		acParts[i] = strconv.Itoa(count)
		if an > 0 {
			afParts[i] = strconv.FormatFloat(float64(count)/float64(an), 'f', 6, 64)
		} else {
			afParts[i] = "0"
		}
	}
	rec.Info["AC"] = strings.Join(acParts, ",")
	rec.Info["AN"] = strconv.Itoa(an)
	rec.Info["AF"] = strings.Join(afParts, ",")
	if totalDP > 0 {
		rec.Info["DP"] = strconv.Itoa(totalDP)
	}

	// MinQual threshold applies a deterministic filter code while keeping
	// any pre-existing concrete codes.
	if cfg.MinQual > 0 && rec.HasQual && rec.Qual < cfg.MinQual {
		rec.Filter = addFilter(rec.Filter, "MinQual")
	} else {
		rec.Filter = removeFilter(rec.Filter, "MinQual")
	}

	if cfg.MinAF > 0 {
		maxAF := 0.0
		for _, value := range afParts {
			if af, err := strconv.ParseFloat(value, 64); err == nil && af > maxAF {
				maxAF = af
			}
		}
		if maxAF < cfg.MinAF {
			rec.Filter = addFilter(rec.Filter, "MinAF")
		} else {
			rec.Filter = removeFilter(rec.Filter, "MinAF")
		}
	}

	if len(rec.Filter) == 0 {
		rec.Filter = []string{"PASS"}
	}

	rec.ConfigID = cfg.ID
}

func formatInt(fields map[string]string, key string) int {
	if fields == nil {
		return 0
	}
	v, ok := fields[key]
	if !ok {
		return 0
	}
	n, err := strconv.Atoi(v)
	if err != nil {
		return 0
	}
	return n
}

func addFilter(filters []string, code string) []string {
	for _, f := range filters {
		if f == code {
			return filters
		}
	}
	out := append([]string{}, filters...)
	out = append(out, code)
	sort.Strings(out)
	return out
}

func removeFilter(filters []string, code string) []string {
	out := filters[:0:0]
	for _, f := range filters {
		if f != code && f != "PASS" {
			out = append(out, f)
		}
	}
	sort.Strings(out)
	return out
}

// Summary renders one human-readable line describing the merged record,
// used by the CLI and audit logs.
func Summary(rec variant.Record) string {
	return fmt.Sprintf("%s:%d %s>%s sources=%d AC=%s AN=%s AF=%s FILTER=%s",
		rec.Chrom, rec.Pos, rec.Ref, strings.Join(rec.Alt, "/"),
		len(rec.Sources), rec.Info["AC"], rec.Info["AN"], rec.Info["AF"],
		strings.Join(rec.Filter, ";"))
}
