package store

import (
	"testing"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

func BenchmarkUpsertSameKey(b *testing.B) {
	st, err := Open(b.TempDir(), func(e, i variant.Record) variant.Record {
		if len(e.Sources) == 0 {
			return i
		}
		e.Sources = append(e.Sources, i.Sources...)
		return e
	}, WithMemLimit(64))
	if err != nil {
		b.Fatal(err)
	}
	defer st.Close()
	rec := variant.Record{Chrom: "c", Pos: 1, Ref: "A", Alt: []string{"T"}}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r := rec
		r.Sources = []variant.SourceRef{{Source: "s", File: "f", Line: int64(i)}}
		if err := st.Upsert(r); err != nil {
			b.Fatal(err)
		}
	}
}
