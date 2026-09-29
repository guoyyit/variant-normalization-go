package store

import (
	"fmt"
	"sync"
	"testing"
	"time"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// TestStore_ConcurrentWriteQuery verifies readers always observe complete,
// self-consistent records while writes+flushes+compaction happen.
func TestStore_ConcurrentWriteQuery(t *testing.T) {
	dir := t.TempDir()
	st, err := Open(dir, func(existing, incoming variant.Record) variant.Record {
		if len(existing.Sources) == 0 {
			out := incoming.Clone()
			out.Info = map[string]string{"ok": "1"}
			return out
		}
		out := existing.Clone()
		out.Sources = append(out.Sources, incoming.Sources...)
		return out
	}, WithMemLimit(64))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()

	const writers = 4
	const perWriter = 40
	var wg sync.WaitGroup

	for w := 0; w < writers; w++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			for i := 0; i < perWriter; i++ {
				rec := variant.Record{
					Chrom: "chr1", Pos: i + 1, Ref: "A", Alt: []string{"T"},
					Sources: []variant.SourceRef{{
						Source: fmt.Sprintf("src-%d", id), File: "f.vcf", Line: int64(i)}},
				}
				if err := st.Upsert(rec); err != nil {
					t.Errorf("upsert: %v", err)
					return
				}
			}
		}(w)
	}

	stop := make(chan struct{})
	for q := 0; q < 4; q++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			for i := 0; i < 200; i++ {
				select {
				case <-stop:
					return
				default:
				}
				page, err := st.Query(PageQuery{Chrom: "chr1", Start: 1, End: perWriter, PageSize: 32})
				if err != nil {
					t.Errorf("query: %v", err)
					return
				}
				for _, rec := range page.Records {
					// Every observed record must have the post-fold marker.
					if rec.Info["ok"] != "1" {
						t.Errorf("observed intermediate state for %s:%d info=%v",
							rec.Chrom, rec.Pos, rec.Info)
						return
					}
				}
				time.Sleep(time.Duration(id+1) * time.Millisecond)
			}
		}(q)
	}

	wg.Wait()
	close(stop)

	count, err := st.Count()
	if err != nil {
		t.Fatal(err)
	}
	if count != perWriter {
		t.Fatalf("count=%d want %d", count, perWriter)
	}

	// Every position must retain all writer provenance after settling.
	page, err := st.Query(PageQuery{Chrom: "chr1", Start: 1, End: perWriter, PageSize: perWriter})
	if err != nil {
		t.Fatal(err)
	}
	for _, rec := range page.Records {
		if len(rec.Sources) != writers {
			t.Fatalf("pos %d sources=%d want %d", rec.Pos, len(rec.Sources), writers)
		}
	}
}
