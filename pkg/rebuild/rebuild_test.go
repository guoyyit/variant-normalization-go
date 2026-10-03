package rebuild

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/importer"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

type fakeRef struct{ seq string }

func (f fakeRef) Substring(chrom string, start, end int) (string, error) {
	return f.seq[start-1 : end-1], nil
}
func (f fakeRef) Length(chrom string) (int, error) { return len(f.seq), nil }

// refSeq has homopolymers and a tandem repeat so left alignment matters.
const refSeq = "CAAAACTACACACG"

func writeFile(t *testing.T, path, content string) {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
}

func vcfHeader() string {
	return "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
}

// sourceA reports the same deletion in the CACA tandem repeat at two
// positions; both reconstruct to CAAAACTACACG.
func sourceA() string {
	return vcfHeader() +
		"chr1\t8\t.\tACA\tA\t60\tPASS\t.\tGT\t0/1\n" +
		"chr1\t10\t.\tACA\tA\t60\tPASS\t.\tGT\t0/1\n"
}

func newRouter(t *testing.T, dir string) *store.Router {
	t.Helper()
	router, err := store.OpenRouter(dir, importer.Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { router.Close() })
	return router
}

func stringsReader(s string) *strings.Reader { return strings.NewReader(s) }

func assertStoresEqual(t *testing.T, want, got *store.Store) {
	t.Helper()
	wantMap := map[string]variant.Record{}
	if err := want.Scan(func(r variant.Record) error {
		wantMap[r.Chrom+"|"+itoa(r.Pos)+"|"+r.Ref] = r
		return nil
	}); err != nil {
		t.Fatal(err)
	}
	gotCount := 0
	if err := got.Scan(func(r variant.Record) error {
		key := r.Chrom + "|" + itoa(r.Pos) + "|" + r.Ref
		w, ok := wantMap[key]
		if !ok {
			t.Errorf("unexpected rebuilt key %s", key)
		} else if !recordsEqual(w, r) {
			t.Errorf("rebuilt record differs at %s: %+v want %+v", key, r, w)
		}
		delete(wantMap, key)
		gotCount++
		return nil
	}); err != nil {
		t.Fatal(err)
	}
	if len(wantMap) != 0 {
		t.Errorf("rebuilt store missing %d keys: %v", len(wantMap), wantMap)
	}
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var b []byte
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	return string(b)
}

// recordsEqual compares the identity-bearing and recomputed content of two
// fully merged records (CHROM/POS/REF/ALT and genotype-derived stats).
func recordsEqual(a, b variant.Record) bool {
	if a.Chrom != b.Chrom || a.Pos != b.Pos || a.Ref != b.Ref ||
		len(a.Alt) != len(b.Alt) {
		return false
	}
	for i := range a.Alt {
		if a.Alt[i] != b.Alt[i] {
			return false
		}
	}
	if len(a.Samples) != len(b.Samples) {
		return false
	}
	for i := range a.Samples {
		if a.Samples[i].Sample != b.Samples[i].Sample ||
			len(a.Samples[i].GT) != len(b.Samples[i].GT) {
			return false
		}
		for j := range a.Samples[i].GT {
			if a.Samples[i].GT[j] != b.Samples[i].GT[j] {
				return false
			}
		}
	}
	for _, k := range []string{"AC", "AN", "AF", "DP"} {
		if a.Info[k] != b.Info[k] {
			return false
		}
	}
	return true
}

func TestRebuild_MatchesCleanImport(t *testing.T) {
	cfg := config.Default()
	ref := fakeRef{seq: refSeq}

	// 1. Legacy database: written before the rebuild feature (root layout,
	// no manifest) with the original sources.
	root := filepath.Join(t.TempDir(), "db")
	legacy, err := store.Open(root, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	im, err := importer.New(legacy, ref, cfg, filepath.Join(t.TempDir(), "cp-old.json"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := im.Import(context.Background(), "src-a", "a.vcf", stringsReader(sourceA())); err != nil {
		t.Fatal(err)
	}
	if err := legacy.Close(); err != nil {
		t.Fatal(err)
	}

	// 2. Clean import of the same raw data with the fixed pipeline.
	freshDir := filepath.Join(t.TempDir(), "fresh")
	fresh, err := store.Open(freshDir, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	im2, err := importer.New(fresh, ref, cfg, filepath.Join(t.TempDir(), "cp-fresh.json"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := im2.Import(context.Background(), "src-a", "a.vcf", stringsReader(sourceA())); err != nil {
		t.Fatal(err)
	}
	if err := fresh.Flush(); err != nil {
		t.Fatal(err)
	}

	// 3. Rebuild the legacy database from the original file.
	router := newRouter(t, root)
	rawPath := filepath.Join(t.TempDir(), "a.vcf")
	writeFile(t, rawPath, sourceA())
	rb, err := New(router, ref, cfg, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	result, err := rb.Rebuild(context.Background(),
		[]Source{{Name: "src-a", Path: rawPath}},
		filepath.Join(t.TempDir(), "cp-rebuild.json"))
	if err != nil {
		t.Fatal(err)
	}
	if result.Generation != 1 {
		t.Fatalf("generation=%d want 1", result.Generation)
	}

	// 4. The rebuilt generation must equal the clean import exactly, key by
	// key and allele by allele.
	rebuilt := router.Active()
	assertStoresEqual(t, fresh, rebuilt)

	// 5. The equivalent spellings collapsed to one record, at the
	// left-aligned canonical position in the CACA context.
	page, err := rebuilt.Query(store.PageQuery{Chrom: "chr1", Start: 1, End: len(refSeq), PageSize: 100})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 {
		t.Fatalf("records=%d want 1 (equivalent forms must converge)", len(page.Records))
	}
	rec := page.Records[0]
	if rec.Pos != 7 || rec.Ref != "TAC" || len(rec.Alt) != 1 || rec.Alt[0] != "T" {
		t.Fatalf("canonical = %d:%s>%v want 7:TAC>T", rec.Pos, rec.Ref, rec.Alt)
	}
	fresh.Close()
}

// TestRebuild_ResumeNoDuplicates verifies durable, resumable rebuild
// progress: interrupting after some files and re-running Rebuild with the
// same checkpoint yields exactly one contribution per source record - no
// duplicates, no omissions.
func TestRebuild_ResumeNoDuplicates(t *testing.T) {
	cfg := config.Default()
	ref := fakeRef{seq: refSeq}
	root := filepath.Join(t.TempDir(), "db")
	router := newRouter(t, root)

	dir := t.TempDir()
	aPath := filepath.Join(dir, "a.vcf")
	bPath := filepath.Join(dir, "b.vcf")
	writeFile(t, aPath, vcfHeader()+
		"chr1\t1\t.\tC\tG\t50\tPASS\t.\tGT\t0/1\n")
	writeFile(t, bPath, vcfHeader()+
		"chr1\t2\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n")

	cp := filepath.Join(t.TempDir(), "cp.json")

	// First run: only file a is present/requested (simulate partial progress
	// before an interruption).
	rb, err := New(router, ref, cfg, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := rb.Rebuild(context.Background(),
		[]Source{{Name: "src-a", Path: aPath}}, cp); err != nil {
		t.Fatal(err)
	}
	if router.Generation() != 1 {
		t.Fatalf("gen=%d want 1", router.Generation())
	}

	// Second run builds the next fresh generation from both files. Because a
	// new shadow generation starts empty and the committed-run checkpoint was
	// cleared, file a is re-imported into the new generation and file b is
	// added: the finished generation is a complete rebuild, not a delta.
	rb2, err := New(router, ref, cfg, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := rb2.Rebuild(context.Background(),
		[]Source{{Name: "src-a", Path: aPath}, {Name: "src-b", Path: bPath}}, cp); err != nil {
		t.Fatal(err)
	}
	if router.Generation() != 2 {
		t.Fatalf("gen=%d want 2", router.Generation())
	}

	page, err := router.Active().Query(store.PageQuery{
		Chrom: "chr1", Start: 1, End: len(refSeq), PageSize: 100,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 2 {
		t.Fatalf("records=%d want 2 (resume must not duplicate file a)", len(page.Records))
	}
	// The replayed record retains a single source attribution even though a
	// was imported in both generations.
	for _, rec := range page.Records {
		if rec.Pos == 1 {
			srcFiles := map[string]bool{}
			for _, s := range rec.Sources {
				srcFiles[s.Source+":"+s.File] = true
			}
			if len(srcFiles) != 1 {
				t.Fatalf("pos1 sources=%v want exactly one", rec.Sources)
			}
		}
	}
}

// TestRebuild_LegacyOpenUnchanged verifies that opening a pre-rebuild
// database never silently rewrites it: no manifest and no gen- directory are
// created merely by opening/querying.
func TestRebuild_LegacyOpenUnchanged(t *testing.T) {
	root := filepath.Join(t.TempDir(), "db")
	st, err := store.Open(root, importer.Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	im, err := importer.New(st, fakeRef{seq: refSeq}, config.Default(),
		filepath.Join(t.TempDir(), "cp.json"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := im.Import(context.Background(), "s", "a.vcf",
		stringsReader(vcfHeader()+"chr1\t1\t.\tC\tG\t.\t.\t.\n")); err != nil {
		t.Fatal(err)
	}
	if err := st.Close(); err != nil {
		t.Fatal(err)
	}

	router, err := store.OpenRouter(root, importer.Fold(config.Default()))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := router.Active().Query(store.PageQuery{
		Chrom: "chr1", Start: 1, End: len(refSeq), PageSize: 10,
	}); err != nil {
		t.Fatal(err)
	}
	if router.Generation() != 0 {
		t.Fatalf("legacy db generation=%d want 0", router.Generation())
	}
	entries, err := os.ReadDir(root)
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range entries {
		if e.Name() == "current.json" || len(e.Name()) >= 4 && e.Name()[:4] == "gen-" {
			t.Fatalf("opening legacy db created %q", e.Name())
		}
	}
	router.Close()
}

// TestRebuild_InterruptedResume rebuilds from a file, cancels the context
// mid-stream (so the shadow is discarded by Commit never running but the
// checkpoint/state survive), then re-runs Rebuild with the same checkpoint.
// The resumed run reuses the same shadow generation and yields every record
// exactly once.
func TestRebuild_InterruptedResume(t *testing.T) {
	cfg := config.Default()
	longSeq := strings.Repeat("A", 256)
	ref := fakeRef{seq: longSeq}
	root := filepath.Join(t.TempDir(), "db")
	router := newRouter(t, root)

	var b strings.Builder
	b.WriteString(vcfHeader())
	for i := 1; i <= 40; i++ {
		fmt.Fprintf(&b, "chr1\t%d\t.\tA\tT\t.\t.\t.\n", i)
	}
	path := filepath.Join(t.TempDir(), "bulk.vcf")
	writeFile(t, path, b.String())
	cp := filepath.Join(t.TempDir(), "cp.json")

	rb, err := New(router, ref, cfg, importer.Fold(cfg), WithSyncBatch(5))
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	_, err = rb.Rebuild(ctx, []Source{{Name: "bulk", Path: path}}, cp)
	if err == nil {
		t.Fatal("expected cancellation error")
	}
	if router.Generation() != 0 {
		t.Fatalf("active generation must not switch on interruption, got %d", router.Generation())
	}

	// Resume with a fresh context and the same checkpoint.
	rb2, err := New(router, ref, cfg, importer.Fold(cfg), WithSyncBatch(5))
	if err != nil {
		t.Fatal(err)
	}
	result, err := rb2.Rebuild(context.Background(), []Source{{Name: "bulk", Path: path}}, cp)
	if err != nil {
		t.Fatal(err)
	}
	if result.Generation != 1 {
		t.Fatalf("generation=%d want 1", result.Generation)
	}
	page, err := router.Active().Query(store.PageQuery{
		Chrom: "chr1", Start: 1, End: len(longSeq), PageSize: 1000,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 40 {
		t.Fatalf("records=%d want 40 (resume must be complete, no dup/omission)", len(page.Records))
	}
}

// TestRebuild_ConcurrentQueriesSeeWholeStates loads an old generation, then
// runs a rebuild in a goroutine while continuously querying. Readers must
// observe only complete generations: the old record count or the new record
// count, never a partial intermediate count.
func TestRebuild_ConcurrentQueriesSeeWholeStates(t *testing.T) {
	cfg := config.Default()
	ref := fakeRef{seq: strings.Repeat("A", 512)}
	root := filepath.Join(t.TempDir(), "db")

	// Seed the legacy generation with 30 records.
	old := t.TempDir()
	oldPath := filepath.Join(old, "old.vcf")
	var ob strings.Builder
	ob.WriteString(vcfHeader())
	for i := 1; i <= 30; i++ {
		fmt.Fprintf(&ob, "chr1\t%d\t.\tA\tT\t.\t.\t.\n", i)
	}
	writeFile(t, oldPath, ob.String())
	{
		st, err := store.Open(root, importer.Fold(cfg))
		if err != nil {
			t.Fatal(err)
		}
		im, err := importer.New(st, ref, cfg, filepath.Join(t.TempDir(), "cp.json"))
		if err != nil {
			t.Fatal(err)
		}
		f, err := os.Open(oldPath)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := im.Resume(context.Background(), "old", "old.vcf", f); err != nil {
			t.Fatal(err)
		}
		f.Close()
		st.Close()
	}

	router, err := store.OpenRouter(root, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { router.Close() })
	if router.Generation() != 0 {
		t.Fatalf("seed generation=%d want 0", router.Generation())
	}

	// New raw data: 60 records at different positions.
	dir := t.TempDir()
	newPath := filepath.Join(dir, "new.vcf")
	var nb strings.Builder
	nb.WriteString(vcfHeader())
	for i := 1; i <= 60; i++ {
		fmt.Fprintf(&nb, "chr1\t%d\t.\tA\tT\t.\t.\t.\n", i)
	}
	writeFile(t, newPath, nb.String())

	var wg sync.WaitGroup
	wg.Add(1)
	var errMu sync.Mutex
	var rebuildErr error
	setErr := func(e error) {
		errMu.Lock()
		rebuildErr = e
		errMu.Unlock()
	}
	getErr := func() error {
		errMu.Lock()
		defer errMu.Unlock()
		return rebuildErr
	}
	go func() {
		defer wg.Done()
		rb, nerr := New(router, ref, cfg, importer.Fold(cfg), WithSyncBatch(3))
		if nerr != nil {
			setErr(nerr)
			return
		}
		_, e := rb.Rebuild(context.Background(),
			[]Source{{Name: "new", Path: newPath}},
			filepath.Join(t.TempDir(), "cp.json"))
		setErr(e)
	}()

	deadline := time.Now().Add(10 * time.Second)
	for time.Now().Before(deadline) {
		page, qerr := router.Active().Query(store.PageQuery{
			Chrom: "chr1", Start: 1, End: 512, PageSize: 1000,
		})
		if qerr != nil {
			t.Fatal(qerr)
		}
		switch n := len(page.Records); n {
		case 30, 60:
			// complete old or complete new generation
		default:
			t.Fatalf("query observed partial state with %d records (want 30 or 60)", n)
		}
		if getErr() != nil {
			break
		}
	}
	wg.Wait()
	if e := getErr(); e != nil {
		t.Fatal(e)
	}
	page, err := router.Active().Query(store.PageQuery{
		Chrom: "chr1", Start: 1, End: 512, PageSize: 1000,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 60 {
		t.Fatalf("final records=%d want 60", len(page.Records))
	}
}

// TestRebuild_IntervalQueryHitsCanonicalPosition verifies that after
// normalization moves a variant, interval queries hit it at its canonical
// position and do not return it for intervals that exclude that position.
// This guards against the previous bug where records kept a stale window
// position while their bases were shortened.
func TestRebuild_IntervalQueryHitsCanonicalPosition(t *testing.T) {
	cfg := config.Default()
	ref := fakeRef{seq: refSeq} // CAAAACTACACACG
	root := filepath.Join(t.TempDir(), "db")
	router := newRouter(t, root)

	// Reported at pos 10 (right-aligned spelling), canonical is pos 7.
	data := vcfHeader() + "chr1\t10\t.\tACA\tA\t.\tPASS\t.\tGT\t0/1\n"
	path := filepath.Join(t.TempDir(), "a.vcf")
	writeFile(t, path, data)

	rb, err := New(router, ref, cfg, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := rb.Rebuild(context.Background(),
		[]Source{{Name: "s", Path: path}},
		filepath.Join(t.TempDir(), "cp.json")); err != nil {
		t.Fatal(err)
	}
	st := router.Active()

	// Interval covering the original reported pos 10 but not canonical pos 7
	// must NOT return it.
	page, err := st.Query(store.PageQuery{Chrom: "chr1", Start: 9, End: 14, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 0 {
		t.Fatalf("stale-position hit: got %d records in [9,14]", len(page.Records))
	}

	// Interval covering canonical pos 7 must return exactly it.
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 1, End: 8, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 || page.Records[0].Pos != 7 {
		t.Fatalf("canonical-position miss: got %d records", len(page.Records))
	}

	// A boundary-touching single-position query returns it.
	page, err = st.Query(store.PageQuery{Chrom: "chr1", Start: 7, End: 7, PageSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 {
		t.Fatalf("exact-position query got %d records", len(page.Records))
	}
}

// TestRebuild_MultiSourceMergeAnnotations verifies two sources reporting the
// same variant in different spellings merge into one record, retain both
// sources, and carry stats recomputed from actual merged genotypes.
func TestRebuild_MultiSourceMergeAnnotations(t *testing.T) {
	cfg := config.Default()
	ref := fakeRef{seq: refSeq}
	root := filepath.Join(t.TempDir(), "db")
	router := newRouter(t, root)

	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	writeFile(t, a,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"+
			"chr1\t8\t.\tACA\tA\t60\tPASS\tDP=10\tGT\t0/1\n")
	writeFile(t, b,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n"+
			"chr1\t10\t.\tACA\tA\t80\tPASS\tDP=5\tGT\t0/1\n")

	rb, err := New(router, ref, cfg, importer.Fold(cfg))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := rb.Rebuild(context.Background(),
		[]Source{{Name: "lab-a", Path: a}, {Name: "lab-b", Path: b}},
		filepath.Join(t.TempDir(), "cp.json")); err != nil {
		t.Fatal(err)
	}

	page, err := router.Active().Query(store.PageQuery{
		Chrom: "chr1", Start: 1, End: len(refSeq), PageSize: 100,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(page.Records) != 1 {
		t.Fatalf("merged records=%d want 1", len(page.Records))
	}
	rec := page.Records[0]

	srcs := map[string]bool{}
	for _, s := range rec.Sources {
		srcs[s.Source] = true
	}
	if !srcs["lab-a"] || !srcs["lab-b"] {
		t.Fatalf("sources=%v want lab-a and lab-b retained", rec.Sources)
	}
	if len(rec.Samples) != 2 {
		t.Fatalf("samples=%d want 2 (S1,S2)", len(rec.Samples))
	}
	if rec.Info["AN"] != "4" || rec.Info["AC"] != "2" {
		t.Fatalf("recomputed AN/AC = %q/%q want 4/2", rec.Info["AN"], rec.Info["AC"])
	}
	if !rec.HasQual || rec.Qual != 80 {
		t.Fatalf("QUAL=%v has=%v want max 80", rec.Qual, rec.HasQual)
	}
}
