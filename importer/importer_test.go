package importer

import (
	"context"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/guoyyit/variant-normalization-go/config"
	"github.com/guoyyit/variant-normalization-go/reference"
	"github.com/guoyyit/variant-normalization-go/store"
)

const header = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`

func seekFile(t *testing.T, content string) (*os.File, func()) {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	return f, func() { _ = f.Close() }
}

func genome() *reference.Memory {
	return reference.NewMemory(map[string]string{
		"chr1": strings.Repeat("A", 2048),
	})
}

func TestImporter_BasicDedup(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(dir)
	if err != nil {
		t.Fatal(err)
	}
	header2 := strings.Replace(header,
		"FORMAT\tS1\n", "FORMAT\tS1\tS2\n", 1)
	body := header2 +
		"chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n"
	f, cleanup := seekFile(t, body)
	defer cleanup()

	job := NewJob("src1", f, st, config.Default(), WithReference(genome()), WithFlushRecords(100))
	prog, err := job.Run(context.Background())
	if err != nil {
		t.Fatal(err)
	}
	if prog.Done != true || prog.Records != 1 {
		t.Fatalf("progress = %+v", prog)
	}

	sn, _ := st.Snapshot(context.Background())
	defer sn.Close()
	got, _, err := sn.Query(store.Query{Chrom: "chr1"})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 1 {
		t.Fatalf("expected dedup to 1 record, got %d", len(got))
	}
	// 统计应基于合并后样本重算
	if got[0].Record.Info["AC"] != "3" || got[0].Record.Info["AN"] != "4" {
		t.Fatalf("stats = %+v", got[0].Record.Info)
	}
	if len(got[0].Sources) != 1 {
		t.Fatalf("sources = %d", len(got[0].Sources))
	}
}

func TestImporter_MultiAlleleSplit(t *testing.T) {
	dir := t.TempDir()
	st, _ := store.Open(dir)
	body := header + "chr1\t10\t.\tA\tT,C\t.\tPASS\t.\tGT\t1/2\n"
	f, cleanup := seekFile(t, body)
	defer cleanup()

	job := NewJob("s", f, st, config.Default(), WithReference(genome()))
	if _, err := job.Run(context.Background()); err != nil {
		t.Fatal(err)
	}
	sn, _ := st.Snapshot(context.Background())
	defer sn.Close()
	got, _, _ := sn.Query(store.Query{Chrom: "chr1"})
	if len(got) != 1 || len(got[0].Record.Alt) != 2 {
		t.Fatalf("expected 1 merged record with 2 ALTs, got %d records", len(got))
	}
}

func TestImporter_RejectInvalid(t *testing.T) {
	dir := t.TempDir()
	st, _ := store.Open(dir)
	body := header +
		"chr1\t10\t.\tA\tZ\t.\tPASS\t.\tGT\t0/1\n" + // 非法碱基
		"chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
	f, cleanup := seekFile(t, body)
	defer cleanup()
	job := NewJob("s", f, st, config.Default(), WithReference(genome()))
	prog, err := job.Run(context.Background())
	if err != nil {
		t.Fatal(err)
	}
	if prog.Rejected != 1 || prog.Records != 1 {
		t.Fatalf("rejected=%d records=%d", prog.Rejected, prog.Records)
	}
}

func TestImporter_RecoveryNoDupNoLoss(t *testing.T) {
	dir := t.TempDir()
	st, _ := store.Open(dir)

	// 构造大量记录，强制小缓冲多次刷盘；中途取消再恢复。
	var sb strings.Builder
	sb.WriteString(header)
	for i := 0; i < 250; i++ {
		pos := 10 + i
		// REF 必须与参考一致；该参考 pos10 以后多为 A，用 A->T SNP。
		if pos <= genome().Length("chr1") {
			sb.WriteString("chr1\t")
			sb.WriteString(itoa(pos))
			sb.WriteString("\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n")
		}
	}
	f, cleanup := seekFile(t, sb.String())
	defer cleanup()

	// 第一次：极小缓冲，但在读到中段时通过 context 取消。
	ctx, cancel := context.WithCancel(context.Background())
	job := NewJob("s", f, st, config.Default(),
		WithReference(genome()), WithFlushRecords(20))
	go func() {
		time.Sleep(20 * time.Millisecond)
		cancel()
	}()
	p1, err := job.Run(ctx)
	if err == nil {
		t.Fatal("expected cancellation error")
	}
	if p1.Done {
		t.Fatal("should not be done after cancel")
	}

	// 恢复：重新打开同一文件（offset 已在检查点中），正常跑完。
	f2, err := os.Open(f.Name())
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	job2 := NewJob("s", f2, st, config.Default(),
		WithReference(genome()), WithFlushRecords(20))
	p2, err := job2.Run(context.Background())
	if err != nil {
		t.Fatalf("resume failed: %v", err)
	}
	if !p2.Done {
		t.Fatal("resume should complete")
	}

	sn, _ := st.Snapshot(context.Background())
	defer sn.Close()
	var all int
	var cursor string
	for {
		got, next, qerr := sn.Query(store.Query{Chrom: "chr1", Limit: 64, AfterKey: cursor})
		if qerr != nil {
			t.Fatal(qerr)
		}
		all += len(got)
		if next == "" {
			break
		}
		cursor = next
	}
	// pos10..259 共 250 条，全部是 A 的位置才合法；统计实际通过数。
	// 参考中 pos10..259 均为 A（见 genome），因此应为 250。
	if all != 250 {
		t.Fatalf("expected 250 unique records after resume, got %d", all)
	}
}

func TestImporter_ConcurrentQueriesNeverSeePartial(t *testing.T) {
	dir := t.TempDir()
	st, _ := store.Open(dir)
	var sb strings.Builder
	sb.WriteString(header)
	for i := 0; i < 300; i++ {
		sb.WriteString("chr1\t")
		sb.WriteString(itoa(10 + i))
		sb.WriteString("\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n")
	}
	f, cleanup := seekFile(t, sb.String())
	defer cleanup()

	var wg sync.WaitGroup
	stop := make(chan struct{})
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			select {
			case <-stop:
				return
			default:
				sn, err := st.Snapshot(context.Background())
				if err == nil {
					got, _, _ := sn.Query(store.Query{Chrom: "chr1"})
					// 在提交完成前，查询只能看到 0 条；提交后看到 300 条。
					if n := len(got); n != 0 && n != 300 {
						t.Errorf("observed intermediate count %d", n)
					}
					_ = sn.Close()
				}
			}
		}
	}()

	job := NewJob("s", f, st, config.Default(),
		WithReference(genome()), WithFlushRecords(25))
	if _, err := job.Run(context.Background()); err != nil {
		t.Fatal(err)
	}
	close(stop)
	wg.Wait()
}

func TestImporter_ConfigOnlyAffectsFuture(t *testing.T) {
	dir := t.TempDir()
	st, _ := store.Open(dir)
	makeBody := func(qual int) string {
		return header + "chr1\t10\t.\tA\tT\t" + itoa(qual) + "\tPASS\t.\tGT\t0/1\n"
	}

	// 第一批：无阈值
	f1, c1 := seekFile(t, makeBody(10))
	j1 := NewJob("a", f1, st, &config.Config{ConfigVersion: "v1", DefaultPriority: 100},
		WithReference(genome()))
	if _, err := j1.Run(context.Background()); err != nil {
		t.Fatal(err)
	}
	c1()

	sn1, _ := st.Snapshot(context.Background())
	g1, _, _ := sn1.Query(store.Query{Chrom: "chr1"})
	if g1[0].ConfigVersion != "v1" || g1[0].Record.Filter[0] != "PASS" {
		t.Fatalf("v1 record = %+v filter=%v", g1[0].ConfigVersion, g1[0].Record.Filter)
	}
	sn1.Close()

	// 第二批：高阈值，只影响新导入；已有 v1 记录保持自洽。
	f2, c2 := seekFile(t, makeBody(5))
	cfg2 := &config.Config{ConfigVersion: "v2", DefaultPriority: 100, MinQual: 30}
	j2 := NewJob("b", f2, st, cfg2, WithReference(genome()))
	if _, err := j2.Run(context.Background()); err != nil {
		t.Fatal(err)
	}
	c2()

	sn2, _ := st.Snapshot(context.Background())
	defer sn2.Close()
	g2, _, _ := sn2.Query(store.Query{Chrom: "chr1"})
	// 同一变异合并后 QUAL 取 max=10，仍低于 30，被标记 LowQual（新配置生效于合并点）。
	var sawLowQual bool
	for _, c := range g2 {
		for _, fl := range c.Record.Filter {
			if fl == "LowQual" {
				sawLowQual = true
			}
		}
	}
	if !sawLowQual {
		t.Fatalf("expected LowQual under new config: %+v", g2[0].Record.Filter)
	}
}

func itoa(n int) string { return strconv.Itoa(n) }
