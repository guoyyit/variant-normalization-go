// varnorm 是 VCF 规范化、多来源合并导入与区间查询的命令行入口。
//
// 子命令：
//
//	import   <store> <source-id> <input.vcf[.gz]|-> [--ref genome.fa] [--config c.json] [--flush N]
//	query    <store> --chrom chr1 [--start N] [--end N] [--limit N] [--after KEY]
//	version
package main

import (
	"context"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/guoyyit/variant-normalization-go/config"
	"github.com/guoyyit/variant-normalization-go/importer"
	"github.com/guoyyit/variant-normalization-go/reference"
	"github.com/guoyyit/variant-normalization-go/store"
)

func main() {
	if len(os.Args) < 2 {
		usage()
		os.Exit(2)
	}
	var err error
	switch os.Args[1] {
	case "import":
		err = runImport(os.Args[2:])
	case "query":
		err = runQuery(os.Args[2:])
	case "version":
		fmt.Println("varnorm dev")
	case "-h", "--help", "help":
		usage()
	default:
		usage()
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

// reorderFlags 把标志参数移动到位置参数之前，兼容 flags 后置的用法。
func reorderFlags(args []string, boolFlags map[string]bool, valued map[string]bool) []string {
	var flags, positional []string
	for i := 0; i < len(args); i++ {
		a := args[i]
		name := strings.TrimLeft(strings.SplitN(a, "=", 2)[0], "-")
		isFlag := strings.HasPrefix(a, "-")
		if !isFlag {
			positional = append(positional, a)
			continue
		}
		flags = append(flags, a)
		if !strings.Contains(a, "=") && valued[name] && i+1 < len(args) {
			i++
			flags = append(flags, args[i])
		}
	}
	return append(flags, positional...)
}

func usage() {
	fmt.Fprintln(os.Stderr, strings.TrimSpace(`
usage:
  varnorm import  <store-dir> <source-id> <input.vcf|-> [--ref genome.fa] [--config c.json] [--flush N]
  varnorm query   <store-dir> --chrom chr1 [--start N] [--end N] [--limit N] [--after KEY]
  varnorm version
`))
}

func runImport(args []string) error {
	fs := flag.NewFlagSet("import", flag.ContinueOnError)
	refPath := fs.String("ref", "", "FASTA reference for left alignment")
	cfgPath := fs.String("config", "", "merge/filter config JSON")
	flush := fs.Int("flush", 10000, "records per in-memory flush")
	if err := fs.Parse(reorderFlags(args, nil, map[string]bool{
		"ref": true, "config": true, "flush": true,
	})); err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 3 {
		return fmt.Errorf("import requires <store-dir> <source-id> <input>")
	}
	storeDir, sourceID, input := rest[0], rest[1], rest[2]

	cfg, err := config.Load(*cfgPath)
	if err != nil {
		return err
	}
	st, err := store.Open(storeDir)
	if err != nil {
		return err
	}
	defer st.Close()

	// 输入必须可 Seek：stdin 先缓冲到临时文件。
	var rs io.ReadSeeker
	var cleanup func()
	if input == "-" {
		tmp, err := os.CreateTemp("", "varnorm-stdin-*.vcf")
		if err != nil {
			return err
		}
		if _, err := io.Copy(tmp, os.Stdin); err != nil {
			return err
		}
		if _, err := tmp.Seek(0, io.SeekStart); err != nil {
			return err
		}
		rs = tmp
		cleanup = func() { _ = tmp.Close(); _ = os.Remove(tmp.Name()) }
	} else {
		f, err := os.Open(input)
		if err != nil {
			return err
		}
		rs = f
		cleanup = func() { _ = f.Close() }
	}
	defer cleanup()

	opts := []importer.Option{importer.WithFlushRecords(*flush)}
	if *refPath != "" {
		genome, err := reference.OpenFASTA(*refPath)
		if err != nil {
			return err
		}
		defer genome.Close()
		opts = append(opts, importer.WithReference(genome))
	}

	job := importer.NewJob(sourceID, rs, st, cfg, opts...)
	prog, err := job.Run(context.Background())
	if err != nil {
		return err
	}
	fmt.Printf("imported source=%s records=%d rejected=%d segments=%d done=%v\n",
		prog.SourceID, prog.Records, prog.Rejected, len(prog.Segments), prog.Done)
	return nil
}

func runQuery(args []string) error {
	fs := flag.NewFlagSet("query", flag.ContinueOnError)
	chrom := fs.String("chrom", "", "chromosome (required)")
	start := fs.Int("start", 0, "start position (1-based, inclusive)")
	end := fs.Int("end", 0, "end position (1-based, inclusive)")
	limit := fs.Int("limit", 0, "page size (0 = all)")
	after := fs.String("after", "", "pagination cursor (Key returned from previous page)")
	if err := fs.Parse(reorderFlags(args, nil, map[string]bool{
		"chrom": true, "start": true, "end": true, "limit": true, "after": true,
	})); err != nil {
		return err
	}
	if fs.NArg() != 1 || *chrom == "" {
		return fmt.Errorf("query requires <store-dir> and --chrom")
	}
	st, err := store.Open(fs.Arg(0))
	if err != nil {
		return err
	}
	defer st.Close()

	sn, err := st.Snapshot(context.Background())
	if err != nil {
		return err
	}
	defer sn.Close()

	got, cursor, err := sn.Query(store.Query{
		Chrom: *chrom, StartPos: *start, EndPos: *end,
		AfterKey: *after, Limit: *limit,
	})
	if err != nil {
		return err
	}
	w := os.Stdout
	for _, c := range got {
		qual := "."
		if c.Record.HasQual {
			qual = fmt.Sprintf("%g", c.Record.Qual)
		}
		fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%s\tQUAL=%s\tFILTER=%s\tsources=%d\t%s\n",
			c.Record.Chrom, c.Record.Pos, strings.Join(c.Record.ID, ","),
			c.Record.Ref, strings.Join(c.Record.Alt, ","),
			qual, strings.Join(c.Record.Filter, ";"),
			len(c.Sources), c.Record.Info)
	}
	if cursor != "" {
		fmt.Fprintf(os.Stderr, "next cursor: %q\n", cursor)
	}
	return nil
}
