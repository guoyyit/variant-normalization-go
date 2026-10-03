package main

import (
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"os/signal"
	"path/filepath"
	"syscall"

	"github.com/guoyyit/variant-normalization-go/pkg/config"
	"github.com/guoyyit/variant-normalization-go/pkg/genome"
	"github.com/guoyyit/variant-normalization-go/pkg/importer"
	"github.com/guoyyit/variant-normalization-go/pkg/normalize"
	"github.com/guoyyit/variant-normalization-go/pkg/stats"
	"github.com/guoyyit/variant-normalization-go/pkg/store"
)

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "varnorm:", err)
		os.Exit(1)
	}
}

func run(args []string) error {
	if len(args) == 0 {
		usage()
		return fmt.Errorf("missing subcommand")
	}
	switch args[0] {
	case "import":
		return importCmd(args[1:], false)
	case "resume":
		return importCmd(args[1:], true)
	case "query":
		return queryCmd(args[1:])
	case "rebuild":
		return rebuildCmd(args[1:])
	case "-h", "--help", "help":
		usage()
		return nil
	default:
		usage()
		return fmt.Errorf("unknown subcommand %q", args[0])
	}
}

func commonFlags(fs *flag.FlagSet) (string, string, *string, *string) {
	db := fs.String("db", "varnorm.db", "store directory")
	cfgPath := fs.String("config", "", "JSON configuration file")
	refPath := fs.String("reference", "", "FASTA reference genome")
	cp := fs.String("checkpoint", "checkpoint.json", "checkpoint file")
	return *db, *cfgPath, refPath, cp
}

func loadConfig(path string) (config.Config, error) {
	if path == "" {
		return config.Default(), nil
	}
	return config.LoadFile(path)
}

func importCmd(args []string, resume bool) error {
	fs := flag.NewFlagSet("import", flag.ContinueOnError)
	dbDir, cfgPath, refPath, cpPath := commonFlags(fs)
	source := fs.String("source", "", "logical source name (required)")
	fs.String("file", "", "kept for CLI symmetry; positional path is used")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 || *source == "" {
		return fmt.Errorf("usage: varnorm import -source NAME [-reference REF.fa] [-config C.json] [-db DIR] FILE.vcf")
	}

	cfg, err := loadConfig(cfgPath)
	if err != nil {
		return err
	}

	var ref interface {
		Substring(string, int, int) (string, error)
		Length(string) (int, error)
	}
	if *refPath != "" {
		g, err := genome.LoadFASTA(*refPath)
		if err != nil {
			return err
		}
		ref = g
	}

	st, err := store.Open(dbDir, importer.Fold(cfg))
	if err != nil {
		return err
	}
	defer st.Close()

	im, err := importer.New(st, ref, cfg, *cpPath)
	if err != nil {
		return err
	}

	ctx, cancel := signalContext()
	defer cancel()

	fileName := filepath.Base(fs.Arg(0))
	f, err := os.Open(fs.Arg(0))
	if err != nil {
		return err
	}
	defer f.Close()

	var report importer.Report
	if resume {
		report, err = im.Resume(ctx, *source, fileName, f)
	} else {
		report, err = im.Import(ctx, *source, fileName, f)
	}
	if err != nil {
		return err
	}

	fmt.Printf("source=%s lines=%d imported=%d rejected=%d\n",
		report.Source, report.Lines, report.Imported, len(report.Rejected))
	for _, rej := range report.Rejected {
		fmt.Fprintf(os.Stderr, "  rejected line %d: %s\n", rej.Line, rej.Reason)
	}
	return nil
}

func queryCmd(args []string) error {
	fs := flag.NewFlagSet("query", flag.ContinueOnError)
	dbDir, cfgPath, _, _ := commonFlags(fs)
	chrom := fs.String("chrom", "", "chromosome (required)")
	start := fs.Int("start", 1, "interval start (1-based, inclusive)")
	end := fs.Int("end", 0, "interval end (1-based, inclusive)")
	pageSize := fs.Int("page-size", 50, "records per page")
	token := fs.String("page-token", "", "opaque continuation token")
	jsonOut := fs.Bool("json", false, "emit JSON instead of text lines")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *chrom == "" || *end == 0 {
		return fmt.Errorf("usage: varnorm query -chrom C -start N -end M [-page-size K] [-page-token T]")
	}

	cfg, err := loadConfig(cfgPath)
	if err != nil {
		return err
	}
	st, err := store.Open(dbDir, importer.Fold(cfg))
	if err != nil {
		return err
	}
	defer st.Close()

	page, err := st.Query(store.PageQuery{
		Chrom:     *chrom,
		Start:     *start,
		End:       *end,
		PageSize:  *pageSize,
		PageToken: *token,
	})
	if err != nil {
		return err
	}

	if *jsonOut {
		enc := json.NewEncoder(os.Stdout)
		return enc.Encode(page)
	}
	for _, rec := range page.Records {
		fmt.Println(stats.Summary(rec))
	}
	if page.NextToken != "" {
		fmt.Fprintf(os.Stderr, "next-token: %s\n", page.NextToken)
	}
	return nil
}

func signalContext() (context.Context, context.CancelFunc) {
	ctx, cancel := context.WithCancel(context.Background())
	ch := make(chan os.Signal, 1)
	signal.Notify(ch, syscall.SIGINT, syscall.SIGTERM)
	go func() {
		<-ch
		cancel()
	}()
	return ctx, cancel
}

func usage() {
	fmt.Fprint(os.Stderr, `varnorm - VCF normalization & multi-source merging

usage:
  varnorm import -source NAME [-reference REF.fa] [-config C.json] [-db DIR]
                 [-checkpoint cp.json] FILE.vcf
  varnorm resume  (same flags; continues from the checkpoint)
	varnorm query  -chrom C -start N -end M [-page-size K] [-page-token T]
                 [-json] [-db DIR]
  varnorm rebuild -source NAME [-reference REF.fa] [-config C.json]
                 [-db DIR] [-checkpoint cp.json] FILE.vcf [FILE.vcf ...]
`)
}

func rebuildCmd(args []string) error {
	fs := flag.NewFlagSet("rebuild", flag.ContinueOnError)
	dbDir, cfgPath, refPath, cpPath := commonFlags(fs)
	source := fs.String("source", "", "logical source name (required)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() < 1 || *source == "" {
		return fmt.Errorf("usage: varnorm rebuild -source NAME [-reference REF.fa] [-config C.json] [-db DIR] [-checkpoint cp.json] FILE.vcf [FILE.vcf ...]")
	}

	cfg, err := loadConfig(cfgPath)
	if err != nil {
		return err
	}

	var ref normalize.Reference
	if *refPath != "" {
		g, err := genome.LoadFASTA(*refPath)
		if err != nil {
			return err
		}
		ref = g
	}

	files := make([]importer.SourceFile, 0, fs.NArg())
	for _, arg := range fs.Args() {
		files = append(files, importer.SourceFile{
			Source: *source,
			File:   filepath.Base(arg),
			Path:   arg,
		})
	}

	ctx, cancel := signalContext()
	defer cancel()

	rep, err := importer.Rebuild(ctx, dbDir, ref, cfg, *cpPath, files)
	if err != nil {
		return err
	}
	fmt.Printf("rebuilt source=%s files=%d imported=%d rejected=%d\n",
		*source, len(files), rep.Imported, len(rep.Rejected))
	for _, rej := range rep.Rejected {
		fmt.Fprintf(os.Stderr, "  rejected line %d: %s\n", rej.Line, rej.Reason)
	}
	return nil
}
