// Package importer 实现流式批量导入与断点恢复。
//
// 内存安全：记录在内存中缓冲到固定条数（FlushRecords）后排序并刷成一个
// 临时有序段，因此常驻内存只与缓冲大小有关，不随输入总量线性增长。
//
// 断点恢复：每刷完一个段就在 jobs/<source>.json 中原子记录下一条待处理
// 数据行的绝对字节偏移、累计记录数与已完成段列表。重启时 Seek 到检查点
// 偏移继续，已刷盘段不重复处理，从而不重复、不遗漏。
package importer

import (
	"bufio"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"github.com/guoyyit/variant-normalization-go/config"
	"github.com/guoyyit/variant-normalization-go/merge"
	"github.com/guoyyit/variant-normalization-go/normalize"
	"github.com/guoyyit/variant-normalization-go/reference"
	"github.com/guoyyit/variant-normalization-go/store"
	"github.com/guoyyit/variant-normalization-go/variant"
	"github.com/guoyyit/variant-normalization-go/vcfio"
)

// 默认缓冲条数：控制单次导入的峰值内存。
const defaultFlushRecords = 10000

// Progress 是导入进度的持久化检查点信息。
type Progress struct {
	SourceID     string   `json:"source_id"`
	Bytes        int64    `json:"bytes"`
	Records      int64    `json:"records"`
	Rejected     int64    `json:"rejected"`
	Segments     []string `json:"segments"`
	NextLine     int64    `json:"next_line"`
	Done         bool     `json:"done"`
	FlushRecords int      `json:"flush_records"`
}

// Job 为一次可断点恢复的流式批量导入。
type Job struct {
	sourceID string
	reader   io.ReadSeeker
	dir      string
	st       *store.Store
	cfg      *config.Config
	genome   reference.Genome
	flush    int

	mu       sync.Mutex
	progress Progress
}

// Option 配置导入任务。
type Option func(*Job)

// WithReference 注入用于规范化左对齐的参考序列。
func WithReference(g reference.Genome) Option {
	return func(j *Job) { j.genome = g }
}

// WithFlushRecords 设置内存缓冲条数（峰值内存控制）。
func WithFlushRecords(n int) Option {
	return func(j *Job) {
		if n > 0 {
			j.flush = n
		}
	}
}

// NewJob 创建导入任务。reader 必须支持 Seek（文件或缓冲到磁盘的输入）。
func NewJob(sourceID string, reader io.ReadSeeker, st *store.Store, cfg *config.Config, opts ...Option) *Job {
	if cfg == nil {
		cfg = config.Default()
	}
	j := &Job{
		sourceID: sourceID,
		reader:   reader,
		dir:      filepath.Join(st.Dir(), "jobs"),
		st:       st,
		cfg:      cfg,
		flush:    defaultFlushRecords,
	}
	for _, o := range opts {
		o(j)
	}
	return j
}

func (j *Job) checkpointPath() string {
	return filepath.Join(j.dir, j.sourceID+".json")
}

func (j *Job) loadProgress() (Progress, error) {
	var p Progress
	data, err := os.ReadFile(j.checkpointPath())
	if errors.Is(err, os.ErrNotExist) {
		return Progress{SourceID: j.sourceID, FlushRecords: j.flush}, nil
	}
	if err != nil {
		return p, err
	}
	if err := json.Unmarshal(data, &p); err != nil {
		return p, fmt.Errorf("decode checkpoint: %w", err)
	}
	if p.FlushRecords <= 0 {
		p.FlushRecords = j.flush
	}
	return p, nil
}

func (j *Job) saveProgress(p Progress) error {
	if err := os.MkdirAll(j.dir, 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(j.dir, ".ckpt-*")
	if err != nil {
		return err
	}
	name := tmp.Name()
	enc := json.NewEncoder(tmp)
	enc.SetIndent("", "  ")
	if err := enc.Encode(p); err != nil {
		_ = tmp.Close()
		_ = os.Remove(name)
		return err
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		_ = os.Remove(name)
		return err
	}
	if err := tmp.Close(); err != nil {
		_ = os.Remove(name)
		return err
	}
	return os.Rename(name, j.checkpointPath())
}

// readHeader 从文件起始逐行解析头部。返回：头部、样本名、第一条数据行原文
// （可能为空，表示恰好在行边界结束）以及数据区起始绝对偏移。
func (j *Job) readHeader() (*vcfio.Header, []string, string, int64, error) {
	if _, err := j.reader.Seek(0, io.SeekStart); err != nil {
		return nil, nil, "", 0, err
	}
	br := bufio.NewReader(j.reader)
	header := vcfio.NewEmptyHeader()
	var samples []string
	var off int64
	for {
		raw, err := br.ReadString('\n')
		off += int64(len(raw))
		line := strings.TrimRight(raw, "\r\n")
		if s := header.AddHeaderLine(line); s != nil {
			samples = s
		}
		if !strings.HasPrefix(line, "#") {
			return header, samples, line, off, nil
		}
		if err != nil {
			if errors.Is(err, io.EOF) {
				return header, samples, "", off, nil
			}
			return nil, nil, "", 0, err
		}
	}
}

func (j *Job) reject(prog *Progress) {
	prog.Rejected++
	j.progress.Rejected = prog.Rejected
}

// Run 执行导入；失败后再次调用从最近检查点继续，不重复不遗漏。
func (j *Job) Run(ctx context.Context) (*Progress, error) {
	j.mu.Lock()
	defer j.mu.Unlock()

	prog, err := j.loadProgress()
	if err != nil {
		return nil, err
	}
	if prog.Done {
		return &prog, nil
	}
	j.progress = prog

	existing := append([]string(nil), prog.Segments...)

	header, samples, firstDataLine, dataStart, err := j.readHeader()
	if err != nil {
		return nil, err
	}
	if prog.Bytes == 0 {
		prog.Bytes = dataStart
	}

	if _, err := j.reader.Seek(prog.Bytes, io.SeekStart); err != nil {
		return nil, fmt.Errorf("seek checkpoint: %w", err)
	}
	br := bufio.NewReader(j.reader)
	off := prog.Bytes
	lineNo := prog.NextLine
	// 首次导入时第一条数据行已在解析头部时被读出，需要补处理一次。
	pendingFirst := prog.Bytes <= dataStart && firstDataLine != ""
	if pendingFirst {
		lineNo = 0
	}

	pipe := newPipeline(normalize.New(j.genome))
	mg := merge.NewMerger(j.cfg)

	buf := make([]*variant.Canonical, 0, prog.FlushRecords)
	flush := func() error {
		if len(buf) == 0 {
			return nil
		}
		merged := dedupMerge(mg, buf)
		sortCanonical(merged)
		segName := fmt.Sprintf("incoming-%s-%d.gob", j.sourceID, len(existing))
		w, err := j.st.NewSegmentWriter(segName)
		if err != nil {
			return err
		}
		for _, c := range merged {
			if err := w.Write(c); err != nil {
				return err
			}
		}
		path, err := w.Finish(j.st, segName)
		if err != nil {
			return err
		}
		existing = append(existing, filepath.Base(path))
		prog.Segments = existing
		prog.Bytes = off
		prog.NextLine = lineNo
		if err := j.saveProgress(prog); err != nil {
			return err
		}
		buf = buf[:0]
		return nil
	}

	handleLine := func(line string) error {
		if strings.TrimSpace(line) == "" {
			return nil
		}
		rec, perr := vcfio.ParseRecord(line, samples)
		if perr != nil {
			if asValidation(perr) {
				j.reject(&prog)
				return nil
			}
			return perr
		}
		rec.Source.ID = j.sourceID
		lineNo++
		rec.Source.Line = lineNo
		produced, perr := pipe.process(rec, header)
		if perr != nil {
			if asValidation(perr) {
				j.reject(&prog)
				return nil
			}
			return perr
		}
		buf = append(buf, produced...)
		prog.Records += int64(len(produced))
		j.progress.Records = prog.Records
		return nil
	}

	// 首次导入：头部阶段已读出的第一条数据行先处理（其字节已计入 dataStart）。
	if pendingFirst {
		if err := handleLine(firstDataLine); err != nil {
			return &j.progress, err
		}
		if len(buf) >= prog.FlushRecords {
			if err := flush(); err != nil {
				return &j.progress, err
			}
		}
	}

	for {
		if err := ctx.Err(); err != nil {
			return &j.progress, err
		}
		raw, err := br.ReadString('\n')
		if len(raw) > 0 {
			off += int64(len(raw))
			if perr := handleLine(strings.TrimRight(raw, "\r\n")); perr != nil {
				return &j.progress, perr
			}
		}
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return &j.progress, err
		}
		if len(buf) >= prog.FlushRecords {
			if err := flush(); err != nil {
				return &j.progress, err
			}
		}
	}
	if err := flush(); err != nil {
		return &j.progress, err
	}

	finalSeg, err := j.compact(mg, existing)
	if err != nil {
		return &j.progress, err
	}
	if err := j.st.Commit([]string{finalSeg}); err != nil {
		return &j.progress, err
	}

	prog.Done = true
	j.progress = prog
	if err := j.saveProgress(prog); err != nil {
		return &j.progress, err
	}
	return &j.progress, nil
}

func asValidation(err error) bool {
	var ve *variant.ValidationError
	return errors.As(err, &ve)
}

// compact 把当前已提交版本的段与本次导入段做 k 路归并、去重合并，
// 输出一个全新有序段（跨来源合并在提交点一次完成，查询看不到中间态）。
func (j *Job) compact(mg *merge.Merger, incoming []string) (string, error) {
	all := append([]string(nil), j.st.CurrentSegmentPaths()...)
	for _, name := range incoming {
		all = append(all, j.st.SegmentPath(name))
	}

	es, err := openEntries(all)
	if err != nil {
		return "", err
	}
	defer es.close()

	outName := fmt.Sprintf("committed-%d.gob", j.st.Generation()+1)
	w, err := j.st.NewSegmentWriter(outName)
	if err != nil {
		return "", err
	}

	var group []*variant.Canonical
	groupKey := func(c *variant.Canonical) string { return merge.GroupKey(&c.Record) }
	flushGroup := func() error {
		if len(group) == 0 {
			return nil
		}
		if err := w.Write(mg.Merge(group)); err != nil {
			return err
		}
		group = group[:0]
		return nil
	}

	var iterErr error
	for c, err := range es.channel() {
		if err != nil {
			iterErr = err
			break
		}
		if len(group) > 0 && groupKey(c) != groupKey(group[0]) {
			if err := flushGroup(); err != nil {
				return "", err
			}
		}
		group = append(group, c)
	}
	if iterErr != nil {
		return "", iterErr
	}
	if err := flushGroup(); err != nil {
		return "", err
	}
	path, err := w.Finish(j.st, outName)
	if err != nil {
		return "", err
	}
	return filepath.Base(path), nil
}

// Progress 返回当前持久化进度。
func (j *Job) Progress() (*Progress, error) {
	j.mu.Lock()
	defer j.mu.Unlock()
	p, err := j.loadProgress()
	if err != nil {
		return nil, err
	}
	return &p, nil
}

func sortCanonical(in []*variant.Canonical) {
	sort.Slice(in, func(i, j int) bool {
		return in[i].Record.Key() < in[j].Record.Key()
	})
}
