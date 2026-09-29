package vcfio

import (
	"bufio"
	"compress/gzip"
	"io"
	"strconv"
	"strings"

	"github.com/guoyyit/variant-normalization-go/variant"
)

// Header 保存 VCF 头部信息，尤其是 FORMAT 字段的 Number 声明。
type Header struct {
	FormatNumber map[string]string
}

// Reader 流式读取 VCF 文本记录。
type Reader struct {
	br        *bufio.Reader
	src       io.Closer
	gz        *gzip.Reader
	header    *Header
	samples   []string
	line      int64
	pending   string
	headerEnd int64
}

// NewReader 创建 VCF Reader（支持明文与 gzip）。
func NewReader(r io.Reader) (*Reader, error) {
	rd := &Reader{header: &Header{FormatNumber: map[string]string{}}}
	br := bufio.NewReader(r)
	peek, err := br.Peek(2)
	if err != nil && err != io.EOF {
		return nil, err
	}
	if len(peek) == 2 && peek[0] == 0x1f && peek[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			return nil, err
		}
		rd.gz = gz
		if c, ok := r.(io.Closer); ok {
			rd.src = c
		}
		br = bufio.NewReader(gz)
	}
	rd.br = br
	if err := rd.readHeader(); err != nil {
		return nil, err
	}
	return rd, nil
}

// readHeader 预读并解析所有以 # 开头的头部行；遇到第一条数据行则暂存。
func (r *Reader) readHeader() error {
	for {
		raw, err := r.br.ReadString('\n')
		r.headerEnd += int64(len(raw))
		if err != nil && raw == "" {
			return nil
		}
		r.line++
		line := strings.TrimRight(raw, "\r\n")
		if !strings.HasPrefix(line, "#") {
			r.pending = line
			return nil
		}
		if err := r.readHeaderLine(line); err != nil {
			return err
		}
		if err != nil {
			return err
		}
	}
}

// HeaderEnd 返回头部（含）结束、第一条数据行开始处的字节偏移。
func (r *Reader) HeaderEnd() int64 { return r.headerEnd }

// Header 返回已解析的头部（读取记录前可用）。
func (r *Reader) Header() *Header { return r.header }

// Samples 返回头部声明的样本名。
func (r *Reader) Samples() []string { return r.samples }

func (r *Reader) readHeaderLine(line string) error {
	if strings.HasPrefix(line, "##FORMAT=<") {
		id, number, ok := parseHeaderMeta(line)
		if ok {
			r.header.FormatNumber[id] = number
		}
		return nil
	}
	if strings.HasPrefix(line, "#CHROM") {
		f := strings.Split(line, "\t")
		if len(f) > 9 {
			r.samples = append(r.samples, f[9:]...)
		}
	}
	return nil
}

// parseHeaderMeta 从 ##FORMAT=<ID=GT,Number=1,...> 提取 ID 与 Number。
func parseHeaderMeta(line string) (id, number string, ok bool) {
	s := strings.TrimPrefix(line, "##FORMAT=")
	s = strings.TrimPrefix(s, "<")
	s = strings.TrimSuffix(s, ">")
	var key, val string
	emit := func(k, v string) {
		switch strings.TrimSpace(k) {
		case "ID":
			id = strings.TrimSpace(v)
		case "Number":
			number = strings.Trim(strings.TrimSpace(v), "\"")
		}
	}
	inQuote := false
	start := 0
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case '"':
			inQuote = !inQuote
		case ',':
			if inQuote {
				continue
			}
			key, val, _ = strings.Cut(s[start:i], "=")
			emit(key, val)
			start = i + 1
		}
	}
	key, val, _ = strings.Cut(s[start:], "=")
	emit(key, val)
	return id, number, id != ""
}

// Read 读取下一条记录，流结束返回 io.EOF。
func (r *Reader) Read() (*variant.Record, error) {
	for {
		if r.pending != "" {
			line := r.pending
			r.pending = ""
			rec, perr := parseRecord(line, r.samples)
			if perr == nil {
				rec.Source.Line = r.line
			}
			return rec, perr
		}
		line, err := r.br.ReadString('\n')
		if err != nil && line == "" {
			return nil, err
		}
		r.line++
		line = strings.TrimRight(line, "\r\n")
		if strings.HasPrefix(line, "##") {
			if err := r.readHeaderLine(line); err != nil {
				return nil, err
			}
		} else if strings.HasPrefix(line, "#") {
			if err := r.readHeaderLine(line); err != nil {
				return nil, err
			}
		} else if line != "" {
			rec, perr := parseRecord(line, r.samples)
			if perr == nil {
				rec.Source.Line = r.line
			}
			return rec, perr
		}
		if err != nil {
			return nil, err
		}
	}
}

// Close 释放 gzip 与底层资源。
func (r *Reader) Close() error {
	var first error
	if r.gz != nil {
		first = r.gz.Close()
	}
	if r.src != nil {
		if err := r.src.Close(); err != nil && first == nil {
			first = err
		}
	}
	return first
}

func parseRecord(line string, sampleNames []string) (*variant.Record, error) {
	f := strings.Split(line, "\t")
	if len(f) < 8 {
		return nil, &variant.ValidationError{Kind: variant.KindInvalidFormat,
			Msg: "expected at least 8 tab-separated columns, got " + strconv.Itoa(len(f))}
	}
	pos, err := strconv.Atoi(f[1])
	if err != nil {
		return nil, &variant.ValidationError{Kind: variant.KindInvalidFormat, Msg: "invalid POS: " + f[1]}
	}
	rec := &variant.Record{
		Chrom:  f[0],
		Pos:    pos,
		ID:     splitList(f[2]),
		Ref:    f[3],
		Alt:    splitList(f[4]),
		Filter: splitList(f[6]),
		Info:   parseInfo(f[7]),
	}
	if f[5] != "." {
		q, err := strconv.ParseFloat(f[5], 64)
		if err != nil {
			return nil, &variant.ValidationError{Kind: variant.KindInvalidFormat, Msg: "invalid QUAL: " + f[5]}
		}
		rec.Qual, rec.HasQual = q, true
	}
	if len(f) >= 9 && f[8] != "." {
		rec.Format = strings.Split(f[8], ":")
		for i := 9; i < len(f); i++ {
			name := strconv.Itoa(i - 8)
			if i-9 < len(sampleNames) {
				name = sampleNames[i-9]
			}
			rec.Samples = append(rec.Samples, parseSample(name, rec.Format, f[i]))
		}
	}
	return rec, nil
}

func splitList(s string) []string {
	if s == "." || s == "" {
		return nil
	}
	return strings.Split(s, ",")
}

func parseInfo(s string) map[string]string {
	info := map[string]string{}
	if s == "." || s == "" {
		return info
	}
	for _, kv := range strings.Split(s, ";") {
		if eq := strings.IndexByte(kv, '='); eq >= 0 {
			info[kv[:eq]] = kv[eq+1:]
		} else {
			info[kv] = ""
		}
	}
	return info
}

func parseSample(name string, format []string, raw string) variant.Sample {
	s := variant.Sample{Name: name, Fields: map[string]string{}}
	parts := strings.Split(raw, ":")
	for i, key := range format {
		if i >= len(parts) {
			break
		}
		s.Fields[key] = parts[i]
		if key == "GT" {
			s.GT = parseGT(parts[i])
		}
	}
	return s
}

// parseGT 解析如 0/1、1|0、./1、2 的基因型。
func parseGT(s string) variant.Genotype {
	g := variant.Genotype{}
	if strings.Contains(s, "|") {
		g.Phased = true
	}
	s = strings.ReplaceAll(s, "|", "/")
	for _, a := range strings.Split(s, "/") {
		if a == "." {
			g.Alleles = append(g.Alleles, -1)
			continue
		}
		v, err := strconv.Atoi(a)
		if err != nil {
			g.Alleles = append(g.Alleles, -1)
			continue
		}
		g.Alleles = append(g.Alleles, v)
	}
	return g
}

// DataReader 在已知头部后、从数据区读取记录（用于断点恢复时跳过头部解析）。
// 调用方需已提供解析好的 Header 与样本名，并把 r 定位到数据起始偏移。
type DataReader struct {
	br      *bufio.Reader
	header  *Header
	samples []string
	line    int64
}

// NewDataReader 创建数据读取器。r 必须已 Seek 到第一条数据行。
func NewDataReader(r io.Reader, header *Header, samples []string) *DataReader {
	return &DataReader{br: bufio.NewReader(r), header: header, samples: samples}
}

// Header 返回共享头部。
func (d *DataReader) Header() *Header { return d.header }

// Read 读取下一条数据记录。
func (d *DataReader) Read() (*variant.Record, error) {
	for {
		line, err := d.br.ReadString('\n')
		if err != nil && line == "" {
			return nil, err
		}
		d.line++
		line = strings.TrimRight(line, "\r\n")
		if line == "" {
			if err != nil {
				return nil, err
			}
			continue
		}
		rec, perr := parseRecord(line, d.samples)
		if perr == nil {
			rec.Source.Line = d.line
		}
		return rec, perr
	}
}

// ParseRecord 解析一条不含换行的 VCF 数据行（供需要自行管理字节偏移的调用方使用）。
func ParseRecord(line string, samples []string) (*variant.Record, error) {
	return parseRecord(line, samples)
}

// NewEmptyHeader 创建空头部（供自行逐行管理偏移的调用方累积元信息）。
func NewEmptyHeader() *Header { return &Header{FormatNumber: map[string]string{}} }

// AddHeaderLine 处理一条 ## 或 #CHROM 头部行；返回样本名列（仅 #CHROM 行非空）。
func (h *Header) AddHeaderLine(line string) []string {
	if strings.HasPrefix(line, "##FORMAT=<") {
		if id, num, ok := parseHeaderMeta(line); ok {
			h.FormatNumber[id] = num
		}
		return nil
	}
	if strings.HasPrefix(line, "#CHROM") {
		f := strings.Split(line, "\t")
		if len(f) > 9 {
			return f[9:]
		}
	}
	return nil
}
