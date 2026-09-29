package reference

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
)

var (
	// ErrUnknownChrom 表示参考中不存在该染色体。
	ErrUnknownChrom = errors.New("unknown chromosome")
	// ErrOutOfBounds 表示请求的坐标超出染色体范围。
	ErrOutOfBounds = errors.New("coordinates out of bounds")
)

// Genome 提供规范化左对齐所需的参考序列随机读取能力。
type Genome interface {
	// Prefix 返回染色体 1-based 区间 [end-len+1, end] 的碱基。
	Prefix(chrom string, end, length int) (string, error)
	// Length 返回染色体长度，坐标为 0 表示未知。
	Length(chrom string) int
}

// Memory 是基于内存字符串的参考实现，主要用于测试。
type Memory struct {
	seq map[string]string
}

// NewMemory 从 染色体 -> 大写碱基序列 的映射创建内存参考。
func NewMemory(seq map[string]string) *Memory {
	cp := make(map[string]string, len(seq))
	for k, v := range seq {
		cp[k] = strings.ToUpper(strings.TrimSpace(v))
	}
	return &Memory{seq: cp}
}

func (m *Memory) Prefix(chrom string, end, length int) (string, error) {
	s, ok := m.seq[chrom]
	if !ok {
		return "", fmt.Errorf("%w: %s", ErrUnknownChrom, chrom)
	}
	if length < 0 || end < length || end > len(s) {
		return "", fmt.Errorf("%w: %s:%d length=%d", ErrOutOfBounds, chrom, end, length)
	}
	return s[end-length : end], nil
}

func (m *Memory) Length(chrom string) int {
	if s, ok := m.seq[chrom]; ok {
		return len(s)
	}
	return 0
}

type fastaIndex struct {
	length    int64 // 序列碱基数
	offset    int64 // 第一条序列行在文件中的字节偏移
	lineBases int64 // 每条完整序列行的碱基数
	lineBytes int64 // 每条完整序列行的字节数（含换行）
}

// FASTA 是支持随机读取的 FASTA 参考（要求每条序列等宽换行，与 samtools faidx 一致）。
type FASTA struct {
	f   *os.File
	idx map[string]fastaIndex
}

// OpenFASTA 打开并扫描索引一个 FASTA 文件（明文）。
func OpenFASTA(path string) (*FASTA, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	fr := &FASTA{f: f, idx: make(map[string]fastaIndex)}
	if err := fr.buildIndex(); err != nil {
		_ = f.Close()
		return nil, err
	}
	return fr, nil
}

func (f *FASTA) buildIndex() error {
	if _, err := f.f.Seek(0, io.SeekStart); err != nil {
		return err
	}
	br := bufio.NewReader(f.f)
	var cur string
	var entry fastaIndex
	var lineNo int64
	var pos int64 // 逻辑流位置（字节）
	flush := func() {
		if cur != "" {
			f.idx[cur] = entry
		}
	}
	for {
		line, err := br.ReadString('\n')
		if err != nil && !errors.Is(err, io.EOF) {
			return err
		}
		if len(line) > 0 {
			pos += int64(len(line))
			trimmed := strings.TrimRight(line, "\r\n")
			switch {
			case strings.HasPrefix(trimmed, ">"):
				flush()
				name := strings.TrimPrefix(trimmed, ">")
				if fi := strings.Fields(name); len(fi) > 0 {
					name = fi[0]
				}
				cur = name
				entry = fastaIndex{offset: pos}
				lineNo = 0
			case cur != "" && len(trimmed) > 0:
				if lineNo == 0 {
					entry.lineBases = int64(len(trimmed))
					entry.lineBytes = int64(len(line))
				}
				entry.length += int64(len(trimmed))
				lineNo++
			}
		}
		if errors.Is(err, io.EOF) {
			break
		}
	}
	flush()
	return nil
}

func (f *FASTA) index(chrom string) (fastaIndex, error) {
	idx, ok := f.idx[chrom]
	if !ok {
		return fastaIndex{}, fmt.Errorf("%w: %s", ErrUnknownChrom, chrom)
	}
	return idx, nil
}

func (f *FASTA) Prefix(chrom string, end, length int) (string, error) {
	idx, err := f.index(chrom)
	if err != nil {
		return "", err
	}
	if length < 0 || int64(end) < int64(length) || int64(end) > idx.length {
		return "", fmt.Errorf("%w: %s:%d length=%d", ErrOutOfBounds, chrom, end, length)
	}
	if length == 0 {
		return "", nil
	}
	first := int64(end - length + 1) // 1-based 起
	last := int64(end)
	startZero := first - 1
	off := idx.offset + (startZero/idx.lineBases)*idx.lineBytes + (startZero % idx.lineBases)
	n := last - first + 1
	// 跨行读取，多读的换行在返回时剔除。
	span := n + (last-1)/idx.lineBases - (first-1)/idx.lineBases
	buf := make([]byte, span*2) // 预留 CRLF 余量
	if _, err := f.f.ReadAt(buf, off); err != nil && !errors.Is(err, io.EOF) {
		return "", err
	}
	out := make([]byte, 0, n)
	for _, b := range buf {
		if b != '\r' && b != '\n' {
			out = append(out, b)
		}
	}
	return strings.ToUpper(string(out[:min(int(n), len(out))])), nil
}

func (f *FASTA) Length(chrom string) int {
	idx, err := f.index(chrom)
	if err != nil {
		return 0
	}
	return int(idx.length)
}

// Close 关闭底层文件。
func (f *FASTA) Close() error { return f.f.Close() }
