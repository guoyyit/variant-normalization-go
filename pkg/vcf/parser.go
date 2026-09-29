package vcf

import (
	"bufio"
	"errors"
	"io"
	"strconv"
	"strings"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// ParseError is a line-annotated parse failure.
type ParseError struct {
	Line   int64
	Reason string
}

func (e *ParseError) Error() string {
	return "vcf line " + strconv.FormatInt(e.Line, 10) + ": " + e.Reason
}

// Parser streams VCF text records. Header lines and comments are skipped;
// only the fixed VCF columns plus GT/DP per sample are materialized.
type Parser struct {
	reader  *bufio.Reader
	counter *countingReader
	line    int64
	offset  int64
	samples []string
}

// countingReader tracks the exact bytes handed to the scanner so checkpoints
// can seek back to record boundaries.
type countingReader struct {
	r      io.Reader
	offset int64
}

func (c *countingReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	c.offset += int64(n)
	return n, err
}

// NewParser creates a streaming parser. The scanner uses a large fixed buffer
// so very long records do not force growth proportional to input size.
func NewParser(r io.Reader) *Parser {
	counter := &countingReader{r: r}
	return &Parser{
		reader:  bufio.NewReaderSize(counter, 64*1024),
		counter: counter,
	}
}

// Samples returns sample names declared by the #CHROM header.
func (p *Parser) Samples() []string { return p.samples }

// Line returns the current 1-based input line number.
func (p *Parser) Line() int64 { return p.line }

// SetLine sets the current line counter, used when resuming a seeked stream.
func (p *Parser) SetLine(n int64) { p.line = n }

// ByteOffset returns the exact byte position in the underlying input up to
// the end of the most recently delivered line.
func (p *Parser) ByteOffset() int64 {
	return p.offset
}

// Next returns the next record or io.EOF when the stream is exhausted.
func (p *Parser) Next() (variant.Record, error) {
	for {
		raw, err := p.reader.ReadBytes('\n')
		if len(raw) > 0 {
			// Logical position excludes bytes still prefetched in buffer.
			p.offset = p.counter.offset - int64(p.reader.Buffered())
		}
		if err != nil && !errors.Is(err, io.EOF) {
			return variant.Record{}, err
		}
		if len(raw) == 0 && errors.Is(err, io.EOF) {
			return variant.Record{}, io.EOF
		}
		p.line++
		line := strings.TrimRight(string(raw), "\r\n")
		if line == "" {
			if errors.Is(err, io.EOF) {
				return variant.Record{}, io.EOF
			}
			continue
		}
		if strings.HasPrefix(line, "##") {
			if errors.Is(err, io.EOF) {
				return variant.Record{}, io.EOF
			}
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			cols := strings.Split(line, "\t")
			if len(cols) > 9 {
				p.samples = cols[9:]
			}
			if errors.Is(err, io.EOF) {
				return variant.Record{}, io.EOF
			}
			continue
		}
		if strings.HasPrefix(line, "#") {
			if errors.Is(err, io.EOF) {
				return variant.Record{}, io.EOF
			}
			continue
		}
		rec, err := parseDataLine(line, p.line, p.samples)
		if err != nil {
			return variant.Record{}, err
		}
		return rec, nil
	}
}

func parseDataLine(line string, lineNo int64, sampleNames []string) (variant.Record, error) {
	cols := strings.Split(line, "\t")
	if len(cols) < 8 {
		return variant.Record{}, &ParseError{Line: lineNo, Reason: "expected at least 8 tab-separated columns"}
	}
	pos, err := strconv.Atoi(cols[1])
	if err != nil {
		return variant.Record{}, &ParseError{Line: lineNo, Reason: "invalid POS: " + cols[1]}
	}

	rec := variant.Record{
		Chrom:  cols[0],
		Pos:    pos,
		Ref:    cols[3],
		Alt:    splitAlt(cols[4]),
		Info:   parseInfo(cols[7]),
		Filter: parseList(cols[6]),
	}
	if cols[2] != "." {
		rec.ID = parseList(cols[2])
	}
	if cols[5] != "." {
		qual, err := strconv.ParseFloat(cols[5], 64)
		if err != nil {
			return variant.Record{}, &ParseError{Line: lineNo, Reason: "invalid QUAL: " + cols[5]}
		}
		rec.HasQual = true
		rec.Qual = qual
	}

	if len(cols) >= 10 {
		format := strings.Split(cols[8], ":")
		for i, values := range cols[9:] {
			name := strconv.Itoa(i + 1)
			if i < len(sampleNames) {
				name = sampleNames[i]
			}
			call, err := parseSample(name, format, strings.Split(values, ":"))
			if err != nil {
				return variant.Record{}, &ParseError{Line: lineNo, Reason: err.Error()}
			}
			rec.Samples = append(rec.Samples, call)
		}
	}
	return rec, nil
}

func splitAlt(field string) []string {
	if field == "." || field == "" {
		return nil
	}
	return strings.Split(field, ",")
}

func parseList(field string) []string {
	if field == "." || field == "" {
		return nil
	}
	return strings.Split(field, ";")
}

func parseInfo(field string) map[string]string {
	if field == "." || field == "" {
		return nil
	}
	out := map[string]string{}
	for _, item := range strings.Split(field, ";") {
		if item == "" {
			continue
		}
		if key, value, ok := strings.Cut(item, "="); ok {
			out[key] = value
		} else {
			out[item] = ""
		}
	}
	return out
}

func parseSample(name string, format, values []string) (variant.SampleCall, error) {
	call := variant.SampleCall{Sample: name, GT: []int{}, Fields: map[string]string{}}
	for i, key := range format {
		if i >= len(values) {
			break
		}
		value := values[i]
		if key == "GT" {
			gt, phased, err := parseGT(value)
			if err != nil {
				return variant.SampleCall{}, err
			}
			call.GT = gt
			call.Phased = phased
			continue
		}
		if value != "." {
			call.Fields[key] = value
		}
	}
	return call, nil
}

func parseGT(field string) ([]int, bool, error) {
	phased := strings.Contains(field, "|")
	parts := strings.FieldsFunc(field, func(r rune) bool { return r == '/' || r == '|' })
	gt := make([]int, 0, len(parts))
	for _, part := range parts {
		if part == "." {
			gt = append(gt, -1)
			continue
		}
		n, err := strconv.Atoi(part)
		if err != nil {
			return nil, false, errors.New("invalid GT value: " + part)
		}
		gt = append(gt, n)
	}
	return gt, phased, nil
}
