package genome

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

// Genome is an in-memory reference keyed by sequence name. It is intended for
// chromosome-scale references used during normalization; the streaming import
// path itself never buffers input, only this reference.
type Genome struct {
	seq map[string]string
}

// LoadFASTA reads a plain or bgzipped-text FASTA file into memory.
func LoadFASTA(path string) (*Genome, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseFASTA(f)
}

// ParseFASTA reads FASTA text and concatenates multi-line sequence blocks.
func ParseFASTA(r io.Reader) (*Genome, error) {
	g := &Genome{seq: map[string]string{}}
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)

	var name string
	var builder strings.Builder
	flush := func() {
		if name != "" {
			g.seq[name] = strings.ToUpper(builder.String())
		}
		builder.Reset()
	}

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			flush()
			name = strings.Fields(line[1:])[0]
			continue
		}
		if name == "" {
			return nil, fmt.Errorf("sequence data before any FASTA header")
		}
		builder.WriteString(line)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	flush()
	return g, nil
}

// Substring returns upper-cased bases in the 1-based half-open range
// [start, end).
func (g *Genome) Substring(chrom string, start, end int) (string, error) {
	seq, ok := g.seq[chrom]
	if !ok {
		return "", fmt.Errorf("unknown chromosome %q", chrom)
	}
	if start < 1 || end > len(seq)+1 || start > end {
		return "", fmt.Errorf("range %d:%d out of bounds for %s (len %d)",
			start, end, chrom, len(seq))
	}
	return seq[start-1 : end-1], nil
}

// Length returns the number of bases for chrom.
func (g *Genome) Length(chrom string) (int, error) {
	seq, ok := g.seq[chrom]
	if !ok {
		return 0, fmt.Errorf("unknown chromosome %q", chrom)
	}
	return len(seq), nil
}
