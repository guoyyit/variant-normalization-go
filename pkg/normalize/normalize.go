package normalize

import (
	"errors"
	"fmt"
	"strings"

	"github.com/guoyyit/variant-normalization-go/pkg/variant"
)

// Reference provides reference sequence for left alignment.
// Positions are 1-based; Substring returns bases in [start, end).
type Reference interface {
	Substring(chrom string, start, end int) (string, error)
	// Length returns the number of bases available for chrom.
	Length(chrom string) (int, error)
}

var (
	// ErrInvalidFormat means the record text/fields are not a valid VCF record.
	ErrInvalidFormat = errors.New("invalid variant format")
	// ErrOutOfBounds means coordinates reference bases outside the chromosome.
	ErrOutOfBounds = errors.New("coordinate out of bounds")
)

// Error is a distinguishable normalization failure.
type Error struct {
	Code   error
	Reason string
}

func (e *Error) Error() string { return e.Code.Error() + ": " + e.Reason }
func (e *Error) Unwrap() error { return e.Code }

// Normalize validates rec and returns its canonical representation:
// alleles are upper-cased, common suffix/prefix anchors are stripped and
// indels are left aligned against the reference. Normalize expects records
// already split to a single ALT (see package split).
func Normalize(rec variant.Record, ref Reference) (variant.Record, error) {
	out := rec.Clone()

	if strings.TrimSpace(out.Chrom) == "" {
		return variant.Record{}, invalid("empty CHROM")
	}
	if out.Pos < 1 {
		return variant.Record{}, invalidF("POS must be >= 1, got %d", out.Pos)
	}
	if len(out.Alt) != 1 {
		return variant.Record{}, invalidF("normalization requires exactly one ALT, got %d", len(out.Alt))
	}

	refAllele := strings.ToUpper(out.Ref)
	if !isBases(refAllele) {
		return variant.Record{}, invalidF("invalid REF allele %q", out.Ref)
	}
	out.Ref = refAllele

	alt := out.Alt[0]
	if strings.TrimSpace(alt) == "" {
		return variant.Record{}, invalid("empty ALT allele")
	}
	if variant.IsSymbolic(alt) || variant.IsBreakend(alt) {
		if err := validateBounds(out, ref, len(refAllele)); err != nil {
			return variant.Record{}, err
		}
		if err := validateGTs(out); err != nil {
			return variant.Record{}, err
		}
		return out, nil
	}
	if variant.IsNoAlt(alt) {
		if err := validateBounds(out, ref, len(refAllele)); err != nil {
			return variant.Record{}, err
		}
		if err := validateGTs(out); err != nil {
			return variant.Record{}, err
		}
		return out, nil
	}

	alt = strings.ToUpper(alt)
	if !isBases(alt) {
		return variant.Record{}, invalidF("invalid ALT allele %q", out.Alt[0])
	}
	if alt == refAllele {
		return variant.Record{}, invalidF("ALT must differ from REF (%s)", alt)
	}
	out.Alt[0] = alt

	if err := validateBounds(out, ref, len(refAllele)); err != nil {
		return variant.Record{}, err
	}
	if err := validateGTs(out); err != nil {
		return variant.Record{}, err
	}

	if ref == nil {
		return variant.Record{}, outOfBoundsF("reference genome required to normalize %s:%d", out.Chrom, out.Pos)
	}

	nRef, nPos, nAlt, err := leftAlign(out.Chrom, out.Pos, refAllele, alt, ref)
	if err != nil {
		return variant.Record{}, err
	}
	out.Ref = nRef
	out.Alt[0] = nAlt
	out.Pos = nPos

	// Normalized records never carry stale pre-normalization annotations.
	delete(out.Info, "AC")
	delete(out.Info, "AN")
	delete(out.Info, "AF")
	return out, nil
}

func isBases(s string) bool {
	if s == "" {
		return false
	}
	for _, b := range s {
		switch b {
		case 'A', 'C', 'G', 'T', 'N':
		default:
			return false
		}
	}
	return true
}

func validateBounds(rec variant.Record, ref Reference, span int) error {
	if ref == nil {
		return nil
	}
	length, err := ref.Length(rec.Chrom)
	if err != nil {
		return outOfBoundsF("cannot resolve chromosome %q: %v", rec.Chrom, err)
	}
	end := rec.Pos + span - 1
	if end > length {
		return outOfBoundsF("%s:%d-%d exceeds chromosome length %d", rec.Chrom, rec.Pos, end, length)
	}
	return nil
}

func validateGTs(rec variant.Record) error {
	for _, sm := range rec.Samples {
		for _, allele := range sm.GT {
			if allele == -1 {
				continue
			}
			if allele < 0 || allele > len(rec.Alt) {
				return invalidF("sample %q genotype allele index %d out of range [0,%d]",
					sm.Sample, allele, len(rec.Alt))
			}
		}
	}
	return nil
}

// trimCommon strips shared suffix then prefix bases, always keeping at least
// one anchor base on each allele.
func trimCommon(r, a string) (string, string) {
	for len(r) > 1 && len(a) > 1 && r[len(r)-1] == a[len(a)-1] {
		r = r[:len(r)-1]
		a = a[:len(a)-1]
	}
	for len(r) > 1 && len(a) > 1 && r[0] == a[0] {
		r = r[1:]
		a = a[1:]
	}
	return r, a
}

// leftAlign implements the canonical left-alignment algorithm (Tan et al.
// 2015, as used by vt/bcftools): walk the allele window one base left until
// it can move no further, trimming shared anchors after each step.
func leftAlign(chrom string, pos int, r, a string, ref Reference) (string, int, string, error) {
	r, a = trimCommon(r, a)
	for {
		if len(r) == len(a) {
			break
		}
		if len(r) == 1 && len(a) == 1 {
			break
		}
		if pos <= 1 {
			break
		}
		leftBase, err := ref.Substring(chrom, pos-1, pos)
		if err != nil {
			return "", 0, "", outOfBoundsF("cannot fetch reference base at %s:%d: %v", chrom, pos-1, err)
		}
		leftBase = strings.ToUpper(leftBase)
		if len(leftBase) != 1 || !isBases(leftBase) {
			return "", 0, "", outOfBoundsF("invalid reference base at %s:%d", chrom, pos-1)
		}
		r = leftBase + r[:len(r)-1]
		a = leftBase + a[:len(a)-1]
		pos--
		r, a = trimCommon(r, a)
	}
	return r, pos, a, nil
}

func invalid(reason string) *Error {
	return &Error{Code: ErrInvalidFormat, Reason: reason}
}

func invalidF(format string, args ...any) *Error {
	return invalid(fmt.Sprintf(format, args...))
}

func outOfBoundsF(format string, args ...any) *Error {
	return &Error{Code: ErrOutOfBounds, Reason: fmt.Sprintf(format, args...)}
}
