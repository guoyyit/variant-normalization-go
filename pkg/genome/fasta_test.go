package genome

import (
	"strings"
	"testing"
)

func TestParseFASTA(t *testing.T) {
	in := ">chr1 description\nACGT\nacgt\n>chr2\nTTTT\n"
	g, err := ParseFASTA(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	if n, _ := g.Length("chr1"); n != 8 {
		t.Fatalf("len=%d", n)
	}
	got, err := g.Substring("chr1", 5, 9)
	if err != nil {
		t.Fatal(err)
	}
	if got != "ACGT" {
		t.Fatalf("bases=%q (must be uppercased)", got)
	}
	if _, err := g.Substring("nope", 1, 2); err == nil {
		t.Fatal("expected unknown chromosome error")
	}
	if _, err := g.Substring("chr1", 9, 10); err == nil {
		t.Fatal("expected out-of-bounds error")
	}
}
