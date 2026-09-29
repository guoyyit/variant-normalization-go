package config

import (
	"os"
	"testing"
)

func TestLoadPriority(t *testing.T) {
	f, _ := os.CreateTemp("", "cfg*.json")
	defer os.Remove(f.Name())
	f.WriteString(`{"id":"v2","source_priority":{"clinvar":10},"default_priority":1,"min_qual":20,"qual_agg":"weighted","filter_mode":"all"}`)
	f.Close()
	cfg, err := LoadFile(f.Name())
	if err != nil {
		t.Fatal(err)
	}
	if cfg.Priority("clinvar") != 10 || cfg.Priority("other") != 1 {
		t.Fatal("priority lookup wrong")
	}
	bad := Default()
	bad.QualAgg = "nope"
	if bad.Validate() == nil {
		t.Fatal("expected validation error")
	}
}
