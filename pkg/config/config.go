package config

import (
	"encoding/json"
	"fmt"
	"os"
)

// Config controls merge retention and recomputed filter status.
type Config struct {
	// ID identifies the configuration version carried by merged records.
	ID string `json:"id"`
	// SourcePriority maps source names to priority values; higher wins.
	SourcePriority map[string]int `json:"source_priority"`
	// DefaultPriority is used for unknown sources.
	DefaultPriority int `json:"default_priority"`
	// MinQual is the threshold below which the merged record is FILTERed.
	MinQual float64 `json:"min_qual"`
	// MinAF is the allele-frequency threshold applied on recompute.
	MinAF float64 `json:"min_af"`
	// QualAgg is one of: "max" (default), "weighted", "min".
	QualAgg string `json:"qual_agg"`
	// FilterMode is "all" (default; PASS only if every source PASSes) or
	// "any" (PASS if at least one source PASSes).
	FilterMode string `json:"filter_mode"`
}

// Default returns the built-in configuration.
func Default() Config {
	return Config{
		ID:              "default",
		SourcePriority:  map[string]int{},
		DefaultPriority: 0,
		MinQual:         0,
		MinAF:           0,
		QualAgg:         "max",
		FilterMode:      "all",
	}
}

// Priority returns the configured priority of source, or the default.
func (c Config) Priority(source string) int {
	if p, ok := c.SourcePriority[source]; ok {
		return p
	}
	return c.DefaultPriority
}

// Validate checks configuration invariants and fills defaults.
func (c *Config) Validate() error {
	if c.QualAgg == "" {
		c.QualAgg = "max"
	}
	switch c.QualAgg {
	case "max", "min", "weighted":
	default:
		return fmt.Errorf("unknown qual_agg %q (want max, min or weighted)", c.QualAgg)
	}
	if c.FilterMode == "" {
		c.FilterMode = "any"
	}
	switch c.FilterMode {
	case "any", "all":
	default:
		return fmt.Errorf("unknown filter_mode %q (want any or all)", c.FilterMode)
	}
	if c.MinQual < 0 {
		return fmt.Errorf("min_qual must be >= 0")
	}
	if c.MinAF < 0 || c.MinAF > 1 {
		return fmt.Errorf("min_af must be within [0,1]")
	}
	if c.SourcePriority == nil {
		c.SourcePriority = map[string]int{}
	}
	if c.ID == "" {
		c.ID = "default"
	}
	return nil
}

// LoadFile reads a JSON configuration file.
func LoadFile(path string) (Config, error) {
	cfg := Default()
	data, err := os.ReadFile(path)
	if err != nil {
		return Config{}, err
	}
	if err := json.Unmarshal(data, &cfg); err != nil {
		return Config{}, fmt.Errorf("parse config %s: %w", path, err)
	}
	if err := cfg.Validate(); err != nil {
		return Config{}, err
	}
	return cfg, nil
}
