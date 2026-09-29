// Package config 定义跨来源合并的保留策略与过滤阈值。
//
// 配置在每次导入开始时读取并随该批次一起持久化；配置变更只影响后续导入，
// 已入库记录不会被回溯改写，因此历史合并结果始终与其当时配置自洽。
package config

import (
	"encoding/json"
	"fmt"
	"os"
)

// Config 为合并保留策略与过滤阈值配置。
type Config struct {
	// SourcePriority 来源 ID 到优先级的映射；数值越小优先级越高。
	// 未列出的来源使用 DefaultPriority。
	SourcePriority map[string]int `json:"source_priority"`

	// DefaultPriority 是未显式配置来源的默认优先级。
	DefaultPriority int `json:"default_priority"`

	// MinQual 质量阈值：合并后记录质量低于该值时 FILTER 标记为 LowQual。
	MinQual float64 `json:"min_qual"`

	// ConfigVersion 配置版本，随每次导入的记录一起持久化，便于审计。
	ConfigVersion string `json:"config_version"`
}

// Default 返回内置默认配置：所有来源等优先级、不过滤。
func Default() *Config {
	return &Config{
		SourcePriority:  map[string]int{},
		DefaultPriority: 100,
		MinQual:         0,
		ConfigVersion:   "default",
	}
}

// Load 从 JSON 文件加载配置；文件不存在时返回默认配置。
func Load(path string) (*Config, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return Default(), nil
		}
		return nil, err
	}
	cfg := Default()
	if err := json.Unmarshal(data, cfg); err != nil {
		return nil, fmt.Errorf("parse config %s: %w", path, err)
	}
	if cfg.SourcePriority == nil {
		cfg.SourcePriority = map[string]int{}
	}
	return cfg, nil
}

// Priority 返回某来源的优先级（数值越小越优先）。
func (c *Config) Priority(sourceID string) int {
	if p, ok := c.SourcePriority[sourceID]; ok {
		return p
	}
	return c.DefaultPriority
}
