# v7 configuration

本目录只存可审计的小型配置和数据源角色，不存原始数据。

- `data_source_registry.tsv`：经 Stage 1 校准过的候选数据角色和有效规模摘要；逐 logical-unit 事实以 `results/v7/registry/` 为准。
- 未来运行配置必须记录输入版本、患者级 split、随机种子、输出目录和失败恢复方式。
