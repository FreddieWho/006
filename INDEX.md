# 006 项目索引

当前活跃版本：**v7**
更新时间：2026-08-31

## 快速导航

| 目的 | 文件 |
|---|---|
| 了解科学问题与创新 | [docs/plan/v7/plan.md](docs/plan/v7/plan.md) |
| 了解阶段、依赖与判定 | [docs/plan/v7/roadmap.md](docs/plan/v7/roadmap.md) |
| 查看当前真实进度 | [results/v7/STATUS.md](results/v7/STATUS.md) |
| 查看 v6 哪些复用/废弃 | [docs/v7/V6_ASSET_DISPOSITION.md](docs/v7/V6_ASSET_DISPOSITION.md) |
| 查看 013 空间资产如何接入 | [docs/v7/SPATIAL_ASSET_INTAKE.md](docs/v7/SPATIAL_ASSET_INTAKE.md) |
| 查看迁移和 Git 边界 | [docs/v7/REPOSITORY_MIGRATION.md](docs/v7/REPOSITORY_MIGRATION.md) |
| 查看科学一致性复核 | [docs/v7/V7_PLAN_SCIENTIFIC_REVIEW.md](docs/v7/V7_PLAN_SCIENTIFIC_REVIEW.md) |
| 查看候选数据角色 | [config/v7/data_source_registry.tsv](config/v7/data_source_registry.tsv) |
| 查看继承清单 | [results/v7/inheritance/asset_manifest.tsv](results/v7/inheritance/asset_manifest.tsv) |
| 查看目录规范 | [STRUCTURE.md](STRUCTURE.md) |
| 查看 agent 规则 | [AGENTS.md](AGENTS.md) |

## 活跃目录

| 路径 | 角色 |
|---|---|
| `docs/plan/v7/` | 科学方案和路线图 |
| `docs/v7/` | 当前状态、资产和迁移说明 |
| `config/v7/` | 小型数据源与运行配置 |
| `scripts/v7/` | 可重放命令入口 |
| `src/v7/` | 可复用实现 |
| `tests/v7/` | v7 接口、泄漏和统计不变量测试 |
| `results/v7/` | 当前状态、manifest 和未来结果 |

## 只读复用与历史

| 路径 | 状态 |
|---|---|
| `scripts/v6_1/`, `scripts/v6_2/` | 历史可重放代码；组件可审计复用 |
| `results/v6_1/`, `results/v6_2/` | 本地历史结果；不作为 v7 结论，不进普通 Git |
| `docs/archive/plans/` | v5/v6 科学方案 |
| `docs/archive/operations/` | v6 合同、审计和 realignment 文档 |
| `scripts/archive/` | 临时诊断与废弃入口 |
| `config/pwy_filter_*`, `src/pathway_filter/` | 早期 pathway pipeline，原位保留以维护历史 import；不是 v7 入口 |
| `mvp/` | 本地 legacy 原型，逻辑归档，不进 Git |

## 外部空间来源

`/home/huyudi/013_spatial` 保持独立仓库和只读来源。v7 通过 registry/manifest 接入，不复制 1.9 TB 数据，不修改其当前运行任务。

## 版本状态

| 版本 | 状态 |
|---|---|
| v5 | archived |
| v6.1 | archived; selected assets reusable |
| v6.2 | superseded; measurement/interface assets reusable |
| v7 | active; scientific computation NOT_RUN |
