# v7 仓库迁移记录

日期：2026-08-31

## 目标结构

- `docs/plan/v7/`：唯一活跃 plan/roadmap；
- `docs/archive/plans/`：v5、v6.1、v6.2 历史方案；
- `docs/archive/operations/`：v6.1/v6.2 历史执行和审计文档；
- `scripts/v6_1/`、`scripts/v6_2/`：保留可重放的历史代码路径；
- `scripts/archive/`：临时诊断或不再作为入口的脚本；
- `config/v7/`、`scripts/v7/`、`src/v7/`、`tests/v7/`：v7 新工作；
- `results/v7/`：v7 轻量状态、manifest 和未来结果；
- `data/`、`results/v6_1/`、`results/v6_2/`：本地大型资产，不进入普通 Git。

## 本次物理移动

| 原位置 | 新位置 | 性质 |
|---|---|---|
| `docs/plan/v5*.md` | `docs/archive/plans/v5/` | 历史方案 |
| `docs/plan/v6_1/` | `docs/archive/plans/v6_1/` | 历史方案 |
| `docs/plan/v6_2/` | `docs/archive/plans/v6_2/` | v7 取代后的历史方案 |
| `docs/v6_1/` | `docs/archive/operations/v6_1/` | 历史运营文档 |
| `docs/v6_2/` | `docs/archive/operations/v6_2/` | 历史 realignment/探索文档 |
| 根目录临时 `check_*`、`debug_*`、`test_*` | `scripts/archive/v6_2_diagnostics/` | 临时诊断代码 |
| 历史下载脚本 | `scripts/archive/v6_data_acquisition/` | v6 数据获取溯源，不作为 v7 入口 |
| 旧目录生成脚本 | `scripts/archive/v5_scaffold/` | 旧结构脚手架，不作为 v7 入口 |
| GSE301741 构建脚本 | `scripts/v6_2/legacy_ingest/` | 保留可重放的数据接入代码 |
| 根目录 `RData/xlsx` | `data/archive/unclassified/` | 本地未分类数据；仍被 Git 忽略 |
| 根目录 anchor 报告 | `docs/archive/reports/` | 历史报告 |
| `reports/pathway_filtering/` | `docs/archive/reports/pathway_filtering/` | 早期 pathway 探索输出 |
| v6.1 横向报告与补充说明 | `docs/archive/reports/v6_1/`、`docs/archive/operations/v6_1/` | 历史证据与数据问题记录 |
| 根目录旧环境清单 | `docs/archive/environments/v6_1/` | 历史环境参考，不作为 v7 lockfile |

## 逻辑归档，不物理移动

- `results/v6_1/` 约 31 GB；
- `results/v6_2/` 约 32 GB；
- `data/` 约 1.1 TB；
- `mvp/` 约 1.5 GB；
- `~/013_spatial` 约 1.9 TB。

原因：大量历史脚本和 manifest 使用绝对/版本化路径，移动会破坏可重放性或干扰独立项目。v7 通过资产清单引用，保留只读来源。

## GitHub 提交边界

提交：代码、测试、Markdown、配置、小型 registry/manifest。
不提交：raw data、h5ad/rds/qs/RData、大型 parquet/CSV、图像、checkpoint、scratch、cache、日志、本地 IDE/agent 索引。

原本单一未推送根提交包含超过 GitHub 100 MB 限制的文件。本次已建立仅本地的 `pre_v7_local_backup`（指向 `a1a9e80`），并从 `main` 索引移除 `results/v6_1/`；3,226 个本地 v6.1 文件仍留在磁盘。`main` 只发布轻量代码、文档和清单，不推送本地大结果，也不推送该备份分支。

## 保留原位的 legacy 代码

`config/pwy_filter_*`、`src/pathway_filter/`、`scripts/preprocessing/`、`scripts/v6_1/`、`scripts/v6_2/` 和 `tools/v6_1/` 暂保原位，以免破坏历史 import、数据入口和结果可重放性。它们都不是 v7 命令入口；未来若复用，先复制最小组件到 `src/v7/` 并补来源记录。
