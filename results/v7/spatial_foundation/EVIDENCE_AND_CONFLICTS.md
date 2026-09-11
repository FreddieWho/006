# Stage 3 证据与冲突

## 扩展资产核对（Stage 4 资产连同一合同）

- Stage 4 选择 60 logical units／71 captures 复用 Stage 3 同一 adapter/身份/分割/坐标合同：`results/v7/spatial_discovery/stage4_selection.tsv`（含 stage1_registry_sha256、r04_manifest_sha256）＋`model_safe_input_manifest.tsv`（model-safe 清单通过禁止字段检查）。
- 新增/改变路径才做技术回放：Stage 4 新增捕获超出 Stage 3 的 11 捕获部分，已在 Stage 4 `unit_discovery_audit.tsv`（71 行 PASS）与 `model_safe_input_manifest.tsv` 中做技术审计；未改变的合同不重复全量回放。
- kNN 邻接≠真实细胞接触：图为 section-local kNN k=6 binary（`config/v7/stage4.yaml:graph`），无图像/分割时不伪造物理邻接；下游不得作接触/屏障解释。

## 证据行

| claim | evidence | locator | 独立单位 | 暴露 |
|---|---|---|---|---|
| 真实输入链路闭合 | 11/11 通过 counts/坐标/identity/adapter QC | `adapter_qc.tsv`、`count_conservation.tsv` | 捕获（技术），患者只汇总 | response-blind |
| 跨平台投影可复核 | 同一 Stage 2 schema 的 score/coverage；1 共同特征非退化 | `feature_coverage.tsv`、`graph_statistics.parquet` | 同上 | 同上 |
| 统计语义稳定 | section 内 exchangeability；同患者不跨 fold | `spatial_null_statistics.parquet`、`split_audit.tsv` | 患者外层 | 同上 |
| 隔离边界通过 | model-safe 禁止字段；GT 封存 locator/hash 未解析值 | `gt_isolation_audit.tsv`、`sealed_validation_manifest.tsv` | NOT_APPLICABLE | GT 值未读 |

对照：99 次 section 内置换 null；1,000 次患者 bootstrap。
冲突：无科学冲突；Xenium 上游 h5ad uns/response 分支已在盲副本移除（原文未改）。
缺失/未运行：image/segmentation/composition 缺失显式为未提供；共同可测特征受 panel 限制；本阶段不回答 TLS/屏障/疗效/repair。
