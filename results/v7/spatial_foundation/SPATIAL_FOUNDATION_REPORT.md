# v7 Stage 3 空间基础回放报告

## 状态

**D3_PASS_WITH_LIMITATIONS**。本次回放覆盖 8 个逻辑 pilot 单元、11 个真实技术捕获、7 个患者身份包，平台为 Visium、Xenium，数据集为 GSE238264、GSE291246、HTAN_VANDERBILT_CRC。置换 null 使用每个 section 内 99 次重排；患者 bootstrap 使用 1,000 次重采样。

## 技术路线

- 先用 Stage1 技术白名单和 R-04 只读定位清单复核 pilot；选择不读取 response、治疗、时间点或 GT 字段。
- HTAN 使用空间 h5ad；GSE238264 使用 10x Visium H5、positions 和 scalefactors；GSE291246 使用已 QC 的 Xenium cell-feature h5ad。
- 所有输入统一为非负整数 raw counts、原生 observation ID、二维坐标和患者—区块—捕获层级；每个 section 单独建稀疏 kNN 图，禁止跨 section 连边。
- 将 Stage2 的 39 个可评分特征投影到空间观测。每个特征同时保存 native score、section 内标准化 score、预期/实际基因和覆盖状态；缺失基因保持缺失，不补零。
- 对连续 score 计算 Moran 型自相关、稀疏 variogram 和 section 内置换 null；再按 section→patient→cohort 等权汇总。没有可信 composition 时输出明确的不可估计状态。

## 结果 / 证据

1. **真实输入链路闭合。** 11/11 个捕获通过 counts、坐标、observation identity 和 adapter QC；原始计数总量与稀疏矩阵审计通过（见 `adapter_qc.tsv`、`count_conservation.tsv`）。
2. **跨平台投影可复核。** Visium 与 Xenium 均产生了同一 Stage2 schema 的 score/coverage 表；1 个特征在两类平台均为完整面板且至少 1 个共同特征生成了非退化空间统计（见 `feature_coverage.tsv`、`graph_statistics.parquet`）。
3. **空间统计语义稳定。** 图构建为 section-local 稀疏图，null 置换保持 section 内 exchangeability；患者 split 中同一患者未跨 fold（见 `spatial_null_statistics.parquet`、`split_audit.tsv`）。
4. **隔离边界通过。** model-safe 清单通过禁止字段检查；HTAN 结构锚点仅以封存 locator/hash 验证，GT 值未解析；GSE175540、USZ TLS 和 ST_CRC_CMS 仍锁定（见 `gt_isolation_audit.tsv`、`sealed_validation_manifest.tsv`）。
5. **平台差异被保留。** Visium 保留 spot/array/pixel 与 scalefactor 状态；Xenium 保留 cell ID、micrometer 坐标和 targeted panel 语义。缺失 image、segmentation、composition 或物理尺度不被伪造为可用。

## 结论

在 response-blind 条件下，本项目已经建立了可把 Visium、Xenium 和 h5ad 真实计数映射到同一空间统计合同的底座，并能从捕获层安全汇总到患者层。这个证据支持 Stage4 开始做结构发现和模型比较，但不支持把本阶段的空间关联称为 TLS、屏障、PD-1 失败机制或 PD1+X 修复证据。

## 已知限制

- Xenium 上游 h5ad 的 uns/response 分支已在响应盲副本中移除；原始文件未修改。
- Xenium 未提供 image/micron 之外的共同图像或真实 composition；对应分析显式为未提供/不适用。
- 共同可测特征数量由平台面板决定，缺失成员保留为 limited_with_uncertainty，不补零。
- 本阶段只完成空间输入与统计底座，不回答 TLS、空间屏障、PD-1 疗效或 PD1+X 修复。

## 可复现入口

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \
  python scripts/v7/spatial/run_stage3.py --permutations 99 --bootstrap-draws 1000
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \
  python scripts/v7/spatial/evaluate_gt_isolation.py
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \
  python scripts/v7/spatial/finalize_stage3.py
```
