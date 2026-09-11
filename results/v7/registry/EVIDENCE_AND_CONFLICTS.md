# Stage 1 证据与冲突

## 资格三列（库存存在 / 本地可运行 / 临床纵向可识别）

| 来源 | 库存存在 | 本地可运行 | 临床/纵向可识别 | 证据行 |
|---|---|---|---|---|
| TASK01（27患者/72样本行） | 是 | 是 | PD-1 anchor 环境内可用；mono/combo 不可交换 | `treatment_response_completeness.tsv`、REGISTRY_REPORT §治疗完整性 |
| TASK02（23/66） | 是 | 是 | PD1+lenvatinib paired/context；不作单药交换对照 | 同上 |
| LAMBRECHT_HCC（44患者/112样本：pre65/post47） | 是 | 是 | 纵向 support；旧 25/22 摘要冲突保留 | 同上＋角色表 reason |
| GSE238264（7 HCC post-only，nivo+cabo；R3/NR4） | 是 | 是 | post-only 空间关联/support；无基线，不分离 PD-1/cabo，不建纵向因果 | 同上 |
| GSE291246（35 BCC sections，17 h5ad 可重放） | 是 | 部分（17/35） | 无 response，只跨癌 support | 同上；患者列以 freeze manifest 为准 |
| GSE211956 | 声明存在 | 待定 | `response_semantics_pending`，未升级 | REGISTRY_REPORT |
| Mendeley HCC（有效 6患者/12样本；本地 h5ad=scRNA） | 是 | 部分 | `no_until_provenance_reconciled` | 同上 |
| 006_bulk 多队列 | 是 | 是 | 按队列审计 treatment/endpoint/mapping 后可用 | 角色表 `external_patient_level_direction_and_clinical_scale` |
| 006_perturbation（GSE133344/GSE193736/GSE306429/GSE90063/XAtlas/L1000） | 是 | 是 | 只作 `intervention_to_program_direction`；cell context 非临床 response | 角色表＋ontology `gse193736_perturbation_design.tsv` |
| 013 HTAN/GSE175540/ST_CRC_CMS/USZ 等 | 是 | 按技术白名单 | discovery/support/validation 按 role ledger；locked external 不进调参 | 角色表＋`role_freeze.tsv` provenance |

独立单位：患者为临床外层；block/section/spot/cell 为嵌套，不增加独立 n；同患者/同块/serial section 进同一 leakage group。
暴露记录：全量 `selection_exposure=not_run_v7_stage1_metadata_audit`（同一 claim 尚未运行选择/阈值/模型）；013 旧 `role_freeze.tsv` 仅 provenance 与历史泄漏边界。
对照：本阶段无科学对照；`complete_for_response=yes` 仅字段完整，不自动授予临床主张资格。

## 冲突

- Lambrecht 旧摘要（25 pre/22 paired）vs 有效主表（44/112）：冲突保留，不沿用旧口径。
- GSE291246 QC 表患者列可能为空：以 metadata-freeze manifest 映射为准。
- 旧 `013_other_spatial` 无角色类别：已消除，全部分配受控 source。
- asset_manifest 11 行 `pending_stage1_manifest`：已由本阶段实际 manifest 闭合。

## 缺失/未运行及结论影响

- 同患者空间 pre/on/post 链：缺失→直接纵向空间 PD1+X rewiring=NOT_IDENTIFIABLE；分子/细胞纵向可继续，空间分支阻塞。
- GSE211956/Mendeley：未闭合→不计入 response-bearing/独立验证；不增加独立患者。
- R-04 K/场结论：未继承（`RESTART_STABILITY_DIAGNOSTIC_COMPLETE_NOT_FORMAL_K_SELECTION`）。
- 科学结论：NOT_RUN；本报告只证明资格、有效样本量与边界。
