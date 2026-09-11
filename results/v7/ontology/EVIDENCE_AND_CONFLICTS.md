# Stage 2 证据与冲突

## 主分析特征连接（覆盖率/carrier/null/稳定性/score 语义）

定位：`measurement_reliability.tsv`（主键 feature_id，49 行含父轴；可执行 39 特征）。

| 列族 | 字段 | 语义 |
|---|---|---|
| 覆盖率 | median_coverage、n_expression_units、n_study_families | 缺失保留，不补零 |
| carrier | carrier_dependency_status、abundance_equal_state_spearman | `CONCORDANT_WITH_RESIDUAL_COMPOSITION_RISK` 保留为风险，不删特征 |
| null | matched_null_ci_high、observed_coherence_ci_low、n_matched_contexts | exact-mask null；不足记 `NOT_TESTABLE` |
| 跨研究稳定性 | leave_family_out_median_spearman±CI | 单队列驱动列为削弱条件 |
| score 语义 | technical_status、D2_class、D2_reason、method_agreement_*、technical_panel_* | 7 measurable／32 uncertainty；区间为算法稳定/置换包络，非生物效应 CI |

三类输入：主要（7 measurable）、敏感性（39 全量）、仅注释（父轴向量/不可校准桥）。
独立单位：表达单元/患者—时间点（非临床独立单位）；临床独立单位为患者（Stage 5 起）。
暴露记录：`selection_exposure.tsv`（response-blind）； bulk/perturbation 不继承 scRNA D2。

## Stage 7 pre/post 可比尺度核对

- Stage 7 输入为同一 Stage 2 分数 schema（`patient_timepoint_scores.parquet`），pre/post 同一尺度。
- Stage 7 仅去技术 nuisance，保留 arm、time、arm×time、patient pairing；composition-inclusive 与 composition-adjusted 并列（见 Stage 7 DECISION_LOG）。核对通过。

## 对照/冲突/缺失

- 对照：随机基因集、exact-mask matched null、leave-study-family-out（见 diagnostics）。
- 冲突：逐细胞先归一化 vs counts 先汇总再归一化为不同估计量，只检排序稳健性（中位 rank concordance 0.758），不要求数值相等。
- 缺失：GSE207422 精确配对仅 P05/P08（Spearman 0.729/0.282，描述性）；GSE193736 只给重复性接口（median Spearman 0.881），不给方向；真实空间投影 NOT_RUN_PENDING_STAGE3。
- 结论影响：bulk 两个 exact pairs 不作特征验证；扰动重复性不作方向证据；D2 uncertainty 不得写成已确认机制。
