# analysis_contract (v6.1)

## 0) Contract objective
本契约用于锁定 v6.1 主分析链路的输入/输出边界，防止“跳步分析”“不可追溯结论”和“旧 MVP 演示结果冒充主证据”。

---

## 1) Step-wise I/O contract

| Step ID | Step name | Allowed inputs (read whitelist) | Required outputs (must produce) |
|---|---|---|---|
| S0 | Cohort/Label/Context Lock | `docs/data_collection.csv`; `mvp/sidecars/*/sample_summary.csv`; 历史 registry 仅用于对照（`mvp/outputs/cohort_registry.csv`, `mvp/outputs/patient_metadata_master.csv`） | `results/v6_1/step0/cohort_registry_v6_1.csv`; `results/v6_1/step0/treatment_context_flags.csv`; `results/v6_1/step0/response_label_dictionary.csv`; `results/v6_1/step0/cohort_qc_audit.csv` |
| S1 | Immune-state Measurement | S0 全部输出；`mvp/sidecars/*/obs.parquet`; `mvp/sidecars/*/var.parquet`; 原始 count matrix 数据源 | `results/v6_1/step1/sample_index_master_v6_1.csv`; `results/v6_1/step1/cell_state_fraction_by_sample.parquet`; `results/v6_1/step1/pseudobulk_by_sample.parquet`; `results/v6_1/step1/immune_state_features.parquet`; `results/v6_1/step1/state_qc_report.csv` |
| S2 | Strong Baseline First | S0/S1 全部输出；预定义 signature/pathway/TF 资源 | `results/v6_1/step2/strong_baseline_results.csv`; `results/v6_1/step2/strong_baseline_cv_predictions.parquet`; `results/v6_1/step2/strong_baseline_model_cards.yaml`; `results/v6_1/step2/negative_control_audit.csv` |
| S3 | Shared/Specific Module Discovery + PD-1 Anchor | S0/S1/S2 全部输出；允许读取 perturb prior（仅辅助约束） | `results/v6_1/step3/shared_specific_modules.csv`; `results/v6_1/step3/immune_module_scores.csv`; `results/v6_1/step3/pd1_anchor_modules.csv`; `results/v6_1/step3/module_stability_report.csv` |
| S4 | PD-1 Anchor -> PD-1+X Repair | S0 `treatment_context_flags`; S3 anchor/module 输出；S2 baseline 输出；perturb prior | `results/v6_1/step4/pd1x_repair_hypotheses.csv`; `results/v6_1/step4/M_signals_v6_1.csv`; `results/v6_1/step4/mechanism_support_precheck.csv` |
| S5 | Spatial/Tissue Biological Adjudication | S3/S4 输出；空间/组织数据（ST/IMC/多重 IHC 等）；LR 数据库 | `results/v6_1/step5/spatial_adjudication_results.csv`; `results/v6_1/step5/lrg_spatial_support.csv`; `results/v6_1/step5/spatial_failure_flags.csv`; `results/v6_1/step5/spatial_adjudication_report.md` |
| S6 | Translational Nomination + External Anchor | S2/S3/S4/S5 全部输出；外部 bulk/clinical anchors | `results/v6_1/step6/PD1X_nomination_cards.csv`; `results/v6_1/step6/external_validation_results.csv`; `results/v6_1/step6/wetlab_priority_list.csv`; `results/v6_1/step6/final_evidence_table.csv` |

---

## 2) No-bypass rules (hard)

1. S1-S6 禁止绕过 S0 直接从原始 metadata 组装 cohort/label/context。
2. S3-S6 禁止绕过 S2 强基线结果直接发布高级模型性能或机制优越性结论。
3. S4 禁止绕过 S3 `pd1_anchor_modules.csv` 直接构建 PD-1+X repair logic。
4. S6 禁止绕过 S5 空间裁判结果发布 L2/L3 机制主结论。
5. S6 禁止直接读取旧 MVP 排名结果生成 nomination（必须从 S3-S5 新链路重建）。
6. 所有主结论必须可追溯到 `output_manifest` 的 step 记录与对应输出文件。

---

## 3) Legacy MVP files: proxy/sanity only
以下文件（或同类文件）只能用于 proxy、sanity check、历史对照，不得作为 v6.1 主证据来源：

- `mvp/outputs/patient_tme_moa_scores.csv`
- `mvp/outputs/patient_tme_moa_scores_refined.csv`
- `mvp/outputs/tme_moa_evidence_table.csv`
- `mvp/outputs/tme_moa_shift_vectors.csv`
- `mvp/outputs/virtual_drug_shift_scores.csv`
- `mvp/outputs/result_ranking.csv`
- `mvp/outputs/bulk_archetype_assignment.csv`
- `mvp/outputs/final_archetype_summary.csv`
- `mvp/outputs/module_archetype_summary.csv`
- `mvp/outputs/showcase_patient_cards.json`
- `mvp/outputs/final_demo_manifest.txt`

说明：
- `mvp/outputs/pseudobulk_by_sample.parquet`、`mvp/outputs/cell_state_fraction_by_sample.parquet` 可用于对照 QC，但 v6.1 主分析必须在 S1 重新导出并写入 `results/v6_1/step1/`。
- `mvp/sidecars/*` 作为数据组织底座可直接继承，但要在 S0/S1 重建 cohort 与样本索引，不可直接复制旧结论表。

---

## 4) Forbidden-in-main-conclusion results
以下结果不得出现在主结论段落（正文主 claim）：

1. 仅基于 marker 均值、简单 marker 组合、或单图可视化直觉的机制判断。
2. 仅基于 LINCS reversal 的药物或靶点提名。
3. 无强基线对照即宣称模型优越。
4. 无空间/组织支持即宣称 TME 机制成立（可降级为 L1 association）。
5. 仅由 Context C（高混杂治疗）驱动的主机制结论。
6. 直接复用 MVP 手工 MOA、heuristic counterfactual、pseudo-archetype bulk assignment、demo ranking 的结论。
7. 无 `PD-1 anchor -> PD-1+X repair` 逻辑链的组合提名。

---

## 5) Evidence gating by step
- L1 必须至少满足：S0/S1/S2/S3 完整、跨队列稳健性与 baseline 对照可复现。
- L2 必须在 L1 基础上新增：S4 机制支持 + S5 空间/组织裁判。
- L3 必须在 L1/L2 基础上新增：S6 外部锚点一致性 + druggability + minimal wet-lab feasibility。

---

## 6) Compliance checklist (pre-release)
- [ ] 每个 Step 均完成 manifest 记录。
- [ ] 主结论中的每一条 claim 均可映射到具体中间表。
- [ ] 所有 nomination 均包含 PD-1 anchor 对应字段。
- [ ] 所有 L2/L3 claim 均附 spatial adjudication 证据。
- [ ] 所有高级模型结论均附强基线对照与负对照审计。
