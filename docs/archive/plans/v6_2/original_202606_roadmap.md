# v6.2 迭代施工 Roadmap（复用导向 + 可迭代提升点标注）

> 本 roadmap 用于在 **v6.1 已完成成果（Phase1–10）+ 2026-06-11 增补层 + `data/` 本地资产** 基础上**增量迭代**到 v6.2。
> 配套：科学定义见 `docs/plan/v6_2/plan.md`（本文件不重复科学论证）；资产路径以 `results/v6_1/INDEX.md`（2026-06-18 实扫版）为准。
> 阅读约定：每步给 **目的 / 复用（现有代码·数据）/ 新建·改造（建议放 `scripts/v6_2/`）/ 交付 / 门 / 🔬迭代提升点**。
> **🔬 = 可通过算法或参数迭代获得性能提升的环节**；每个 🔬 标注「调什么 · 提升哪个指标 · 如何评测」，是方法学创新的主战场。
> 每步内容已细化到可直接据以撰写该 phase 的施工 prompt。

---

## A. 当前起点（增量基线）

| 资产 | 状态 | v6.2 角色 |
|---|---|---|
| Step1 registry / metadata / splits | 完成 | 直接复用，仅增量登记新队列 |
| Step2 hotfix 特征矩阵 | 完成 | 直接复用为主输入 |
| Step3 可解释基线 | 完成（强 ML 回滚） | 复用 + 补 scVI latent |
| Phase4 模块 + `phase4_7_module_scores_by_sample.csv` | 完成 | **= PASCAR Layer A/B 输入矩阵** |
| Phase4 `module_displacement_if_paired.csv` | 完成 | **反事实层位移雏形** |
| Phase5 机制组 + PD1X/shared/HCC 裁定表 | 完成 | 因果先验 + 机制节点定义 |
| Phase6 扰动/靶点预筛 | 完成 | perturb 先验 + 靶点候选 |
| Phase7 IMbrave150 投影 + GSE238264 规则空间 + 候选证据矩阵 | 完成（规则级） | 升级为定量 niche + L1/L2/L3 |
| 增补层 TCR 表 / 救援 pseudobulk / GSE238264 原始 | 暂存（未并入） | 回灌主线 |
| `data/imm/{GSE123813,GSE176021,GSE149614,checkmate}` | **本地** | anchor 救援 + HCC 深描，无需下载 |
| `data/combo/{GSE200996,GSE235863,GSE207422,bi_2021_rcc}` | **本地** | anchor / PD1X / 外锚 |
| `data/db/{tcga,gtex,icgc,l1000}` | **本地** | 外部 bulk 锚（接入即可） |
| `data/perturb_seq/{GSE90063,133344,193736,306429,XAtlas}` | **本地** | perturb 先验矩阵 |
| `data/{tcr,lr_interaction,tf_regulon}` | **本地** | TCR master / niche LRG / TF |

> **复用原则**：v6.2 不重跑 Step1–Phase7 的稳定产物；只在「机制内核」（PASCAR A/B/C）与「数据回灌/激活」处增量。任何 phase 的施工 prompt 都应先 `读取并校验上游产物`，再 `仅产出本步增量`。

---

## B. 跨步骤公共件（先建，所有 🔬 迭代依赖它）

### B1. 统一评测器 `scripts/v6_2/eval_harness.py`（**P0，迭代的前提**）
- 目的：让每个 🔬 点的"改进"可量化、可复现、可比较——这是 v6.1 缺失项，也是顶刊 robustness 卖点。
- 功能：统一 split（复用 `step1/frozen_patient_split_v6_1.csv`）、LOCO / leave-center-out、置换零分布、校准（calibration）、自助置信区间、负对照；输出标准 `eval_record.yaml`（含 metric、seed、参数哈希）。
- 交付：`scripts/v6_2/eval_harness.py` + `results/v6_2/eval/registry.csv`（每次迭代一行，便于横比）。
- 🔬 迭代提升点：评测指标集本身可扩展（加 decision curve、niche-FDR、ITE 校准），但每次扩展须记录到 registry。

### B2. 迭代登记 `results/v6_2/ITERATION_LOG.md`
- 每个 🔬 点的每次尝试记：基线值 → 新值 → 改了什么（算法/参数）→ 是否采纳。防止"调了但说不清提升多少"。

---

## Step 0：控制层 / 范围 / 双轨协议

- **目的**：锁 v6.2 scope、PASCAR 三层接口契约、双轨判决阈值。
- **复用**：`docs/v6_1/{scope_lock,analysis_contract,terminology_dictionary}.md`（增量改写，非重写）。
- **新建**：`docs/v6_2/{scope_lock_v6_2.md, analysis_contract_v6_2.md, dual_track_decision_protocol.md, output_manifest_template.yaml}`；术语新增 module-ITE / counterfactual repair / spatial niche / triangulation / ISP。
- **交付**：上述 5 文件 + `scripts/v6_2/` 目录骨架 + B1/B2 公共件。
- **门 G0**：契约与双轨协议锁定，B1 评测器可跑通空样例。
- 🔬：无（治理步）。

---

## Step 1：数据资产激活（**重在复用，不重盘点**）

- **目的**：把 `data/` 已本地资产与增补层产物登记为 v6.2 可用队列；仅对真正缺失者发下载请求。
- **复用**：
  - `results/v6_1/INDEX.md §6` 数据盘点；`step1/cohort_registry_v6_1.csv`、`dataset_role_assignment.csv`、`treatment_context_flags.csv`、`frozen_patient_split_v6_1.csv`。
  - 脚本 `scripts/v6_1/{step1_data_governance.py, register_project_cohort_package.py, build_project_data_pool.py}`（直接复用其登记逻辑）。
- **新建·改造**：`scripts/v6_2/step1b_activate_cohorts.py` —— 把 `data/imm/{GSE123813,GSE176021,GSE149614,checkmate}`、`data/combo/{GSE200996,GSE235863,GSE207422,bi_2021_rcc}`、`data/db/{tcga,gtex,icgc}`、增补层 4 数据集登记进 `cohort_registry_v6_2.csv` 并赋角色；标 access=local/need_download/controlled。
- **交付**：`cohort_registry_v6_2.csv`、`external_cohort_activation_plan_v6_2.md`、`missing_data_request_v6_2.md`（仅含 GSE215011/140901/279750/78220/91061/PRJEB23709 等 absent 项）。
- **门 G1**：每个 v6.2 角色（anchor / HCC 深描 / PD1X / 外锚 / 空间 / 扰动 / TCR）≥1 个 local 队列可用。
- 🔬：无（治理步）。

---

## Step 2：Immune-state measurement（增量：niche-ready 注释 + 回灌 + 可选 scFM）

- **目的**：复用 Step2 主矩阵；只补「niche-ready 髓系细分」「救援队列回灌」「可选 scFM 表征」。
- **复用**：
  - **主矩阵**：`step2_repair/repair_B2_unified_feature_matrix_by_sample.hotfix.parquet`（直接用，不重算）。
  - 注释/髓系/特征工厂脚本：`scripts/v6_1/{step2_5_consensus_annotation_v1.py, step2_6_myeloid_adjudication_v1.py, step2_7_cell_fraction…, step2_8_pseudobulk…, step2_9_signature_pathway_tf…, step2_9_tf_activity_rescue_v1.py, step2_10_unified_feature_matrix_v1.py}`。
  - 救援 pseudobulk：增补层 `pseudobulk_addendum.GSE120575_rescue.tsv.gz`、`…GSE229772_rescue.tsv.gz`（已建，回灌即可）。
  - `data/tf_regulon/`、`data/pathway/`。
- **新建·改造**：
  - `scripts/v6_2/step2b_niche_ready_annotation.py`：在现有 consensus 注释上**细化髓系亚型**（S100A9⁺CD14⁺单核 / TREM2⁺巨噬 / C1QC⁺ / CAF / 内皮），为 Step8 niche 映射做准备。
  - `scripts/v6_2/step2c_backfill_rescue_features.py`：把救援 pseudobulk 接 `step2_8/2_9` 流程，再生 hotfix 矩阵；移除幻影患者 `gse120575_2`。
  - （可选）`scripts/v6_2/step2d_scfm_embedding.py`：scFM/scVI embedding → `data/features/cell_embeddings/`（当前为空），配 leave-cancer 混杂审计。
- **交付**：`cell_state_annotation_v6_2.csv`（niche-ready）、回灌后的 hotfix 矩阵 v2、`myeloid_QC_report_v6_2.md`、（可选）`scfm_embedding_*.parquet` + `scfm_confounding_audit.md`。
- **门 G2**：髓系 QC 通过；救援队列无零表达；可选 scFM 须过 leave-cancer 审计方可下游。
- 🔬 **迭代提升点**：
  - 🔬 **髓系注释分辨率**（聚类 resolution、marker panel、参考映射）→ 提升 niche 信号的细胞态纯度 → 评测：注释一致性 / 与已知 marker 的 AUROC。
  - 🔬 **QC 阈值**（MAD 倍数、doublet 阈值、ambient 校正）→ 提升下游模块稳定性 → 评测：bootstrap 模块支持度。
  - 🔬 **scFM 选型与微调**（模型、是否 fine-tune、embedding 维度）→ 提升跨癌种 batch 对齐 → 评测：leave-cancer 下免疫态保真 vs 癌种泄漏（kBET/iLISI 类）。

---

## Step 3：强基线包（增量：补 scVI latent，streamline）

- **目的**：复用 Phase3.4 基线；补缺的 scVI latent simple classifier；统一进 B1 评测器。
- **复用**：`step3_qc_aware_strong_baseline/07_milestone_B_step4_handoff/{step3_baseline_results_master.csv, step3_model_evidence_grading.csv, step3_robustness_summary.csv, step3_negative_control_summary.csv}`；`04_interpretable_baselines/`；**不复用** `05_strong_ml_baselines/`（已回滚）。
- **新建·改造**：`scripts/v6_2/step3b_scvi_latent_baseline.py`（scVI latent + 简单分类器，仅作 baseline）；把全部基线重跑进 B1 评测器统一口径。
- **交付**：`baseline_results_v6_2.csv`、`baseline_comparison_report_v6_2.md`（含 scVI latent）。
- **门 G3**：六类基线齐备且同口径（fraction/signature/pathway/EN·XGB/static module/scVI latent）；后续因果/反事实须超此基线方可主叙事。
- 🔬 **迭代提升点**：
  - 🔬 **基线模型族与正则**（EN 的 α/l1_ratio、XGB 深度/学习率、scVI latent dim/层数）→ 抬高"基线天花板"（基线越强，复杂方法的增量越可信）→ 评测：LOCO AUC/AUPRC + 校准。
  - 🔬 **特征族选择**（哪些 feature family 进基线）→ 评测：去一族的 ΔAUC（特征重要性归因）。

---

## Milestone B：Baseline + Causal 合练（Step 2–5 后）
检查 immune-state → 基线 → 模块 → 因果图链路自洽；输出 `milestone_B_baseline_causal_rehearsal.md`。

---

## Step 4：Shared/HCC-specific module discovery（增量：复用 + 稳定性增强）

- **目的**：直接复用 Phase4 模块作为因果图节点；仅在稳定性/粒度上迭代。
- **复用**：`phase4_role_aware_module_discovery_20260611/09_final_report/{phase4_module_master_table.csv, phase4_phase5_handoff_main_modules.csv}`、`07_module_scoring/phase4_7_module_scores_by_sample.csv`（=主矩阵）、`02_response_blind_modules/phase4_2_module_detection_parameters.yaml`、`03_module_stability/*`、`05_module_classification/phase4_5_HCC_specific_module_candidates.csv`；脚本 `scripts/v6_1/run_phase4_role_aware_module_discovery.py`。
- **新建·改造**：仅当 Step2 回灌改变了特征时，按需重跑 Phase4 的模块评分；否则直接复用。
- **交付**：`module_score_matrix_v6_2.csv`（= 复用或回灌后再生）、`module_stability_report_v6_2.md`。
- **门 G4a**：模块矩阵 join key 与 anchor/treatment_context 对齐无泄漏。
- 🔬 **迭代提升点**：
  - 🔬 **模块发现方法与粒度**（WGCNA 软阈值 β、模块数 k、合并阈值；或换 graph-free 模块）→ 提升模块稳定性与可解释性 → 评测：bootstrap 成员支持度、cross-cohort 复现率。
  - 🔬 **模块得分算法**（AUCell vs GSVA vs 平均 z）→ 提升与 response 的关联效应量 → 评测：cohort-adjusted 效应量 + LOCO 稳定性。

---

## Step 5：PASCAR Layer A — 共享/特异因果模块图（**核心创新**）

- **目的**：把模块从关联升级为因果结构：跨癌种不变 shared 骨架 + HCC 特异残差边。
- **复用**：`module_score_matrix_v6_2.csv`（Step4）；`phase5_4_{shared_mechanism,HCC_context,PD1X_repair_logic}_adjudication.csv`（作因果先验/边方向先验）；`data/perturb_seq/*` + `phase6_2_perturbation_direction_programs.csv`（**扰动方向作边定向软约束**）；`data/tf_regulon/`（TF→target 先验边）。
- **新建·改造**（`scripts/v6_2/step5_causal_module_graph.py`）：
  - v0：模块级 graphical lasso（GGM），shared/HCC 池差集 = specific 边 + 置换零控。
  - **v1（主推）**：模块级不变因果发现——跨癌种作"环境"，ICP/IRM 风格筛 shared invariant 父集；HCC 特异边用环境交互项；扰动方向软约束定向。
  - v2（野心）：可微因果发现（NOTEARS/DCDI 变体）+ 扰动作软干预。
- **交付**：`shared_specific_causal_edges.csv`、`causal_module_graph.graphml`、`causal_module_graph_report.md`。
- **门 G4**：shared backbone LOCO 稳定；≥1 HCC 特异因果屏障可报告；因果方法在"机制解释/外部一致"上不弱于关联基线。
- 🔬 **迭代提升点（方法学主战场）**：
  - 🔬 **因果方法选择**（GGM → ICP/IRM → NOTEARS/DCDI）→ 提升因果边 precision/recall → 评测：对已知免疫调控关系（TF-regulon/扰动真值）的 precision/recall。
  - 🔬 **不变性判据与环境定义**（ICP 显著性 α、"环境"按癌种 vs 按中心、不变性阈值）→ 提升 shared/specific 分离度 → 评测：shared 边的跨癌种复现率 vs specific 边的特异性。
  - 🔬 **扰动软约束权重**（先验边定向的惩罚强度 λ_prior）→ 提升边定向准确率 → 评测：与扰动方向程序一致的边占比。
  - 🔬 **正则强度**（glasso λ、稀疏度）→ 提升图的稳定性与可解释性 → 评测：bootstrap 边支持度。

---

## Step 6：PD-1 monotherapy anchor 救援（双轨判决门）

- **目的**：用本地 anchor 队列建立 PD-1 因果敏感性/原发耐药轴；触发双轨判决。
- **复用**：
  - **本地 anchor 数据（无需下载）**：`data/imm/GSE123813`、`data/imm/GSE176021`、`data/combo/{GSE200996,GSE207422}`、`data/combo/bi_2021_rcc`。
  - Step1/Step2 治理与特征脚本（anchor 分支重用）；`step3_step4_handoff_support_PD1_anchor.csv`（现有 anchor 支持证据）；`step1/PD1_anchor_label_availability_report.md`。
- **新建·改造**：`scripts/v6_2/step6_anchor_rescue.py` —— 整合上述队列（剥离 GSE200996 联合臂）→ anchor 分支 Step1/2 → within-cohort meta-analysis → 混杂审计（复用增补层 `cohort_response_confounding_audit` 逻辑）。
- **交付**：`PD1_anchor_modules.csv`、`PD1_sensitivity_axis.csv`、`PD1_primary_resistance_axis.csv`、`PD1_anchor_patient_scores.csv`、`dual_track_decision_record.md`。
- **门 G5（双轨判决）**：加权 cohort-response 混杂 < HIGH 且 within-cohort 方向一致 → **R1**；否则 → **R2**（反事实流形切 HCC responder-like）。
- 🔬 **迭代提升点**：
  - 🔬 **anchor 队列组合与谐调**（纳入哪些癌种、是否用 scVI/scFM 跨队列谐调）→ 直接决定 G5 混杂等级 → 评测：加权 max-cohort fraction、cohort-response χ² p、within-cohort 方向一致率。
  - 🔬 **meta-analysis 加权**（固定 vs 随机效应、按样本量/质量加权）→ 提升合并方向稳健性 → 评测：留一队列后方向稳定性。
  - 🔬 **响应标签映射阈值**（RECIST/MPR 二分界）→ 评测：response-mapping 敏感性分析下结论稳定性。

---

## Step 7：PASCAR Layer B — 反事实修复 / module-ITE（**核心创新**）

- **目的**：在因果图上量化"修复模块把患者推向 responder-like 多远"，落地正式 M-signals。
- **复用**：`causal_module_graph.graphml`（Step5）；`PD1_anchor_patient_scores.csv`（Step6，R1）或 HCC responder-like（R2）；**`phase4_7_module_displacement_if_paired.csv`（已有配对位移雏形，作 ITE 初始化/校验）**；`phase6_4_{PD1X_repair_axis_summary,mechanism_intervention_axis_map}.csv`（X-class 先验）；`phase5_4_PD1X_repair_logic_adjudication.csv`。
- **新建·改造**（`scripts/v6_2/step7_counterfactual_repair.py`）：
  - v0：规则版 mIMS/MRI/mICS（保底产文件）。
  - **v1（主推）**：条件最优传输反事实（adapt CINEMA-OT / scDRP）——responder/non-responder 流形间求位移场；`do(repair_M)` → module-ITE = 沿位移到 responder-like 的可达增益 = MRI；mICS = ITE × druggability × 三角证据。
  - v2：基于 Layer A 因果图的 structural counterfactual（do-calculus 路径传播）。
- **交付**：正式 `mIMS_scores.csv`、`MRI_scores.csv`、`mICS_scores.csv`；`PD1X_repair_logic_table.csv`、`residual_barrier_modules.csv`、`X_mechanism_class_mapping.yaml`、`counterfactual_repair_report.md`。
- **门 G6**：每 nomination 连回 anchor（R1）/HCC responder-like（R2）；不凭单证据源；M-signals 文件存在。
- 🔬 **迭代提升点（方法学主战场）**：
  - 🔬 **OT 代价与正则**（cost metric、entropic-OT 的 ε、unbalanced-OT 边际松弛）→ 提升反事实匹配质量与 ITE 校准 → 评测：在配对 pre/post（`module_displacement_if_paired`）或留出扰动上的 ITE 预测误差。
  - 🔬 **隐空间/解构维度**（scDRP 式 disentangled latent 维度、是否扣除混杂方向）→ 提升 ITE 的因果可识别性 → 评测：负对照（随机标签）下 ITE→0 的程度。
  - 🔬 **do(repair) 算子定义**（模块整体平移 vs 沿因果父集传播）→ 提升与生物学先验一致性 → 评测：与扰动方向程序/已知 X(anti-VEGF 等)效应的一致率。
  - 🔬 **mICS 权重组合**（ITE/druggability/centrality/niche/external 的加权）→ 提升 Top 候选的下游验证命中率 → 评测：与 Phase7 外部/空间支持的一致性、湿验命中。

---

## Step 8：空间 niche 定量裁判（**升为共主轴**）

- **目的**：把空间从规则赋标升级为 niche 级定量；回答 plan Q4（屏障是否成 niche、X 修复是否预测 niche 解构）。
- **复用**：
  - **空间原始数据已在本地**：`new_data_integration_addendum_20260611/scratch/GSE238264_raw/GSM7661255_HCC1R … 7661261_HCC7NR.tar.gz`（7 例，文件名带 R/NR）；`spatial_sample_manifest.GSE238264.tsv`。
  - 现有规则空间结果作对照：`phase7_4_spatial_mechanism_adjudication.csv`、`phase7_4_candidate_spatial_support_matrix.csv`。
  - `data/lr_interaction/`（配体-受体）；`cell_state_annotation_v6_2.csv`（Step2 niche-ready）；模块基因来自 `phase4_2_module_feature_membership.csv`。
- **新建·改造**（`scripts/v6_2/step8_spatial_niche_adjudication.py`）：niche 构建（ISCHIA 风格局部细胞群落）+ 细胞互作图 + 模块基因/signature → spot/niche 投影 + T 细胞排斥/髓系-肿瘤邻近/APC niche/CAF 屏障/LRG 共定位 + 置换零控。
- **交付**：`spatial_niche_scores.csv`、`myeloid_tumor_neighborhood_scores.csv`、`T_cell_exclusion_scores.csv`、`APC_niche_scores.csv`、`CAF_barrier_scores.csv`、`LRG_niche_support.csv`、`spatial_niche_adjudication_report.md`。
- **门 G7**：≥1 核心因果屏障获 niche 级定量支持（预定义方向 + 零控 FDR）；否则降级。
- 🔬 **迭代提升点（方法学主战场）**：
  - 🔬 **niche 数 k 与邻域定义**（niche 聚类数、邻域半径/最近邻 K）→ 提升 niche 检出力与稳定性 → 评测：niche 跨样本复现率、轮廓系数。
  - 🔬 **互作检验与零模型**（置换策略、空间随机化方式）→ 控制 niche-FDR → 评测：零模型下假阳率。
  - 🔬 **模块→spot 投影方法**（marker 富集 vs 反卷积 RCTD vs signature 打分）→ 提升空间信号信噪比 → 评测：与组织学/已知 niche 的一致性。

---

## Step 9：扰动先验 + in-silico 扰动三角验证

- **目的**：构建 perturb-prior 矩阵；（可选）scFM ISP；与 OT 反事实三角验证；构建 TCR master。
- **复用**：
  - **本地扰动全齐**：`data/perturb_seq/{GSE90063,GSE133344,GSE193736,GSE306429,XAtlas}` + `data/db/l1000`。
  - 现有扰动映射：`phase6_2_{perturbation_mapping_scores,perturbagen_evidence_table,perturbation_direction_programs}.csv`。
  - **TCR**：增补层 `tcr_clone_addendum.GSE236581.tsv`（1.1M 记录）+ `data/tcr/` + anchor 队列自带 scTCR（GSE123813/176021/200996）。
- **新建·改造**：`scripts/v6_2/step9a_perturb_prior_matrix.py`（本地扰动→矩阵）；（可选）`step9b_scfm_isp.py`（scFM in-silico 扰动）；`step9c_tcr_master.py`（→`tcr_clone_master_v6_2.csv`，clonal expansion 连耗竭/细胞毒模块）；`step9d_triangulation.py`（OT vs scFM vs prior 方向收敛）。
- **交付**：`perturbation_prior_matrix.csv`、`perturbation_isp_triangulation.csv`、`tcr_clone_master_v6_2.csv`、（可选）`isp_benchmark_report.md`。
- **门 G8**：scFM ISP 若不收敛/不稳→退可选佐证，三角化以 OT 反事实 + 扰动先验为准。
- 🔬 **迭代提升点**：
  - 🔬 **perturb-prior 聚合**（按细胞类型/剂量加权、方向显著性阈值）→ 提升先验方向可靠性 → 评测：与已知调控方向一致率。
  - 🔬 **scFM ISP 设置**（模型、fine-tune、扰动算子）→ 提升 ISP-vs-真值/与 OT 收敛度 → 评测：在已测扰动上的方向 AUROC、三引擎一致率。
  - 🔬 **三角化聚合规则**（≥2 一致的加权 vs 投票、冲突惩罚）→ 提升最终 nomination 的稳健性 → 评测：与外部/空间支持的吻合度。

---

## Milestone C：Mechanism 合练 + 双轨头条裁决（Step 6–9 后）
首次跑通 PASCAR 链（因果图→anchor→反事实→niche→三角雏形）。**依 G5 + HCC/泛癌信号 + 时间锁头条 R1/R2/R3**。输出 `milestone_C_mechanism_rehearsal.md`、`headline_route_decision.md`。

---

## Step 10：PASCAR Layer C — 多模态证据图三角化 + L1/L2/L3 矩阵（**核心创新**）

- **目的**：把所有证据重组为可裁判 L1/L2/L3 + 三角化把关。
- **复用**：`phase7_5_candidate_evidence_matrix.csv`（已有候选证据矩阵，**升级为 L1/L2/L3**）；`phase5_2_{mechanism_evidence_matrix,support_conflict_matrix}.csv`；各层产物（因果/ITE/niche/perturb-ISP/外锚）。
- **新建·改造**（`scripts/v6_2/step10_evidence_graph.py`）：v0 加权聚合+置换标定；**v1（主推）正交证据贝叶斯融合**（各模态条件独立证据→后验+CI，置换/null 校准似然）+ 三角化规则；v2（可选）异质 GNN 仅作敏感性对照。
- **交付**：`candidate_module_evidence_matrix.csv`（L1/L2/L3 + triangulation count + 降级规则）、`evidence_graph_report.md`、`conflict_downgrade_log.md`。
- **门 G9**：每 candidate 有显式 L1/L2/L3 + 三角化命中数；关键冲突已降级记录。
- 🔬 **迭代提升点**：
  - 🔬 **证据权重/先验**（各模态似然权重、先验概率）→ 提升 L 等级与已知生物学/留出真值的吻合 → 评测：对正/负对照机制的 L 等级 precision。
  - 🔬 **似然校准**（置换 null 拟合方式）→ 提升后验校准 → 评测：calibration（可靠性曲线）。
  - 🔬 **三角化阈值**（升 L2/L3 需几条一致、冲突如何降级）→ 平衡敏感性/特异性 → 评测：保留候选的下游验证命中率 vs 数量。

---

## Step 11：Bulk/external anchor 激活（**本地 TCGA/GTEx/ICGC 直接接入**）

- **目的**：把 bulk 外锚从趋势图升级为模块级投影验证。
- **复用**：
  - **本地外锚**：`data/db/{tcga(26G),gtex,icgc}`、`data/combo/{IMbrave150,GSE235863}`、`data/imm/checkmate`。
  - 现有投影：`phase7_2_{external_projection_scores,candidate_external_support_matrix,external_conflict_log}.csv` + 脚本 `run_phase7_validation_package.py`（升级其投影为模块级+CI）。
- **新建·改造**：`scripts/v6_2/step11_external_anchor.py`（模块级投影 + 置信区间 + 失败案例；HCC 用 TCGA-LIHC/ICGC-LIRI + IMbrave150 + checkmate；需下载者见 Step1 缺失清单）。
- **交付**：`bulk_module_projection_scores.csv`、`external_anchor_results.csv`、`external_consistency_flags.csv`、`HCC_external_support_report.md`。
- **门 G10**：核心机制在 ≥1 独立外锚复现方向。
- 🔬 **迭代提升点**：
  - 🔬 **投影方法**（top-gene 均值 → 模块级 ssGSEA/反卷积 → 生存/response 模型）→ 提升外锚一致性灵敏度 → 评测：方向一致率 + 生存 HR 显著性。
  - 🔬 **批次/纯度校正**（TCGA bulk 的肿瘤纯度、批次）→ 降低假阳 → 评测：纯度敏感性分析。

---

## Milestone D：Evidence 合练（Step 10–11 后）
为每候选定 L1/L2/L3，判进主文/补充。输出 `milestone_D_evidence_rehearsal.md`、`candidate_drop_or_keep_decision.md`。

---

## Step 12：Evidence cards + nomination + 最小湿验

- **目的**：整理 evidence cards 与 PD-1+X nomination（非临床推荐），设计最小湿验。
- **复用**：`phase7` 证据卡 + `phase8/05_validation_execution_package` + `phase6_3_target_candidate_registry.csv`（靶点）。
- **新建·改造**：`scripts/v6_2/step12_nomination_cards.py`（14 字段证据卡，见 plan §8）；`minimal_wetlab_validation_plan_v6_2.md`（IHC/mIF 髓系 niche → 共培养/极化 → CRISPR/小分子，按平台可得性）。
- **交付**：`PD1X_nomination_cards.csv`、`final_evidence_table.csv`、`top_{modules,targets}_summary.md`、`minimal_wetlab_validation_plan_v6_2.md`。
- **门 G11**：每 nomination 连 anchor/HCC；不凭单证据源；湿验可执行。
- 🔬：无（整合步）；湿验靶点选择隐含 mICS 权重迭代（见 Step7 🔬）。

---

## Step 13：论文主图 + 复现包（闭合 v6.1 三阻塞项）

- **目的**：压成 6 主图主线 + 投稿就绪复现包。
- **复用**：`phase8/04_figure_construction_plan`、`phase9/03_main_figure_draft_package`、`phase10/{02_main_figure_production,04_manuscript_polishing}`（图模板与手稿骨架直接演进）。
- **新建·改造**：`scripts/v6_2/step13_figures_and_repro.py`；补 Supp.S1（精确计数）、Phase4/Step3 参数附录、环境锁（含 scFM/OT/因果方法版本与种子）。
- **交付**：`figure_plan_v6_2.md`（Fig3 因果图 / Fig4 niche / Fig5 反事实为重头）、`main_results_narrative.md`、`reproducibility_manifest_v6_2.md`、`table_figure_source_mapping.csv`。
- **门 G12**：external_submission_ready=true（三阻塞项闭合）。
- 🔬：无（出版步）。

---

## Milestone E：Paper Story 合练（Step 12–13 后）
合成最终故事线，定主路线 R1/R2/R3、六图结论、主文 vs 补充、降级 claim、下一轮三优先、湿验优先级、20/40 分钟汇报大纲。输出 `milestone_E_*`。

---

## C. 🔬 迭代提升点总览（方法学创新主战场速查）

| 步骤 | 🔬 提升点 | 性能指标 | 评测器手段 |
|---|---|---|---|
| Step2 | 髓系分辨率 / QC 阈值 / scFM 选型 | 注释纯度 / 模块稳定性 / 跨癌种对齐 | annotation AUROC / bootstrap / kBET-iLISI |
| Step3 | 基线模型族·正则 / 特征族 | 基线天花板 AUC | LOCO + 校准 |
| Step4 | 模块方法·粒度 / 打分算法 | 模块稳定性·效应量 | bootstrap 支持度 / cohort-adjusted ES |
| **Step5** | **因果方法 / 不变性判据 / 扰动约束权重 / 正则** | **因果边 P/R·shared-specific 分离** | **vs TF-regulon/扰动真值 P/R, 跨癌种复现率** |
| Step6 | anchor 组合·谐调 / meta 加权 / 标签阈值 | G5 混杂等级 | χ²p, max-cohort frac, 方向一致率 |
| **Step7** | **OT 代价·正则 / 解构维度 / do 算子 / mICS 权重** | **ITE 校准·因果可识别性** | **配对位移/留出扰动 ITE 误差, 负对照 ITE→0** |
| **Step8** | **niche k·邻域 / 零模型 / 投影方法** | **niche 检出力·FDR·信噪比** | **复现率, 零模型假阳率, 组织学一致** |
| Step9 | perturb 聚合 / scFM ISP / 三角规则 | 先验方向可靠·三引擎一致 | 方向 AUROC, 一致率 |
| **Step10** | **证据权重·先验 / 似然校准 / 三角阈值** | **L 等级 precision·后验校准** | **正负对照 precision, 可靠性曲线** |
| Step11 | 投影方法 / 批次纯度校正 | 外锚一致性 | 方向一致率, 生存 HR, 纯度敏感性 |

> **加粗的 Step5/7/8/10 是顶刊方法学增量的四个核心可迭代环节**；建议每个都先建 v0 基线值入 `ITERATION_LOG.md`，再逐轮迭代 v1/v2，用 B1 评测器量化每次提升。

---

## D. 压缩版施工顺序

0 控制层 → 1 资产激活 → 2 immune-state(+回灌/niche-ready) → 3 基线(+scVI) → **MsB** → 4 模块复用 → **5 Layer A 因果图🔬** → 6 anchor 救援(双轨门 G5) → **7 Layer B 反事实🔬** → **8 niche🔬** → 9 扰动+ISP三角 → **MsC(锁头条)** → **10 Layer C 证据图🔬** → 11 外锚激活 → **MsD** → 12 nomination+湿验 → 13 主图+复现 → **MsE**

---

## E. 与 v6.1 roadmap 的关系
v6.1 roadmap 为"从零 11 步"；本版为"在 v6.1 成果 + 本地资产上的增量 14 步"，差异：①每步显式列复用产物与脚本；②机制层替换为 PASCAR A/B/C 三创新层并标 🔬 迭代点；③空间升共主轴、外锚用本地 TCGA/GTEx/ICGC；④anchor 用本地 GSE123813/176021/200996/207422 救援；⑤新增 B1 统一评测器使迭代可量化。
