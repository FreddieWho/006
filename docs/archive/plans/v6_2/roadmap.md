# v6.2 Agent-orchestrated Roadmap

> 当前执行模式（2026-08-27）：按探索性研究推进；旧的预注册式阈值、冻结 gate 和 response lock 仅保留作历史溯源，不约束新的探索分支。详见 `docs/v6_2/EXPLORATORY_MODE_DIRECTIVE_20260827.md`。

# PD-1-anchored immune failure mechanism and X-class repair framework

中文名：

# v6.2 施工路线图：以 PD-1 为锚的免疫治疗失败机制与 X-class 修复框架

---

## 0. Roadmap 定位

本文件是 v6.2 的施工路线图。它不重复 scientific plan 的全部科学论证，也不包含具体 phase prompt。它定义：

1. 阶段顺序；
2. 各 lane 的并行关系；
3. 每阶段核心输入、任务、输出；
4. Go / No-Go / downgrade gate；
5. 哪些地方可用 AutoResearch-style 迭代升级；
6. 每个关键技术节点允许探索的技术路线范围；
7. 哪些结果可进入 Result Layer A/B/C；
8. 何时启动独立 agent / session；
9. 何时需要 integration rehearsal 和 hostile reviewer audit。

Revision note, 2026-06-24:
added explicit Spatial Data Standardization Gate to make existing Phase 2 / Phase 11 spatial requirements executable without changing scientific scope.

Revision note, 2026-06-24 (2) — 生物学补强挂钩（依据 `results/v6_2/PLAN_ROADMAP_BIOLOGICAL_REVIEW_20260624.md`，用户批准）：
已并入 B1 + B2；B3–B6 落档为补强审阅挂钩，进行到对应 phase 时**必须唤起**该挂钩做补强审阅。

| 挂钩 | 内容 | 唤起 phase | 是否已并入 |
|---|---|---|---|
| B1 | HCC 病因分层：**可插拔、低权重、HCC-context 限定**；默认不进 pan-cancer 主线、不主导 IO 重点；覆盖不足则关闭 | Phase 2（审覆盖）、Phase 4（打标签）、Phase 9（HCC residual sensitivity 协变量） | 是 |
| B2 | 模块本体补强：TLS/B、neutrophil-NET、tumor-intrinsic WNT exclusion、progenitor/terminal exhaustion（仍 response-blind） | Phase 4（refinement）、Phase 6（注释/稳定性雷达）、Phase 11（TLS/APC-T niche） | 是 |
| B3 | 当代 X-class + 临床失败组合（LEAP-002、COSMIC-312）作 held-out 负样本；anti-TIGIT 作正向用例 | Phase 12、Phase 13 | 否（挂钩） |
| B4 | primary / adaptive / acquired resistance 分型（数据允许时） | Phase 9、Phase 10 | 否（挂钩） |
| B5 | TCR 克隆替换 vs 再激活 作 responder 签名（有 paired TCR 时前移为推荐项） | Phase 9 | 否（挂钩） |
| B6 | 显式声明 perturbation 先验对 myeloid/stromal/TLS/中性粒的覆盖偏倚（缺支持≠证据为负） | Phase 12、Phase 14 | 否（挂钩） |

本 roadmap 的基本原则是：

**先让数据、模块和评测系统可靠，再训练 SRB；先证明 dominant barrier 可辨识，再做 repair；先有 spatial/perturbation/external 的证据接口，再让 ZSL 排序；先形成 evidence card，再决定论文 claim。**

---

## 1. 总体阶段结构

v6.2 采用三阶段结构。

### Stage 1：Rescue / Prebuild

目标：确认输入可信、资产可用、MVP 污染切断、关键数据覆盖足够，建立所有后续迭代依赖的评测和冻结机制。

输出目标：

* source-of-truth 锁定；
* data go/no-go；
* frozen input family；
* evaluation harness；
* iteration log；
* immune-state feature contract；
* anchor / spatial / perturbation / external bulk 覆盖审计；
* Stage 2 是否可启动的 decision。

### Stage 2：Core Build

目标：完成 v6.2 的核心机制模型和证据链第一轮闭环。

输出目标：

* response-blind modules；
* barrier identifiability gate；
* SRB Stage-A latent；
* anti-PD1 anchor verdict；
* shared/context-modulated modules；
* SRB transport / CRTN repair direction；
* SNGM spatial niche support；
* perturbation coverage gate；
* X-class descriptor ontology；
* ZSL alignment first-pass；
* 至少一个无关键冲突的 candidate evidence card v0。

Stage 2 的硬底线：

**如果无法形成至少一个无关键冲突的 candidate evidence card v0，不进入 Stage 3。**

### Stage 3：Manuscript Package

目标：将可防守的 candidate 组织成论文主线、主图、补图、证据卡、复现包和 claim boundary。

输出目标：

* Result Layer A/B/C 中至少一个层级完整成立；
* 主图故事；
* candidate evidence cards；
* external / spatial / perturbation validation package；
* hostile reviewer audit；
* manuscript-ready claim boundary；
* minimal reproducibility package。

---

## 2. 工作 lane 设计

v6.2 不按单线程 pipeline 执行，而按 lane 并行推进。

### Lane A：Control / Governance Lane

角色：Controller agent 常驻。

负责：

* source-of-truth；
* frozen input；
* branch agent 启动/关闭；
* handoff 审查；
* gate verdict；
* integration rehearsal；
* hostile reviewer audit；
* 关键科学决策提交给用户确认。

### Lane B：Asset / Data Lane

角色：data asset agent。

负责：

* cohort registry；
* metadata；
* response label environment；
* data go/no-go；
* anchor/spatial/perturbation/external bulk/TCR 覆盖审计；
* patch proposal。

### Lane C：Immune-state / Module Lane

角色：immune-state agent + module agent。

负责：

* immune-state construction；
* feature family audit；
* response-blind module discovery；
* module annotation；
* barrier identifiability gate。

### Lane D：SRB / ML Core Lane

角色：SRB agent。

负责：

* module-factorized latent；
* SRB Stage-A pretraining；
* context head；
* anchor head；
* transport head；
* alignment head；
* ablation；
* calibration。

### Lane E：Spatial Lane

角色：spatial asset / SNGM agent。

负责：

* spatial asset audit；
* niche construction；
* statistical niche baseline；
* SNGM if data supports；
* spatial support/conflict verdict。

### Lane F：Perturbation / X-class / ZSL Lane

角色：perturbation-ZSL agent。

负责：

* X-class descriptor ontology；
* perturbation coverage gate；
* target / pathway / drug family mapping；
* Combination ZSL；
* Mechanism ZSL；
* ZSL leakage / abstention audit。

### Lane G：Evidence / Manuscript Lane

角色：integration / manuscript agent。

负责：

* evidence cards；
* conflict log；
* claim tier；
* result layer assembly；
* figure source tables；
* manuscript skeleton；
* hostile reviewer audit.

---

## 3. 公共基础设施

### Phase 0.1：Project Control Brief

目的：

建立所有 agent / session 共同遵守的轻量项目宪法。

核心内容：

* v6.2 project identity；
* 当前 scientific plan 为最高科学依据；
* MVP-derived 结果不继承，只能作为 raw data locator / failure audit；
* HCC 是 pre-specified index cancer，不是典型癌种；
* ZSL 是 formal extrapolation layer，不是机制证明；
* SRB 是主模型，SNGM 独立空间验证；
* graph attention 不等于 causal edge；
* perturbation prior 是 soft condition + coverage gate；
* no clinical recommendation；
* branch agent 不可修改 frozen input；
* spike 默认不能自动升主线。

交付：

* `AGENT.md`
* `source_of_truth_and_no_mvp_inheritance.md`
* `claim_boundary.md`

Gate：

* AGENT.md 必须短而硬；
* 不写具体 prompt；
* 不写复杂局部 AGENT.md；
* 能被所有 branch agent 读懂。

Agent 编排：

* Controller agent 执行；
* 用户确认后冻结。

---

### Phase 0.2：Evaluation Harness

目的：

建立统一评测系统，防止每个 agent 用自己的指标宣布胜利。没有这个，AutoResearch-style 迭代会变成随机调参祭祀。

核心功能：

* patient-level split；
* leave-one-cohort；
* leave-one-cancer；
* leave-center-out；
* response-label environment split；
* permutation null；
* negative control；
* bootstrap CI；
* calibration；
* abstention metric；
* external direction consistency；
* spatial permutation null；
* perturbation direction consistency；
* model ablation registry。

交付：

* `eval_harness_spec.md`
* `eval_registry.csv`
* `iteration_log_template.md`
* `negative_control_catalog.md`

Gate：

* 所有后续模型必须向 eval harness 写入结果；
* 任何模型如果没有 baseline 和 negative control，不允许进入 evidence card；
* 所有迭代必须记录 seed、input version、parameter hash、metric、adopt/reject decision。

[AUTO-RESEARCH 入口]

从这一阶段开始，所有可迭代技术点必须使用同一循环：

1. 提出可证伪假设；
2. 生成一个小改动；
3. 在固定 frozen input + eval harness 上运行；
4. 写入 iteration log；
5. 与 v0 baseline 比较；
6. 只采纳能提升预注册指标且不增加关键冲突的改动。

技术探索范围：

主路线：

* 自实现统一评测器；
* 指标和 split 固定；
* 每次迭代显式登记。

降级路线：

* 如果完整 eval harness 过重，先实现 core split + bootstrap + permutation + calibration。

禁止：

* 临时换 split；
* 调完才定义指标；
* 只报告最佳 seed；
* 不记录失败实验。

---

## 4. Stage 1：Rescue / Prebuild

### Phase 1：Source-of-truth and Inheritance Lock

目的：

防止旧方案、MVP、v6.1 偏移结果和新 v6.2 方案互相污染。

任务：

1. 汇总当前可用 plan、review、roadmap、ZSL 报告和 v6.1 完成结果；
2. 标记哪些文件是 source-of-truth，哪些只是背景；
3. 明确 v6.1 资产需要重新审计后才能进入 v6.2；
4. 明确 MVP-derived feature、archetype、heuristic MOA、virtual shift、demo counterfactual 禁止进入 v6.2；
5. 明确可接受 raw locator / failure audit 的范围。

交付：

* `source_of_truth_matrix.md`
* `inheritance_policy.md`
* `legacy_asset_quarantine_list.csv`

Gate：

* 所有后续 phase 只能读取 approved frozen asset；
* 旧文件中与当前 v6.2 plan 冲突者自动降级为 historical context；
* MVP 任何 derived result 不得作为主输入。

Agent：

* Controller agent；
* 如文件很多，可启动 asset-audit branch agent；
* 合并后由 Controller 冻结。

技术探索范围：

无模型探索。这里要保守、清晰、可审计。

---

### Phase 2：Data Go / No-Go Audit

目的：

在模型训练前确认关键数据是否足以支撑 SRB、CRTN、SNGM、ZSL 和外部验证。这个阶段不做漂亮图，只回答“能不能做”。

审计对象：

1. clean anti-PD1 / ICI anchor；
2. expanded / support anchor；
3. paired pre/post ICI data；
4. observed PD1+X / combination data；
5. HCC index cancer data；
6. spatial / tissue-niche data；
7. perturbation / target prior data；
8. external bulk / clinical anchor；
9. TCR support；
10. response-label environment metadata。

任务：

* 样本量审计；
* response label 质量审计；
* treatment context 审计；
* patient/sample duplication 审计；
* timepoint 审计；
* platform/cohort/batch 审计；
* X-class coverage 审计；
* spatial metadata 审计；
* perturbation coverage 审计；
* （B1 挂钩）HCC etiology 元数据覆盖审计（HBV/HCV/alcohol/NASH-MAFLD/其他）；仅决定可插拔低权重 sensitivity 层能否启用，覆盖不足则关闭并记缺口，不影响主线；
* 数据缺口与可补数据列表。

交付：

* `data_go_no_go_report.md`
* `cohort_registry_v6_2_candidate.csv`
* `response_label_environment_table.csv`
* `data_coverage_matrix.csv`
* `data_gap_register.csv`
* `stage1_data_gate.yaml`

Gate：

* 若 clean anchor 不足，R1 降级为 R1-lite 或 R2；
* 若 paired displacement 不足，CRTN 最高档禁用；
* 若 spatial 数据不足，SNGM 降级为 statistical niche adjudication；
* 若 perturbation coverage 低，ZSL 和 transport 中 perturbation loss 降权；
* 若 response label environment 不清，所有 supervised claim 降级。

Agent：

* Data asset agent 独立执行；
* Controller 审核并冻结 `data_coverage_matrix`；
* 高风险数据补充需要用户确认。

[AUTO-RESEARCH 标注]

这里不是算法迭代，而是可进行 **data acquisition search iteration**：

* 检索新增 spatial HCC ICI 数据；
* 检索 PD1+X paired / pre-post 数据；
* 检索 perturb-seq / LINCS target coverage；
* 每次新增数据必须通过同一 go/no-go schema。

技术探索范围：

主路线：

* 只补能直接修复 evidence gap 的数据。

并行探索：

* spatial 补强；
* PD1+X observed combination 补强；
* external bulk HCC ICI 补强；
* perturbation coverage 补强。

禁止：

* 为了增加样本数引入 metadata 不明的高混杂数据；
* 把高混杂联合治疗并入 clean anchor；
* 因数据稀缺而修改科学问题。

---

### Phase 2.5：Spatial Data Standardization Gate

目的：

把 downloaded spatial assets 转成可审计、可复现、可降级的 v6.2 spatial input family。该阶段只定义空间数据能否进入 Phase 11，不改变 SRB / SNGM / Evidence Card 的科学边界。

任务：

1. 建立 raw archive / extracted file 索引；
2. 按 Visium、Xenium、MERFISH、CosMx、GeoMx / region-level 等平台分类；
3. 建立 spatial sample metadata 最小表；
4. 对齐 expression matrix、coordinates、image、scale factor 和 sample ID；
5. 标记 patient、tissue、cancer、treatment、pre/post、response / failure label 的可用性；
6. 执行样本级 QC；
7. 记录 platform-specific fallback；
8. 冻结可进入 Phase 11 的 standardized spatial input。

交付：

* `spatial_raw_asset_manifest.csv`
* `spatial_sample_metadata_minimal.csv`
* `spatial_standardization_manifest.csv`
* `spatial_qc_report.csv`
* `spatial_standardization_failure_log.md`

Gate：

* 缺少 expression + coordinates 的样本不能进入 Phase 11；
* metadata 不足的样本只能作为 localization support；
* 缺少 response / treatment label 的样本不能支持 response / failure claim；
* image / coordinate 不齐或空间质量不足的样本不能训练 SNGM；
* GeoMx / region-level 数据进入 region-localization fallback，不作为 full spatial graph；
* 任何 raw download 目录不得绕过本 gate 直接进入 Phase 11。

Agent：

* Spatial asset agent 执行；
* Controller 审核并冻结 `spatial_standardization_manifest.csv`；
* SNGM agent 只能读取通过本 gate 的 standardized spatial input。

技术探索范围：

主路线：

* 保守解包、最小 schema、样本级 QC、显式降级。

可靠降级：

* 只生成 raw asset registry 和 localization support list；
* 对无法标准化的数据只进入 failure log。

禁止：

* 为了训练 SNGM 强行补齐缺失 metadata；
* 把 image / coordinate 不齐的样本作为 graph evidence；
* 把没有 response / treatment label 的空间数据写成 failure evidence；
* 让空间标准化结果改变 v6.2 的科学问题。

---

### Phase 3：Frozen Input Family

目的：

把 Stage 2 所需输入冻结，防止多 agent 并行时各自用不同 metadata / feature matrix。

任务：

1. 生成或审计 cohort registry；
2. 生成 sample/patient metadata master；
3. 建立 response-label environment；
4. 建立 canonical patient split；
5. 建立 feature family eligibility；
6. 建立 raw-to-feature provenance；
7. 标记 frozen input version。

交付：

* `cohort_registry.frozen_v0.csv`
* `patient_metadata_master.frozen_v0.csv`
* `sample_metadata_master.frozen_v0.csv`
* `patient_split.frozen_v0.csv`
* `feature_family_eligibility.frozen_v0.csv`
* `frozen_input_manifest.yaml`

Gate：

* branch agent 只能读取 frozen input；
* 修改必须提交 `patch_proposal`；
* Controller 才能发布 frozen_v0.1 / v1.0；
* 所有输出必须声明 input version。

Agent：

* Controller + data asset agent；
* 后续所有 branch agent 启动前必须读取 frozen manifest。

技术探索范围：

无高级模型。只允许 schema 和 QC 改进。别在这里发明宇宙，宇宙已经够乱了。

---

---

### Phase 3.5：Single-cell Processing QC Gate

目的：

把 Phase3 冻结的 scRNA / scTCR / processed object / rescue expression 资产转成可审计、可追溯、可进入整合的标准分析对象。该阶段不做生物学结论，不做模块发现，不做 supervised modeling，不做完整 atlas claim。它只回答一个问题：

**这些单细胞对象是否足够干净、可追踪、可整合，能否安全进入 Phase4A / Phase4B。**

本阶段是 Phase3 Frozen Input Family 与 Phase4 immune-state construction 之间的工程质量门。Phase3 只冻结 metadata、role、split、manifest 和使用边界；Phase3.5 冻结 expression object、QC、HVG、layer decision、barcode/sample/patient 对齐和 integration readiness。

输入：

* `frozen_input_manifest.yaml`
* `cohort_registry.frozen_v0.csv`
* `sample_metadata_master.frozen_v0.csv`
* `patient_metadata_master.frozen_v0.csv`
* `patient_split.frozen_v0.csv`
* `response_label_environment.frozen_v0.csv`
* `dataset_role_and_feature_eligibility.frozen_v0.csv`
* `label_inventory_and_semantic_review.csv`
* Phase3 formal intake / spatial gate / expression rescue ledger
* downloaded or local scRNA objects：raw count、processed Seurat RDS、h5ad、MTX、H5、metadata、TCR files、normalized expression objects

任务：

1. 建立 raw / processed object processing manifest。
   对每个 cohort 记录原始文件、对象类型、路径、hash、assay/layer、metadata 来源、是否可读取、是否与 `frozen_v0` 对齐。

2. 建立 matrix layer decision。
   对每个 cohort 明确后续使用哪一层：raw counts、filtered raw counts、normalized expression、processed Seurat assay、h5ad layer、pseudobulk rescue、support-only matrix。
   禁止 raw count 和 normalized matrix 混用而不记录。

3. 统一 gene identifier。
   完成 Ensembl ID / gene symbol / duplicated gene symbol / mitochondrial-ribosomal-stress genes 的标准化记录。
   输出 gene overlap、gene loss、duplicated gene handling 和 immune-core gene coverage。

4. 对齐 sample / patient / cell barcode。
   确保每个 cell 可以追溯到 `cohort_id`、`sample_id`、`patient_id`、`sample_key`、`patient_key`。
   对多 library / multi-timepoint / sorted population / pre-post 数据建立 library-to-sample mapping。
   不允许 cell barcode 离开 sample context 后直接合并。

5. 执行 sample-level 和 cell-level QC。
   至少记录每个样本的 cell 数、gene 数、UMI 数、mitochondrial fraction、ribosomal fraction、hemoglobin / platelet / dissociation / stress signal、low-quality cell fraction、empty / ambient RNA 风险、doublet 风险。
   允许使用已有 QC-pass object，但必须记录其 QC 来源和不可重复风险。

6. 建立 filtering decision log。
   对每个 cohort / sample 记录过滤前后细胞数、过滤阈值、过滤原因、是否采用原作者 QC、是否追加本项目 QC。
   低细胞数、低基因覆盖、严重 batch 或 metadata 不完整的样本进入 downgraded / support-only，不得悄悄删除。

7. 建立 TCR / GEX join QC。
   对有 TCR 的 cohort 记录 GEX-TCR join key、join 率、matched cell 数、unmatched TCR 数、unmatched GEX 数、clonotype 字段可用性。
   TCR join 不清楚的 cohort 不进入 TCR feature 主线，只能作为 support。

8. 建立 cell annotation source inventory。
   记录每个 cohort 的原始 celltype 字段、annotation 层级、是否有 immune / stromal / tumor coarse label、是否有 myeloid / T / NK / CAF / endothelial 细分、annotation 来源、可信度和是否需要 Phase4A 重新映射。

9. 冻结 HVG strategy。
   HVG 选择必须 response-blind。
   至少输出三类 HVG：

   * pan-cancer global HVG；
   * major-cell-state-stratified HVG；
   * caution / exclusion HVG list。
     caution list 应标记 mitochondrial、ribosomal、hemoglobin、stress、cell-cycle、dissociation、platform-dominated、single-cohort-dominated genes。
     HVG 不得按 responder / non-responder 选择。

10. 执行 integration readiness audit。
    对每个 cohort 和每个 major cell compartment 判断是否可进入 full integration、coarse integration、reference mapping only、pseudobulk only、support-only 或 excluded。
    审计内容包括 gene overlap、cell count、sample count、batch covariates、platform、sample composition、annotation availability、raw/normalized layer clarity。

11. 建立 batch covariate registry。
    冻结后续 integration 可用的 covariates：cohort、study、platform、chemistry、sample、patient、cancer type、tissue、timepoint、treatment context、library batch、sorting strategy。
    response label 只能作为 evaluation / stratification / leakage audit 变量，不能作为 integration covariate 直接参与校正。

12. 生成 Phase4A / Phase4B handoff。
    明确哪些对象进入 integration，哪些只进 pseudobulk，哪些只做 support，哪些被排除。
    所有 downstream 分析必须读取 Phase3.5 标准对象或 handoff manifest，不能直接读取 raw download 目录。

交付：

* `sc_processing_manifest.frozen_v0.csv`
* `matrix_layer_decision_table.frozen_v0.csv`
* `gene_id_harmonization_report.md`
* `gene_overlap_and_immune_core_coverage.csv`
* `cell_qc_summary_by_sample.csv`
* `cell_filtering_decision_log.csv`
* `sample_cell_barcode_index.csv`
* `gex_tcr_join_qc_report.csv`
* `cell_annotation_source_inventory.csv`
* `hvg_global_pan_cancer.tsv`
* `hvg_by_major_cell_state.tsv`
* `hvg_caution_or_exclusion_list.tsv`
* `batch_covariate_registry.csv`
* `integration_readiness_report.md`
* `analysis_object_registry.frozen_v0.csv`
* `phase3_5_to_phase4_handoff.yaml`

Gate：

* 样本 / 病人 key 无法与 `frozen_v0` 对齐的对象不能进入 Phase4A / Phase4B。
* 表达层类型不清楚的对象不能进入 full integration。
* raw count / normalized layer 混不清的对象只能 support-only。
* QC 后细胞数过低的样本只能 support-only 或 excluded。
* gene overlap 或 immune-core gene coverage 过低的 cohort 不能进入 full pan-cancer integration。
* HVG 被单一 cohort / platform / stress genes 主导时，必须重做 HVG 或降级为 compartment-specific HVG。
* TCR join 不清楚或 join 率过低时，不生成主线 TCR features。
* cell annotation 来源不清的 cohort 不能进入 fine-grained annotation claim，只能进入 coarse-level 或 Phase4A 重新映射。
* response label、outcome、split、survival 字段不得进入 expression feature construction。
* 所有输出必须声明 `input_version=frozen_v0`。

Agent：

* Data processing agent 主执行。
* Immune-state agent 参与 QC 标准和 HVG 设计。
* Controller 审核 `phase3_5_to_phase4_handoff.yaml`。
* 如 QC / annotation / layer 决策争议较大，可启动短期 adversarial data QC reviewer。

[AUTO-RESEARCH 标注]

这里不是生物学模型迭代，而是 processing / QC strategy iteration。
允许小范围比较：

* QC threshold；
* HVG strategy；
* doublet / ambient filtering strategy；
* raw-count vs normalized-object reuse strategy；
* gene ID harmonization strategy。

评价指标：

* retained cell fraction；
* retained sample fraction；
* immune-core gene coverage；
* cross-cohort gene overlap；
* HVG cohort dominance；
* QC covariate balance；
* downstream integration readiness；
* annotation recoverability；
* TCR join quality。

技术探索范围：

主路线：

* 保守 QC；
* response-blind HVG；
* 明确 layer decision；
* 样本级可追溯；
* cohort-level 降级优先于强行合并。

可靠降级：

* full integration → coarse integration；
* coarse integration → reference mapping only；
* reference mapping only → pseudobulk only；
* pseudobulk only → support-only；
* support-only → registry-only。

禁止：

* 直接从 raw download 目录绕过 Phase3.5 进入 Phase4。
* 为了增加样本量强行合并低质量或 layer 不明对象。
* 用 response / outcome 指导 QC、HVG 或 annotation。
* 在没有 barcode / sample / patient 对齐的对象上生成 patient-level feature。
* 把作者已处理对象直接当作可信输入而不记录 assay、layer、QC 来源。
* 为了追求整合漂亮而过度校正 cancer / treatment / biological state。

---

### Phase 4A：Single-cell Integration and Cell-state Harmonization

目的：

在 Phase3.5 通过 QC 和 integration-readiness 的对象上，构建 response-blind 的统一 cell-state annotation 和 integration reference。该阶段回答：

**不同队列的细胞类型和细胞状态能否被放到同一个可比较的生物学坐标系里。**

Phase4A 不生成最终 patient-level modeling matrix，不做 response association，不做 module discovery。它只负责 cell-state integration、annotation harmonization、细胞状态可靠性和后续 Phase4B 的 cell-state ontology。

输入：

* `phase3_5_to_phase4_handoff.yaml`
* `analysis_object_registry.frozen_v0.csv`
* `cell_qc_summary_by_sample.csv`
* `cell_annotation_source_inventory.csv`
* `hvg_global_pan_cancer.tsv`
* `hvg_by_major_cell_state.tsv`
* `hvg_caution_or_exclusion_list.tsv`
* `batch_covariate_registry.csv`
* Phase3 frozen metadata / split / response environment
* integration-ready scRNA objects

任务：

1. 建立 cell-state ontology v6.2。
   至少包含 coarse → mid → selected fine 三层。
   coarse 层覆盖 T/NK、B/plasma、myeloid、DC/APC、CAF/stromal、endothelial、tumor-like、cycling、unknown / low-confidence。
   mid / fine 层覆盖 CD8 cytotoxic、CD8 dysfunctional、CD8 progenitor-like、CD8 terminal exhaustion-like、Treg、NK、inflammatory myeloid、suppressive myeloid、macrophage、C1QC-like macrophage、TREM2-like macrophage、S100A9-like monocyte、DC/APC、TLS/B cell、CAF、endothelial、vascular / angiogenic-like、tumor-intrinsic WNT-exclusion-like 等。
   B2 挂钩必须纳入 TLS/B、neutrophil/TAN/NET、tumor-intrinsic WNT exclusion、CD8 progenitor / terminal exhaustion 的注释与稳定性雷达，但仍保持 response-blind。

2. 执行 marker-based harmonization。
   用 curated marker registry 和 canonical immune / stromal / tumor markers 进行粗到中层映射。
   marker 只能用于 annotation，不得按 response 方向手动挑选。

3. 执行 reference mapping。
   使用合适的 public reference / project reference / existing annotations 进行 label transfer。
   对 annotation conflict 生成 conflict matrix，而不是直接覆盖。

4. 执行 integration model comparison。
   可比较 scVI、scANVI、Harmony、Seurat integration、Scanorama 或其他轻量路线。
   integration 目标是 cell-state comparability，不是消除所有 cohort difference。
   cancer type、treatment context、timepoint 中可能包含真实生物差异，不能无脑校正掉。

5. 评估 integration 质量。
   同时评估 batch mixing、celltype preservation、marker enrichment、known lineage separation、rare-state preservation、cancer/cohort overcorrection risk。
   不允许只用 UMAP 视觉判断。

6. 评估 annotation 可靠性。
   对每个 cell state 输出 reliability score，包括 marker support、reference agreement、cross-cohort reproducibility、sample coverage、cell count、conflict rate。
   低可靠性细分状态必须降级到 coarse ontology。

7. HCC index-context refinement。
   对 HCC 相关 myeloid、macrophage、CAF、endothelial、vascular、tumor-like 状态进行单独可靠性审计。
   如果 myeloid / stromal / endothelial 细分不稳定，后续 HCC context modulation 不允许写细粒度 barrier claim，只能写 broader myeloid-stromal / vascular barrier。

8. pan-cancer annotation consistency audit。
   检查主要 cell-state label 是否跨 cancer / cohort 可复现。
   如果某个 cell state 只存在于单一 cohort 或由单一平台驱动，标记为 cohort-specific / sensitivity-only。

9. 生成 cell-level annotation master。
   每个 cell 记录 `cell_key`、`cohort_id`、`sample_key`、`patient_key`、原始 annotation、harmonized annotation、ontology level、reliability、integration batch、低质量标记、是否进入 Phase4B aggregation。

10. 生成 Phase4B aggregation eligibility。
    明确每个 sample / cell state 是否可用于 fraction、pseudobulk、signature/pathway/TF、TCR feature、HCC refinement、support-only analysis。

交付：

* `cell_state_ontology_v6_2.yaml`
* `cell_state_annotation_master.csv`
* `cell_state_annotation_conflict_matrix.csv`
* `integration_method_comparison_report.md`
* `integration_qc_metrics.csv`
* `cell_state_reliability_report.md`
* `hcc_context_cell_state_refinement_report.md`
* `myeloid_stromal_endothelial_reliability_radar.csv`
* `b2_module_ontology_annotation_radar.csv`
* `phase4a_to_phase4b_aggregation_eligibility.csv`
* `phase4a_decision_manifest.yaml`

Gate：

* 若 integration 破坏 lineage / major cell type separation，则不得采用该 integration 作为主 annotation 坐标。
* 若 batch mixing 好但 marker preservation 差，视为 overcorrection 风险，不得进入主线。
* 若 marker-based 和 reference-mapping 严重冲突，必须降级为 coarse ontology 或进入 adversarial annotation review。
* 若 rare cell state 只在单一 cohort / platform 出现，不得作为 pan-cancer shared state。
* 若 HCC myeloid / CAF / endothelial 细分不可靠，不得写细粒度 HCC barrier claim。
* 若 cell annotation reliability 不足，Phase4B 只能聚合到上一级 ontology。
* response label 不得参与 annotation、integration、cell-state naming。

Agent：

* Immune-state agent 主执行。
* Data processing agent 提供 object / QC 支持。
* Controller 审核是否满足 Phase4B aggregation。
* 对 myeloid / stromal / exhaustion / TLS 注释争议，可启动 adversarial annotation reviewer。

[AUTO-RESEARCH 标注]

可迭代点：

* integration algorithm；
* integration covariates；
* HVG subset；
* ontology resolution；
* marker panel；
* reference atlas；
* myeloid / CAF / endothelial / exhaustion subtyping granularity。

评价指标：

* marker enrichment；
* reference agreement；
* cross-cohort label consistency；
* batch mixing；
* celltype preservation；
* rare-state preservation；
* overcorrection risk；
* spatial mappability；
* downstream Phase4B feature completeness。

技术路线探索范围：

主路线：

* marker-based harmonization + reference mapping 双轨；
* coarse-to-mid ontology 优先；
* fine-grained subtype 必须可靠性过门；
* response-blind annotation。

可靠降级：

* fine subtype → mid-level state；
* mid-level state → coarse lineage；
* integrated annotation → reference mapping only；
* full atlas → per-cancer / per-compartment atlas；
* cell-level fine state → patient-level broad fraction。

并行探索：

* scVI / scANVI latent as auxiliary；
* scFM embedding as optional alignment support；
* HCC-specific sub-atlas；
* myeloid / stromal compartment-specific integration。

禁止：

* 用 response label 命名 cell state。
* 手动挑 marker 迎合预期机制。
* 把 UMAP 好看当成 integration 成功。
* 为了 pan-cancer 一致性抹掉 cancer / organ context biology。
* 在低可靠性细胞亚型上写强机制 claim。
* 未经过 Phase3.5 的对象不得进入 Phase4A。

---

### Phase 4B：Patient-timepoint Immune-state Feature Construction

目的：

把 Phase4A 可靠的 cell-state annotation 和 Phase3.5 标准分析对象聚合为 patient-timepoint 级 immune-state feature matrix，为 Phase5 strong baseline、Phase6 response-blind module discovery、Phase8 SRB Stage-A 和后续 anchor / repair / ZSL 提供统一输入。

该阶段回答：

**每个病人、每个时间点的免疫-肿瘤状态，如何被稳定、可解释、可建模地表示。**

Phase4B 不做最终 responder 发现，不做 dominant barrier 判定，不做 X-class repair。它只构建 feature family，并记录 missingness、confounding、eligibility 和 response-label environment 绑定。

输入：

* `phase4a_decision_manifest.yaml`
* `phase4a_to_phase4b_aggregation_eligibility.csv`
* `cell_state_annotation_master.csv`
* `cell_state_reliability_report.md`
* `analysis_object_registry.frozen_v0.csv`
* `matrix_layer_decision_table.frozen_v0.csv`
* `gex_tcr_join_qc_report.csv`
* Phase3 frozen metadata / split / response environment
* curated gene sets：immune signatures、pathway sets、TF regulons、LR pairs、immune-core gene universe
* optional TCR assets passing Phase3.5 QC

任务：

1. 建立 analysis universe registry。
   不按过细数据角色拆碎主分析，而是建立统一主分析池并附加条件标签。
   至少包含：

   * `primary_pan_cancer_scRNA_universe`
   * `anchor_calibration_subset`
   * `endpoint_sensitivity_subset`
   * `treatment_cleanliness_sensitivity_subset`
   * `hcc_index_context_subset`
   * `support_validation_pool`
   * `excluded_or_registry_only_pool`
     数据角色作为 sample/cohort 标签与敏感性条件保留，不作为默认拆分主分析的理由。

2. 生成 cell-state fraction features。
   按 patient-timepoint / sample 计算 coarse、mid、reliable fine cell states 的 fraction。
   对低细胞数或低可靠性 cell state 进行 shrinkage / censoring / reliability tag。
   输出 raw fraction、normalized fraction、composition-aware transform 和 QC covariates。

3. 生成 cell-state-specific pseudobulk。
   只对通过细胞数、样本数、gene coverage 和 annotation reliability 的 cell state 生成 pseudobulk。
   pseudobulk 必须保留 sample / patient / timepoint / cell-state provenance。
   低细胞数 cell state 只能 broad-bin pseudobulk 或 excluded。

4. 生成 pathway / TF / signature activity。
   在合适层级上计算 pathway、TF regulon、immune signature、exhaustion / cytotoxic / IFN / APC / antigen presentation / myeloid suppression / Treg / CAF / endothelial / TLS / WNT exclusion 等 scores。
   评分方法可以比较 AUCell、GSVA、ssGSEA、decoupler / PROGENy / DoRothEA 或同类方法。
   所有 signature 必须 response-blind 定义。

5. 生成 optional TCR features。
   仅对通过 GEX-TCR join QC 的 cohort 生成 TCR features。
   可包括 clone expansion、clonotype diversity、expanded clone fraction、T cell state × clone coupling、paired pre/post clone persistence / replacement。
   若 paired TCR 不足，则 TCR features 标记为 optional / sensitivity-only。

6. 生成 tumor / stromal / immune QC covariates。
   包括 tumor purity proxy、immune fraction proxy、stromal fraction proxy、cell count、gene count、UMI、mitochondrial fraction、batch、platform、cohort、cancer type、tissue source、timepoint、treatment context 等。
   这些 covariates 后续用于 confounding audit，不应被混作生物机制主特征。

7. 绑定 response-label environment。
   将 response endpoint、RECIST / pathologic response、treatment cleanliness、monotherapy / combo、pre/post、timepoint schema、cohort、cancer type 作为 environment / stratification / sensitivity variables。
   response label 不进入 response-blind module discovery 的 feature construction。

8. 生成 feature dictionary。
   每个 feature 必须有 family、source、level、cell state、gene set、method、input layer、eligible universe、missingness、confounding risk、allowed downstream use。
   允许 downstream use 至少包括：baseline、module discovery、SRB input、anchor calibration、sensitivity-only、support-only、excluded。

9. 生成 feature missingness and confounding audit。
   对每个 feature family 和 analysis universe 计算 missingness、cohort dominance、cancer dominance、platform dependence、cell-count dependence、sample-quality dependence。
   高 missingness 或高 cohort dominance 的 feature family 不进入主线，只能 sensitivity-only。

10. 生成 patient-timepoint immune-state feature matrix。
    输出统一矩阵，行单位为 patient-timepoint 或 sample-timepoint，列为通过 eligibility 的 feature。
    同时输出主矩阵、feature-family 矩阵和 sensitivity 矩阵。
    所有矩阵必须可回溯到 Phase3 frozen input、Phase3.5 processing 和 Phase4A annotation。

11. 生成 Phase5 handoff。
    明确哪些 feature family 和 universe 可进入 strong baseline / confounding audit，哪些只能 sensitivity，哪些被 blocked。
    后续 Phase5 不得直接绕过 Phase4B 读取 cell-level object 构造新 feature。

交付：

* `analysis_universe_registry_v6_2.csv`
* `cell_state_fraction_matrix.csv`
* `cell_state_specific_pseudobulk_registry.csv`
* `cell_state_specific_pseudobulk_matrix.parquet`
* `pathway_tf_signature_activity_matrix.csv`
* `optional_tcr_feature_matrix.csv`
* `sample_patient_timepoint_qc_covariates.csv`
* `immune_state_feature_matrix.frozen_v0.csv`
* `immune_state_feature_matrix_by_family/`
* `feature_dictionary_v6_2.csv`
* `feature_missingness_report.csv`
* `feature_confounding_tags.csv`
* `response_environment_binding_table.csv`
* `phase4b_to_phase5_handoff.yaml`
* `PHASE4_IMMUNE_STATE_FEATURE_CONSTRUCTION_REPORT.md`

Gate：

* 主矩阵必须使用 patient-level split，不允许同一 patient 跨 train / val / test。
* feature column 中不得包含 response、outcome、survival、split、treatment response 等泄漏字段。
* 高 missingness feature family 不进入 primary universe。
* 强 cohort / cancer / platform dominance 的 feature 必须标记 confounding risk。
* 低可靠性 cell state 不生成 fine-level patient feature，只能生成 coarse feature。
* 低细胞数 pseudobulk 不进入主线。
* TCR join 不过门的 cohort 不生成 TCR 主特征。
* pathologic response 与 RECIST 不得合并为同质监督标签，只能通过 endpoint sensitivity 处理。
* support-only / bulk / spatial / combo 数据不得替代 clean scRNA anchor。
* 所有 supervised downstream model 只能通过 Phase5 使用 patient-level split 和 Phase4B handoff。

Agent：

* Immune-state agent 主执行。
* Data processing agent 维护 matrix/provenance。
* Controller 审核 feature dictionary、universe registry 和 handoff。
* Phase5 baseline agent 只允许读取 `phase4b_to_phase5_handoff.yaml` 中允许的矩阵与 universe。

[AUTO-RESEARCH 标注]

可迭代点：

* aggregation level；
* pseudobulk minimum cell threshold；
* signature/pathway/TF scoring method；
* feature family inclusion rule；
* reliability-weighted feature construction；
* compositional transform；
* TCR feature definition；
* missingness / confounding threshold。

评价指标：

* feature completeness；
* missingness；
* cohort dominance；
* cancer dominance；
* cell-count dependence；
* marker/pathway biological sanity；
* cross-cohort reproducibility；
* spatial mappability；
* downstream baseline utility；
* negative-control leakage risk。

技术路线探索范围：

主路线：

* patient-timepoint aggregation；
* response-blind feature construction；
* feature dictionary + confounding tags；
* unified primary analysis pool with conditional sensitivity labels；
* strict handoff to Phase5。

可靠降级：

* fine cell-state features → mid/coarse features；
* cell-state-specific pseudobulk → lineage-level pseudobulk；
* pathway/TF/signature family → curated signature only；
* TCR main features → TCR sensitivity-only；
* primary universe → support/sensitivity universe。

并行探索：

* reliability-weighted aggregation；
* compartment-specific feature matrices；
* HCC index-context feature layer；
* anchor calibration subset feature layer；
* endpoint sensitivity feature layer。

禁止：

* 用 response label 构造 feature。
* 为了提高下游 AUC 临时增加 response-derived feature。
* 把 support-only 数据混入 clean anchor supervised universe。
* 在 feature construction 阶段做模块发现或机制 claim。
* 让 Phase5 / Phase6 绕过 Phase4B 直接从 cell-level object 重新造 feature。
* 把数据角色拆成互相孤立的小分析，破坏 pan-cancer 主线。

---

## 5. Stage 2：Core Build

### Phase 5：Strong Baseline and Confounding Audit

目的：

建立所有复杂模型的压力 benchmark。它不是为了最高分，而是判断信号是否已经被简单 feature、cell fraction 或 cohort/cancer confounding 解释。

任务：

1. cell fraction baseline；
2. curated signature baseline；
3. pathway / TF baseline；
4. pseudobulk module baseline；
5. elastic net / logistic；
6. limited XGBoost / random forest；
7. scVI/patient latent baseline；
8. negative controls；
9. label leakage audit；
10. cancer/cohort/batch/cell fraction confounding audit。

交付：

* `baseline_results.csv`
* `baseline_model_manifest.yaml`
* `confounding_audit_report.md`
* `negative_control_report.md`
* `feature_family_downstream_eligibility.csv`

Gate：

* 如果 baseline 显示信号主要来自 cohort/cancer/batch，相关 feature 降级；
* 如果 cell fraction baseline 已解释主要结果，后续机制模型必须证明增量不是 cell fraction surrogate；
* 如果 response leakage，停止监督分析并修复；
* 如果简单 baseline 已达到极强但机制不清，SRB 仍可做解释模型，但不能只靠预测性能宣传。

Agent：

* Baseline agent；
* 独立 hostile statistics reviewer；
* Controller 决定哪些 feature family 进入 module / SRB。

[AUTO-RESEARCH 标注]

可迭代点：

* regularization；
* feature family；
* split strategy；
* calibration；
* missingness handling；
* class imbalance handling。

评价指标：

* LOCO performance；
* leave-center-out performance；
* calibration；
* negative control collapse；
* cell-fraction sensitivity；
* response-environment robustness。

技术路线探索范围：

主路线：

* interpretable baseline + limited strong ML。

可靠降级：

* signature/pathway/fraction-only baseline；
* cohort-specific descriptive statistics。

并行探索：

* scVI latent baseline；
* simple MLP patient classifier as stress test；
* domain-adversarial baseline if cohort leakage severe。

禁止：

* baseline 变成模型动物园；
* 只看 AUC；
* 在 patient split 外混用同患者样本；
* 用 baseline 高分直接做机制 claim。

---

### Phase 6：Response-blind Module Discovery

目的：

构建 SRB 的 module ontology 和 latent 分块单位。module 是后续 barrier gate、SRB、SNGM、ZSL、evidence card 的共同语义单位。

任务：

1. response-blind module discovery；
2. module membership；
3. module scoring；
4. carrier cell-state annotation；
5. pathway / TF / LR annotation；
6. module stability；
7. spatial mappability；
8. perturbation mappability；
9. module version freeze。

交付：

* `module_membership.v0.csv`
* `module_score_matrix.v0.csv`
* `module_attribute_table.v0.csv`
* `module_stability_report.md`
* `module_mappability_report.md`
* `module_freeze_manifest.yaml`

Gate：

* 模块必须 response-blind 产生；
* 模块必须可在 patient level 评分；
* 不稳定高分辨率模块合并；
* 不能用 response-derived module 进入 SRB；
* 模块数量过多会导致 SRB 不稳，必须压缩到可解释范围。

Agent：

* Module agent；
* Parallel method comparison agent 可单独执行 WGCNA/NMF/sparse AE 对比；
* Controller 决定 module freeze。

[AUTO-RESEARCH 标注]

可迭代点：

* WGCNA-like resolution；
* NMF rank；
* graph community resolution；
* sparse autoencoder bottleneck；
* gene universe；
* cell-state-specific vs pan-cell module。

评价指标：

* bootstrap membership stability；
* cross-cohort reproducibility；
* pathway coherence；
* carrier cell-state clarity；
* spatial mappability；
* perturbation mappability；
* downstream SRB reconstruction utility。

技术路线探索范围：

主路线：

* WGCNA/NMF/graph community consensus；
* sparse dictionary refinement as upgrade。

可靠降级：

* broader immune programs；
* curated immune programs；
* pathway-level modules。

并行探索：

* sparse autoencoder modules；
* topic modeling；
* cell-state-specific module discovery。

禁止：

* response-supervised module discovery；
* 为得到漂亮机制强行拆小模块；
* 不稳定模块进入 dominant barrier。

---

### Phase 7：Barrier Identifiability Gate

目的：

在 dominant barrier 和 repair transport 前，判断 modules 是否可分。这里是 v6.2 防止“归因 artifact”的关键 gate。

任务：

1. bootstrap stability；
2. module assignment flip rate；
3. collinearity / condition number；
4. multi-resolution consistency；
5. deconvolution consistency；
6. carrier cell-state consistency；
7. perturbation coverage；
8. spatial mappability；
9. coarse barrier merge。

交付：

* `barrier_identifiability_table.csv`
* `dominant_barrier_eligibility.csv`
* `coarse_barrier_merge_map.csv`
* `unstable_barrier_log.md`

Gate：

* 未通过 gate 的 module 不能进入 dominant barrier；
* 共线模块必须合并，不硬拆；
* HCC-specific residual claim 必须建立在可辨识 barrier 上；
* barrier gate 失败不等于项目失败，但会限制 repair 和 X-class claim。

Agent：

* Barrier gate agent；
* Hostile reviewer 重点审查；
* Controller 冻结 eligible barrier set。

[AUTO-RESEARCH 标注]

可迭代点：

* module resolution；
* bootstrap parameters；
* deconvolution strategy；
* collinearity threshold；
* merge rule；
* carrier-cell consistency threshold。

评价指标：

* flip rate；
* condition number；
* bootstrap support；
* cross-resolution agreement；
* downstream repair stability；
* spatial/perturbation mappability。

技术路线探索范围：

主路线：

* deterministic gate + bootstrap；
* explicit coarse merge。

可靠降级：

* broad barrier class；
* no dominant barrier claim，仅做 module-level evidence。

并行探索：

* Bayesian latent class / mixture model for barrier assignment；
* graph-based barrier clustering。

禁止：

* 对不可分模块做精细机制 claim；
* 用 NN attribution 替代可辨识性审计；
* 把 myeloid/Treg/stromal 共线硬说成三个独立 barrier。

---

### Phase 8：SRB Stage-A Representation

目的：

训练 Structured Repair Backbone 的稳定 module-factorized latent。此阶段不急于 repair 或 ZSL，先让 latent 可靠。

任务：

1. module-factorized encoder；
2. reconstruction head；
3. graph attention regularization；
4. environment-aware representation；
5. deconfounding / adversarial audit；
6. Stage-A latent freeze；
7. comparison against non-neural embeddings。

交付：

* `srb_stageA_model_manifest.yaml`
* `module_factorized_latent.csv`
* `srb_reconstruction_report.md`
* `srb_leakage_audit.md`
* `srb_stageA_eval_record.yaml`

Gate：

* latent 必须可重建 module features；
* latent 不能主要编码 cohort/cancer/batch；
* module latent 必须有解释性；
* 若 SRB Stage-A 不超过 PCA/NMF/scVI 在关键稳定性指标上的表现，降级为 static module matrix + simpler heads。

Agent：

* SRB agent；
* Baseline comparison agent；
* Controller 审查是否冻结 Stage-A latent。

[AUTO-RESEARCH 标注]

可迭代点：

* latent dimension；
* module block size；
* shared vs module-specific encoder；
* graph prior；
* reconstruction loss；
* adversarial weight；
* dropout / regularization；
* pretraining objective。

评价指标：

* reconstruction quality；
* cohort leakage；
* cancer leakage；
* response label leakage；
* module interpretability；
* downstream context/anchor utility；
* calibration and stability。

技术路线探索范围：

主路线：

* module-factorized MLP encoder；
* reconstruction + graph regularization；
* environment-aware audit。

可靠降级：

* module score matrix；
* PCA/NMF embedding；
* scVI latent + module projection。

并行探索：

* attention pooling over cell states；
* variational module encoder；
* contrastive patient-timepoint encoder。

禁止：

* 黑箱 global embedding 取代 module-factorized latent；
* 没有 leakage audit 就进入 supervised heads；
* 为了提高 response AUC 牺牲模块解释性。

---

### Phase 9：Context Modulation and Anti-PD1 Anchor

目的：

在 SRB latent 上识别 shared anti-PD1 failure modules、context-modulated modules，并建立 anti-PD1 responder-compatible direction。

任务：

1. context head；
2. anchor head；
3. clean / expanded / support anchor 分层；
4. response-label environment 处理；
5. HCC index-context residual effect；（B1 挂钩：etiology 作低权重可插拔 sensitivity 协变量，覆盖足时启用；B4 挂钩：数据允许时区分 primary/adaptive/acquired resistance；B5 挂钩：有 paired TCR 时纳入克隆替换 vs 再激活）
6. anchor reliability；
7. R1 / R1-lite / R2 routing。

交付：

* `shared_antipd1_modules.csv`
* `context_modulated_modules.csv`
* `hcc_index_context_residuals.csv`
* `anchor_reliability_report.md`
* `responder_compatible_direction.csv`
* `route_verdict.yaml`

Gate：

* clean anchor 强：R1；
* anchor partial：R1-lite；
* anchor failed but HCC/context strong：R2；
* anchor failed and context weak：只能做 Result Layer A / resource-like story；
* response label environment 解释主要效应时，anchor claim 降级。

Agent：

* SRB agent；
* Anchor rescue branch agent；
* Controller 做 route verdict；
* 用户确认主路线切换。

[AUTO-RESEARCH 标注]

可迭代点：

* clean/expanded/support anchor weighting；
* response environment treatment；
* context modulation penalty；
* anchor head architecture；
* label harmonization；
* HCC residual interaction modeling。

评价指标：

* cross-cohort direction consistency；
* leave-one-cancer stability；
* response-label environment robustness；
* external bulk direction；
* TCR support if available；
* confounding audit。

技术路线探索范围：

主路线：

* SRB context + anchor heads；
* mixed-effect / meta-analysis cross-check。

可靠降级：

* responder centroid direction；
* module-level meta-analysis；
* R1-lite / R2 route。

并行探索：

* invariant risk minimization style penalty；
* domain-adversarial context head；
* TCR-informed anchor if data supports。

禁止：

* support anchor 替代 clean anchor；
* 高混杂 PD1+X 队列反推 clean PD1 anchor；
* response label source 当作普通预测 feature。

---

### Phase 10：CRTN / Transport Head

目的：

估计 dominant barrier repair 是否将 anti-PD1 failure state 推向 responder-compatible direction。

任务：

1. dominant barrier selection；
2. transport head training；
3. CRTN L1 validation；
4. optional L2 held-out combination validation；
5. compare against OT / linear / permutation / signature reversal；
6. perturbation coverage-gated consistency；
7. CRTN credibility level assignment。

交付：

* `crtn_model_manifest.yaml`
* `repair_displacement_matrix.csv`
* `barrier_repair_effects.csv`
* `xclass_conditioned_repair_effects.csv`
* `crtn_baseline_comparison.csv`
* `crtn_credibility_report.md`

Gate：

* L1 不达标：CRTN 降级为 descriptive / OT repair；
* L2 只有真实 PD1+X 或 pre-specified held-out combination 支持时才可写；
* 无湿实验不写 proven causal therapeutic effect；
* repair direction 与 spatial / perturbation / external 冲突时降级。

Agent：

* SRB transport agent；
* Perturbation agent 提供 coverage；
* Controller 整合 CRTN credibility；
* Hostile reviewer 检查是否 overclaim。

[AUTO-RESEARCH 标注]

可迭代点：

* single-step residual MLP；
* conditional flow matching；
* entropic OT；
* loss weights；
* graph locality；
* perturbation consistency；
* condition descriptors；
* uncertainty estimation。

评价指标：

* held-out displacement prediction；
* responder-direction movement；
* permutation null；
* negative control repair effect near zero；
* graph locality；
* perturbation consistency；
* spatial consistency；
* external direction consistency。

技术路线探索范围：

主路线：

* residual MLP transport；
* CRTN L1 validation。

可靠降级：

* entropic OT + barycentric projection；
* linear displacement；
* nearest-responder displacement。

并行探索：

* conditional flow matching；
* neural ODE / flow spike；
* causal graph-constrained transport。

禁止：

* 把 model-based counterfactual 写成真实疗效；
* 未过 barrier gate 的 module 进入 transport；
* 不做 baseline 比较就宣传 CRTN。

---

### Phase 11：Spatial Asset Module and SNGM

目的：

把空间数据从补图升级为组织生态裁判，但不强行并入 SRB 主干。

任务：

1. 读取 Phase 2.5 标准化输出，不直接读取 raw download 目录；
2. spatial asset audit；
3. module spatial projection；
4. cell-state deconvolution / segmentation；
5. statistical niche baseline；
6. myeloid-tumor niche；
7. T cell exclusion；
8. APC-T cell niche；
9. CAF/vascular barrier；
10. LR spatial proximity；
11. SNGM if data supports；
12. spatial support/conflict verdict。

交付：

* `spatial_asset_registry.csv`
* `spatial_module_projection.csv`
* `niche_score_table.csv`
* `spatial_permutation_report.md`
* `sngm_model_manifest.yaml`, if used
* `spatial_support_conflict_verdict.csv`

Gate：

* 未通过 Phase 2.5 的样本不得进入 SNGM；
* spatial metadata 不足：只做 localization support；
* 空间质量不足：不训练 GNN；
* spatial 与 repair 冲突：candidate 降级；
* spatial 支持强但其他证据弱：保留为 niche hypothesis。

Agent：

* Spatial asset agent；
* SNGM agent 可独立；
* Controller 不允许 SNGM 结果自动改主线。

[AUTO-RESEARCH 标注]

可迭代点：

* neighborhood radius / k；
* niche definition；
* permutation null；
* deconvolution method；
* LR filtering；
* graph architecture；
* contrastive alignment；
* region label integration。

评价指标：

* niche reproducibility；
* permutation FDR；
* tissue-region consistency；
* module-niche coherence；
* response/failure association where valid；
* consistency with SRB dominant barrier；
* histology interpretability。

技术路线探索范围：

主路线：

* statistical niche baseline；
* SNGM only if spatial data supports。

可靠降级：

* deconvolution-only niche；
* module localization；
* region enrichment。

并行探索：

* hetero-GNN；
* scRNA-spatial contrastive alignment；
* GraphTME/ISCHIA-style niche construction；
* mIF/IHC marker panel design.

禁止：

* SNGM 并入 SRB 端到端主干；
* spatial proximity 单独证明功能互作；
* metadata 弱的空间数据支持 response claim。

---

### Phase 12：Perturbation Coverage and X-class Descriptor Ontology

目的：

构建 X-class 机制描述空间，并判断每个 X-class 是否有足够 perturbation / target / pathway support。

任务：

1. X-class ontology；（B3 挂钩：纳入当代组合 anti-VEGF/TKI/anti-CTLA4/anti-TIGIT/anti-LAG3/anti-CSF1R/anti-TGFβ/adenosine；并把临床失败组合 LEAP-002、COSMIC-312 作 ZSL/CRTN held-out 负样本，anti-TIGIT 作正向用例。B6 挂钩：显式声明 perturbation 先验对 myeloid/stromal/TLS/中性粒 barrier 的覆盖偏倚——缺支持≠证据为负）
2. mechanism descriptor；
3. perturbation signature mapping；
4. target family mapping；
5. pathway / TF / drug family mapping；
6. perturbation coverage gate；
7. direction support/conflict；
8. X-class eligibility for alignment head。

交付：

* `xclass_descriptor_ontology.yaml`
* `xclass_attribute_table.csv`
* `perturbation_coverage_matrix.csv`
* `target_module_mapping.csv`
* `perturbation_direction_support.csv`
* `xclass_eligibility_table.csv`

Gate：

* descriptor 不得包含 response/outcome leakage；
* coverage 低时不允许强 perturbation claim；
* perturbation 与 CRTN 方向冲突时降级；
* X-class descriptor 必须是 soft attribute bundle，不是硬标签。

Agent：

* Perturbation/ZSL agent；
* Domain expert reviewer 可短期介入；
* Controller 审核 descriptor leakage。

[AUTO-RESEARCH 标注]

可迭代点：

* X-class 粒度；
* marker set；
* pathway descriptor；
* target family mapping；
* perturbation aggregation；
* coverage weighting；
* conflict scoring。

评价指标：

* coverage；
* direction AUROC where labels exist；
* module-X compatibility；
* expert plausibility；
* semantic leakage audit；
* downstream ZSL calibration。

技术路线探索范围：

主路线：

* soft mechanism class descriptors；
* rule-based perturbation coverage gate。

可靠降级：

* coarse X-class；
* curated pathway/target mapping；
* evidence-card-only X-class scoring。

并行探索：

* text/ontology embedding；
* perturbation encoder；
* drug family prototype embedding。

禁止：

* descriptor 中写入 response 信息；
* 药名 taxonomy 取代机制 taxonomy；
* LINCS reversal 单独决定 X-class 提名。

---

### Phase 13：ZSL Alignment Head

目的：

检验 learned mechanism space 是否支持未见 X-class / 未见 PD1+X 组合的组合外推。ZSL 不是锦上添花，它是 v6.2 的正式外推层，但不替代机制证据。

任务：

1. Mechanism Prototype ZSL；
2. Combination ZSL；
3. Graph/Niche-regularized Mechanism ZSL；
4. Cancer ZSL validation only；
5. leave-one-X-class；
6. leave-one-combination；
7. open-world abstention；
8. calibration；
9. simple zero-shot baselines；
10. semantic leakage audit。

交付：

* `zsl_split_manifest.yaml`
* `mechanism_zsl_results.csv`
* `combination_zsl_results.csv`
* `zsl_baseline_comparison.csv`
* `zsl_calibration_abstention_report.md`
* `zsl_leakage_audit.md`
* `zsl_promotion_verdict.yaml`

Gate：

* 未超过 zero-shot-capable baseline：降级为 auxiliary；
* calibration 失败：不进主文 claim；
* abstention 失败：不写 open-world generalization；
* semantic leakage：ZSL 结果无效；
* Cancer ZSL 默认只是 validation，不是主 claim。

Agent：

* ZSL agent；
* Independent leakage reviewer；
* Controller 决定 ZSL 是否进入 Result Layer C。

[AUTO-RESEARCH 标注]

可迭代点：

* prototype descriptor；
* contrastive loss；
* pairwise ranking loss；
* graph/niche regularization；
* abstention threshold；
* calibration method；
* descriptor embedding；
* additive vs dual-encoder vs tri-encoder。

评价指标：

* leave-one-X-class Top-k；
* leave-one-combination pairwise ranking；
* calibration；
* abstention AUROC；
* evidence consistency；
* conflict rate；
* external direction consistency。

技术路线探索范围：

主路线：

* SRB alignment head；
* mechanism prototype ZSL；
* combination ZSL。

可靠降级：

* nearest prototype；
* attribute matching；
* additive PD1+X baseline；
* evidence-card-only ranking。

并行探索：

* graph/niche-regularized ZSL；
* text/ontology embedding；
* contrastive mechanism alignment；
* retrieval-augmented mechanism prototype.

禁止：

* 把 ZSL 写成发现全新生物学；
* 隐藏简单 baseline；
* closed-world accuracy 伪装 zero-shot；
* Cancer ZSL 失败否定主框架。

---

### Phase 14：Evidence Card Integration

目的：

将 module、barrier、SRB、CRTN、SNGM、perturbation、ZSL、external validation 统一到 candidate evidence cards 中。

任务：

1. candidate module evidence；
2. candidate X-class evidence；
3. evidence tier assignment；
4. conflict-first rule；
5. claim boundary；
6. candidate ranking；
7. validation priority；
8. result layer classification。

交付：

* `candidate_evidence_cards.csv`
* `evidence_items.csv`
* `conflict_log.md`
* `claim_tier_table.csv`
* `candidate_ranking.csv`
* `result_layer_assignment.yaml`

Gate：

* evidence conflict 优先于 ranking score；
* barrier 未过 identifiability gate，不能进强 claim；
* repair 与 spatial / perturb / external 冲突，降级；
* 单一证据不能升 L3；
* Stage 2 若无 evidence card v0，不进 Stage 3。

Agent：

* Integration agent；
* Controller 审查；
* Hostile reviewer audit；
* 用户确认 Top candidate / story lock。

[AUTO-RESEARCH 标注]

可迭代点：

* evidence weights；
* conflict penalties；
* calibration；
* Bayesian fusion vs rule tier；
* candidate score aggregation；
* negative controls；
* reliability curves。

评价指标：

* conflict detection；
* tier stability；
* bootstrap consistency；
* external validation agreement；
* reviewer interpretability；
* candidate prioritization robustness。

技术路线探索范围：

主路线：

* rule-based evidence tier + conflict-first；
* ranking 和 claim tier 分离。

可靠降级：

* manual evidence card；
* L1/L2/L3 simplified tier。

并行探索：

* Bayesian evidence fusion；
* learned calibration model if enough historical candidates；
* graph evidence propagation.

禁止：

* 总分抵消关键冲突；
* ZSL 高分自动升 claim；
* 把 evidence card 写成装饰表。

---

## 6. Stage 3：Manuscript Package

### Phase 15：External Validation and Tissue Validation Priority

目的：

为 Result Layer A/B/C 提供外部方向支持和组织/功能验证路径。

任务：

1. external bulk module projection；
2. HCC index-context external validation；
3. treatment-context sensitivity；
4. response-label sensitivity；
5. failure-case reporting；
6. tissue marker selection；
7. IHC/mIF/spatial validation candidate design；
8. optional functional validation shortlist。

交付：

* `external_validation_results.csv`
* `bulk_projection_report.md`
* `hcc_context_validation_report.md`
* `tissue_validation_candidate_panel.csv`
* `wetlab_priority_table.csv`
* `failure_case_report.md`

Gate：

* external direction 反向：candidate 降级；
* tissue marker 不可验证：L3 降级；
* HCC 只作为 index context，不写代表泛癌；
* 没有实验验证时不写 proven therapeutic effect。

Agent：

* External validation agent；
* Tissue validation planning agent；
* Controller 审核 evidence card 更新。

[AUTO-RESEARCH 标注]

可迭代点：

* bulk projection method；
* purity correction；
* treatment-context stratification；
* marker panel selection；
* spatial/tissue validation design。

评价指标：

* direction consistency；
* survival/response association where valid；
* purity sensitivity；
* marker specificity；
* tissue feasibility；
* candidate rank robustness。

技术路线探索范围：

主路线：

* module-level projection；
* external direction support；
* tissue marker feasibility。

可靠降级：

* external support only；
* failure-case reporting；
* no L3 claim.

并行探索：

* multiple bulk projection methods；
* deconvolution-adjusted projection；
* public mIF/spatial tissue reanalysis。

禁止：

* bulk projection 当机制证明；
* survival trend 当 anti-PD1 修复证明；
* 无组织验证路径却写强 translational claim。

---

### Phase 16：Result Layer Assembly

目的：

决定论文最终可落在哪个结果层级，而不是无限补分析。

Result Layer A：anti-PD1 failure module atlas and context modulation

需要：

* shared anti-PD1 failure modules；
* context-modulated modules；
* HCC index-context residuals；
* barrier identifiability map；
* anchor reliability；
* baseline / robustness validation。

Result Layer B：counterfactual barrier repair and spatial niche support

需要：

* SRB transport；
* dominant barrier repair direction；
* SNGM or statistical spatial niche support；
* perturbation coverage gate；
* CRTN credibility level；
* repair-compatible mechanism hypotheses。

Result Layer C：X-class zero-shot extrapolation and PD1+X hypothesis nomination

需要：

* X-class descriptor ontology；
* alignment head / compositional ZSL；
* unseen X-class / combination ranking；
* evidence cards；
* experimentally testable PD1+X hypotheses。

交付：

* `result_layer_decision.yaml`
* `main_storyline.md`
* `figure_plan.md`
* `claim_boundary_final.md`

Gate：

* 只能选择证据最强、冲突最少的层级作为主故事；
* Layer C 不稳时，不强行包装；
* Layer A 如果很强，也可成为保底 manuscript direction；
* 主故事必须由 evidence card 支撑。

Agent：

* Manuscript integration agent；
* Controller + hostile reviewer；
* 用户确认主故事。

技术探索范围：

无新模型。这里是科学判断，不是继续堆算法。

---

### Phase 17：Figure Source and Manuscript Package

目的：

把证据链转化成主图、补图、表格、方法和复现包。

任务：

1. 主图源表；
2. 补图源表；
3. evidence card figure；
4. SRB model schematic；
5. SNGM spatial figure；
6. CRTN repair trajectory figure；
7. ZSL generalization figure；
8. external validation figure；
9. failure/downgrade figure；
10. reproducibility package。

交付：

* `main_figure_source_tables/`
* `supplementary_figure_source_tables/`
* `figure_package_plan.md`
* `methods_skeleton.md`
* `result_summary_tables/`
* `minimal_reproducibility_package/`
* `manuscript_outline.md`

Gate：

* 每张图必须可追溯到 frozen input 和 evidence card；
* 不允许图先行、证据后补；
* 方法图不能掩盖降级；
* 所有模型必须有 baseline/ablation 位置。

Agent：

* Manuscript agent；
* Figure package agent；
* Reproducibility agent；
* Controller final audit。

技术探索范围：

只允许可视化和复现工程优化，不再改变核心模型。

---

### Phase 18：Hostile Reviewer Audit and Freeze

目的：

模拟敌对审稿人，从数据、统计、机制、空间、ZSL、临床边界和复现角度打穿项目。

审计问题：

1. 是否存在 cohort/cancer/batch confounding；
2. 是否存在 response-label leakage；
3. 是否把 cell fraction surrogate 当机制；
4. barrier 是否不可辨识；
5. graph attention 是否被过度写成因果；
6. CRTN 是否只是 signature reversal；
7. SNGM 是否只是空间展示；
8. perturbation prior 是否被当真值；
9. ZSL 是否存在 semantic leakage；
10. external validation 是否方向一致；
11. HCC 是否被错误写成泛癌代表；
12. 是否有 direct clinical treatment claim；
13. 是否所有 figure source 可复现。

交付：

* `hostile_reviewer_audit.md`
* `claim_overreach_checklist.md`
* `go_no_go_decision.yaml`
* `revision_required_register.csv`
* `final_freeze_manifest.yaml`

Gate：

* hard blocker 未解决，不 freeze；
* claim overreach 必须降级；
* 关键冲突必须进入 manuscript limitation 或 candidate downgrade；
* freeze 后只允许勘误，不允许重写核心模型。

Agent：

* Independent hostile reviewer agent；
* Controller；
* 用户最终确认。

技术探索范围：

不允许新增主模型。只允许修复证据、降级 claim、补充审计。

---

## 7. AutoResearch-style 迭代总表

以下节点允许并建议使用 AutoResearch-style 迭代。每个节点必须以 v0 baseline 开始，再进行 v1/v2 探索。

| 节点                        | 可迭代内容                                               | 主指标                                              | 降级标准                        | 是否允许并行 agent    |
| ------------------------- | --------------------------------------------------- | ------------------------------------------------ | --------------------------- | --------------- |
| Immune-state construction | annotation resolution, mapping, pathway scoring     | label consistency, module stability, mappability | coarse ontology             | 是               |
| Strong baseline           | feature family, regularization, calibration         | LOCO, calibration, leakage                       | interpretable baseline only | 是               |
| Module discovery          | WGCNA/NMF/sparse AE, resolution                     | stability, enrichment, mappability               | broad immune programs       | 是               |
| Barrier gate              | merge rule, collinearity threshold                  | flip rate, condition number, repair stability    | coarse barrier              | 是               |
| SRB Stage-A               | latent dim, graph prior, loss weights               | reconstruction, leakage, downstream utility      | module matrix/PCA           | 是               |
| Context/Anchor            | environment penalty, anchor weighting               | direction stability, external support            | R1-lite/R2                  | 是               |
| CRTN                      | MLP, OT, flow matching, graph locality              | displacement prediction, null, consistency       | OT/linear repair            | 是               |
| SNGM                      | niche radius, graph architecture, null model        | FDR, reproducibility, tissue coherence           | statistical niche           | 是               |
| Perturbation              | coverage aggregation, target mapping                | direction support, conflict rate                 | rule-based support          | 是               |
| ZSL                       | prototype, dual encoder, graph/niche regularization | leave-one-X, calibration, abstention             | prototype/attribute only    | 是               |
| Evidence card             | weights, conflict penalty, calibration              | tier stability, external agreement               | rule tier                   | Controller only |

AutoResearch-style 的边界：

* 每轮只改一个主要因素；
* 每轮必须写入 iteration log；
* 每轮必须与 v0 和当前 best 比较；
* 不能因为单指标提升而接受产生更大混杂或更强冲突的方案；
* 探索结果默认是 proposal，除非通过 gate，否则不进主线。

---

## 8. Agent / Session 切换规则

### 8.1 何时开独立 agent

适合独立 agent：

* data asset audit；
* immune-state construction；
* module discovery comparison；
* baseline/confounding audit；
* SRB model training；
* spatial/SNGM；
* perturbation/ZSL；
* external validation；
* hostile reviewer audit。

不适合独立 agent：

* source-of-truth 裁决；
* claim tier 升级；
  -主故事锁定；
* Result Layer 选择；
* final freeze。

这些必须由 Controller 汇总后交给用户确认。

### 8.2 每个 branch agent 的输出要求

每个 branch agent 至少输出：

* `handoff_manifest.yaml`
* `decision_log.md`
* `evidence_items.csv`
* `conflict_or_failure_log.md`
* `next_stage_readiness.yaml`

如果涉及模型迭代，还必须输出：

* `eval_record.yaml`
* `iteration_log_entry.md`
* `baseline_comparison.csv`
* `adopt_reject_recommendation.md`

### 8.3 Integration rehearsal

必须在以下节点进行 integration rehearsal：

1. Stage 1 结束；
2. Module + barrier gate 结束；
3. SRB Stage-A 结束；
4. CRTN + SNGM + perturbation 初步完成；
5. ZSL 初步完成；
6. Evidence card v0 完成；
7. Manuscript package freeze 前。

Integration rehearsal 的核心问题：

* 各 branch 是否使用同一 frozen input；
* 证据是否能进入 evidence card；
* 有无方向冲突；
* 哪些 branch 需要降级；
* 是否能进入下一阶段。

---

## 9. Go / No-Go 总规则

### Stage 1 Go

需要：

* frozen input 可用；
* data coverage 足够支撑至少 Result Layer A；
* clean/expanded/support anchor 分层清楚；
* spatial/perturbation/external coverage 已审计；
* eval harness 可用；
* MVP-derived 结果已隔离。

否则：

* 修复数据；
* 降级目标；
* 不进入 SRB Core Build。

### Stage 2 Go

需要：

* L1 robust modules；
* barrier identifiability gate 至少部分通过；
* anchor strong 或 partial，或 R2 route 明确；
* SRB Stage-A latent 可用；
* CRTN 或降级 repair direction 可解释；
* 至少一个 evidence card v0 无关键冲突。

否则：

* 不进入 manuscript package；
* 回到数据、module、barrier 或 anchor 问题复盘。

### Stage 3 Go

需要：

* 至少一个 Result Layer 成立；
* 主故事和 claim tier 匹配；
* hostile reviewer audit 无 hard blocker；
* figure source 可追踪；
* minimal reproducibility package 可运行；
* 所有 overclaim 已降级。

否则：

* 不 freeze manuscript；
* 降级 claim 或回补证据。

---

## 10. 最终施工顺序压缩版

Phase 0.1 Project Control Brief
Phase 0.2 Evaluation Harness
Phase 1 Source-of-truth and Inheritance Lock
Phase 2 Data Go / No-Go Audit
Phase 2.5 Spatial Data Standardization Gate
Phase 3 Frozen Input Family
Phase 3.5 Single-cell Processing QC Gate
Phase 4A Single-cell Integration and Cell-state Harmonization
Phase 4B Patient-timepoint Immune-state Feature Construction
Phase 5 Strong Baseline and Confounding Audit
Phase 6 Response-blind Module Discovery
Phase 7 Barrier Identifiability Gate
Phase 8 SRB Stage-A Representation
Phase 9 Context Modulation and Anti-PD1 Anchor
Phase 10 CRTN / Transport Head
Phase 11 Spatial Asset Module and SNGM
Phase 12 Perturbation Coverage and X-class Descriptor Ontology
Phase 13 ZSL Alignment Head
Phase 14 Evidence Card Integration
Phase 15 External Validation and Tissue Validation Priority
Phase 16 Result Layer Assembly
Phase 17 Figure Source and Manuscript Package
Phase 18 Hostile Reviewer Audit and Freeze

核心路线可以压缩成一句：

**冻结输入 → 标准化空间资产 → 构建 immune-state → 建强基线 → 发现模块 → 过 barrier gate → 训练 SRB → 建 anchor/context → 做 repair transport → 做 spatial niche → 做 perturbation/ZSL → 合成 evidence card → 选择 Result Layer → 生成 manuscript package。**

---

## 11. Roadmap 的关键判断

v6.2 的施工不应该追求“所有模型全开”。真正的优先级是：

1. 数据和评测可靠；
2. 模块和 barrier 可辨识；
3. SRB latent 可解释；
4. repair direction 可验证；
5. spatial niche 能裁判；
6. ZSL 能证明机制空间具备组合外推；
7. evidence card 能限制 claim；
8. manuscript story 能在失败分支下仍然成立。

这条路线保留方法学野心，但不让项目死于方法学膨胀。
SRB 是主干，SNGM 是空间裁判，ZSL 是外推层，evidence card 是 claim 控制器。
任何一个模型失败，都应该降级，而不是拖着全项目一起殉葬。
