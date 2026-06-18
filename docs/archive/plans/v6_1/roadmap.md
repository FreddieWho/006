# v6.1 分步施工总览

整体分为 11 步：

0. 项目控制层与范围锁定
1. 数据资产盘点、cohort registry 与 metadata 主表
2. Immune-state measurement
3. 强基线包
4. Shared/HCC-specific immune module discovery
5. PD-1 monotherapy anchor 建立
6. PD-1+X repair logic 构建
7. Spatial/tissue adjudication
8. Perturbation support 与 mICS/MRI 计算
9. Bulk/external anchor
10. Evidence cards 与 nomination package
11. 论文主图与复现包整理

中间设置 5 个合练里程碑：

* **Milestone A：Scope/Data 合练**，完成 Step 0-1 后
* **Milestone B：Baseline 合练**，完成 Step 2-3 后
* **Milestone C：Mechanism 合练**，完成 Step 4-6 后
* **Milestone D：Evidence 合练**，完成 Step 7-9 后
* **Milestone E：Paper Story 合练**，完成 Step 10-11 后

---

# Step 0：项目控制层与范围锁定

## 任务目的

把 v6.1 的研究目标、数据边界、治疗结构、命名规范和不可做事项先锁死。这个步骤不是写漂亮话，而是为了防止后续 agent 一看到“空间、扰动、PD-1+X、大模型”就开始像受惊的章鱼一样往八个方向喷墨。

## 目标

形成一个可执行的 `scope_lock_v6_1.md`，明确：

* 主问题：PD-1-centered TME mechanism nomination
* 主训练逻辑：PD-1/ICI monotherapy anchor → PD-1+X repair logic
* 主数据：ICI scRNA / bulk external / spatial tissue / perturb prior
* 主输出：shared module、HCC-specific module、PD-1 anchor、X repair mechanism、candidate target
* 不做：临床推荐系统、网页、全模态大一统模型、RL 顺序用药主线、空间 GNN 主模型、未验证 counterfactual 药物推荐
* MVP 继承规则：只继承数据组织和表示底座，不继承手工 MOA 为主证据

## 建议输出

* `docs/v6_1/scope_lock_v6_1.md`
* `docs/v6_1/decision_log.md`
* `docs/v6_1/analysis_contract.md`
* `docs/v6_1/output_manifest_template.yaml`
* `docs/v6_1/terminology_dictionary.md`

## Step 0 详细 prompt

```text
你现在是 v6.1 项目的总控 agent。你的任务不是开始计算，而是锁定研究范围、输出契约和执行边界。

项目总目标：
构建一个 PD-1-centered TME mechanism nomination framework。核心逻辑是：先利用 PD-1/ICI monotherapy 或较干净 ICI 队列学习基础免疫敏感性与原发耐药模块，再在 PD-1+X 队列或外部证据中解释 X 修复了什么残余 TME 障碍，最终输出可被空间/组织证据、扰动先验、bulk 外部锚点和最小实验支持的 PD-1+X 机制假设与候选靶点。

请严格遵循以下原则：

1. 治疗结构采用 PD-1-centered continuum，而不是把 PD-1 单药和 PD-1+X 拆成两个独立课题。
   - PD-1/ICI monotherapy：学习基础 ICI 敏感性与耐药轴。
   - PD-1+X：解释在 PD-1 基础上还需要修复什么 TME 障碍。
   - 高混杂联合治疗：只做敏感性分析或外部支持，不进入主训练标签核心。

2. 证据等级只使用三档：
   - L1 Robust association：跨队列、跨中心、强基线、稳健性支持。
   - L2 Mechanistic support：pathway/TF、perturb prior、ligand-receptor、spatial/tissue 支持。
   - L3 Translational nomination：外部 bulk/clinical anchor、drugability、spatial adjudication、最小实验可行性支持。

3. 强基线必须前置，任何复杂模型都必须和强基线比较。强基线至少包括：
   - cell fraction baseline
   - curated immune signature baseline
   - pseudobulk pathway baseline
   - elastic net / logistic / Cox / XGBoost baseline
   - static module baseline
   - scVI latent + simple classifier baseline

4. 空间/组织数据是 biological adjudicator，不是主模型装饰。空间层主要回答：
   - T cell exclusion 是否成立
   - myeloid-tumor neighborhood 是否支持核心模块
   - APC niche 是否与 responder-like module 相关
   - LRG 是否有空间共定位基础
   - HCC-specific barrier 是否具有组织生态意义

5. MVP 只允许审慎继承：
   - 可以继承 cohort registry、sidecar-first、patient/sample metadata、pseudobulk、scVI atlas 作为 QC/表征底座、patient-level 表示经验、scope lock 与输出契约。
   - 只能降级继承 archetype、trajectory、axis、TME-MOA、virtual drug shift，作为展示层、baseline 或 sanity check。
   - 不允许把手工 MOA、heuristic counterfactual、pseudo-archetype bulk assignment、demo-oriented ranking 作为主结果。

6. 禁止事项：
   - 不搭建网页。
   - 不做临床推荐系统。
   - 不把 LLM/agent 放进核心建模。
   - 不强行使用 foundation model、GNN、RL。
   - 不用 LINCS reversal 单独提名药物。
   - 不凭 marker 均值单独提出机制结论。
   - 不把空间图当作装饰性验证。

请输出以下文件内容草案：
1. docs/v6_1/scope_lock_v6_1.md
2. docs/v6_1/decision_log.md
3. docs/v6_1/analysis_contract.md
4. docs/v6_1/output_manifest_template.yaml
5. docs/v6_1/terminology_dictionary.md

scope_lock_v6_1.md 必须包含：
- project title
- main scientific question
- primary analysis
- secondary analysis
- excluded analysis
- treatment continuum definition
- evidence levels
- MVP inheritance policy
- main output list
- closure/fallback logic

analysis_contract.md 必须包含：
- 每一步允许读取什么输入
- 每一步必须输出什么表
- 后续步骤不得直接绕过哪些中间表
- 哪些旧 MVP 文件只能做 proxy / sanity check
- 哪些结果不得出现在主结论中

output_manifest_template.yaml 必须为后续所有步骤提供统一字段：
- step_id
- input_files
- output_files
- scripts
- parameters
- random_seed
- cohort_scope
- notes
- downstream_dependencies

terminology_dictionary.md 必须统一以下术语：
- PD-1 anchor
- PD-1+X repair logic
- shared immune module
- HCC-specific barrier module
- mIMS
- MRI
- mICS
- ESR
- LRG-spatial
- biological adjudication
- translational nomination
```

---

# Step 1：数据资产盘点、cohort registry 与 metadata 主表

## 任务目的

把所有可用队列从“我大概有一些数据”整理成可训练、可分层、可追踪的数据资产表。这个阶段不追求分析结果，只追求别把样本、患者、治疗、时间点、响应标签搞成一锅临床粥。

## 目标

建立三个主表：

* `cohort_registry_v6_1.csv`
* `sample_metadata_master_v6_1.csv`
* `patient_metadata_master_v6_1.csv`

同时建立治疗上下文和数据用途标记：

* `treatment_context_flags.csv`
* `dataset_role_assignment.csv`
* `data_leakage_risk_log.md`

## 关键字段

`cohort_registry_v6_1.csv`：

* cohort_id
* disease
* cancer_group
* source_type
* data_modality
* sample_type
* platform
* center_or_study
* treatment_regimen_raw
* treatment_context
* timepoint_schema
* response_label_type
* response_label_quality
* paired_available
* spatial_available
* bulk_available
* perturb_available
* usable_for_PD1_anchor
* usable_for_PD1X_extension
* usable_for_HCC_specific
* usable_for_shared_module
* usable_for_external_anchor
* usable_for_spatial_adjudication
* exclusion_reason

`treatment_context` 建议只用四类，别又开始造分类学帝国：

* `PD1_ICI_anchor`
* `PD1X_extension`
* `high_confounding_support`
* `external_anchor_only`

## Step 1 详细 prompt

```text
你现在执行 v6.1 Step 1：数据资产盘点、cohort registry 与 metadata 主表构建。

本步骤只做数据治理，不做模型训练，不做机制解释，不做图。请优先保证 patient/sample/timepoint/treatment/response 的一致性。任何 metadata 不清楚的地方必须显式标记 unknown 或 ambiguous，不允许凭感觉补全。你不是算命摊，虽然很多公共数据集确实像命理材料。

项目治疗结构采用 PD-1-centered continuum：

1. PD1_ICI_anchor
   用于学习基础 ICI/PD-1 敏感性和耐药轴。
   包括 PD-1 单药、PD-L1 单药、较干净 ICI 队列。
   CTLA-4 可作为辅助 ICI 参考，但不得混淆为 PD-1 主轴。

2. PD1X_extension
   用于解释在 PD-1 基础上还需要修复什么 TME 障碍。
   包括 PD-1 + anti-VEGF/TKI、PD-1 + 其他较明确联合治疗。
   如果是 locoregional therapy + ICI，需要标记为谨慎使用。

3. high_confounding_support
   用于敏感性分析或支持性解释，不进入主训练标签核心。
   包括 TACE/HAIC/radiotherapy 强混杂联合治疗、多线治疗后、timepoint 不清楚 post-treatment、治疗方案过度混杂队列。

4. external_anchor_only
   用于 bulk/spatial/clinical 外部锚定，不参与主模型训练。

请完成以下任务：

A. 扫描并整理所有可用数据源
- 包括当前项目已有 MVP 数据、HCC 主干数据、多癌种 ICI reference、bulk external anchor、spatial/tissue 数据、perturb prior。
- 不要求下载新数据，但如果发现缺失，请记录在 missing_data_request.md。
- 对每个数据源判断其角色，而不是只列名字。

B. 建立 cohort_registry_v6_1.csv
每个 cohort 至少包含以下字段：
- cohort_id
- cohort_name
- disease
- cancer_group
- source_type
- data_modality
- sample_type
- platform
- center_or_study
- raw_data_path
- processed_data_path
- sidecar_path
- clinical_metadata_path
- treatment_regimen_raw
- treatment_context
- timepoint_schema
- response_label_type
- response_label_quality
- paired_available
- spatial_available
- bulk_available
- perturb_available
- usable_for_PD1_anchor
- usable_for_PD1X_extension
- usable_for_HCC_specific
- usable_for_shared_module
- usable_for_external_anchor
- usable_for_spatial_adjudication
- exclusion_reason
- notes

C. 建立 sample_metadata_master_v6_1.csv
每个 sample 至少包含以下字段：
- cohort_id
- patient_id
- sample_id
- original_sample_id
- disease
- tissue_source
- sample_type
- timepoint_raw
- timepoint_normalized
- treatment_regimen_raw
- treatment_context
- response_raw
- response_binary
- response_ordered
- survival_available
- spatial_available
- bulk_available
- scRNA_available
- usable_in_main_analysis
- usable_in_sensitivity_analysis
- exclusion_reason

timepoint_normalized 只允许：
- pre
- on_treatment
- post
- long_post
- unknown

response_binary 只允许：
- responder
- non_responder
- unknown
- not_applicable

D. 建立 patient_metadata_master_v6_1.csv
每个 patient 至少包含：
- cohort_id
- patient_id
- disease
- treatment_context_primary
- has_pre_sample
- has_on_treatment_sample
- has_post_sample
- has_paired_pre_post
- best_response_raw
- best_response_binary
- survival_available
- number_of_samples
- number_of_modalities
- use_for_PD1_anchor
- use_for_PD1X_extension
- use_for_HCC_specific
- use_for_external_anchor
- notes

E. 建立 treatment_context_flags.csv
每一行是 sample 或 patient 的治疗上下文标记，至少包含：
- cohort_id
- patient_id
- sample_id
- treatment_regimen_raw
- treatment_context
- confidence_level
- reason_for_assignment
- ambiguity_note

confidence_level 只允许：
- high
- medium
- low

F. 建立 dataset_role_assignment.csv
每个 cohort 必须被分配至少一个角色：
- main_scRNA_training
- PD1_anchor_training
- PD1X_extension_support
- HCC_specific_deep_dive
- spatial_adjudication
- bulk_external_anchor
- perturb_prior
- sensitivity_only
- excluded

G. 建立 data_leakage_risk_log.md
必须检查并记录：
- 同一 patient 多样本是否可能跨 train/test
- 同一 cohort/center 是否支配某一 response group
- pre/post 是否会泄漏 response
- bulk external anchor 是否与 scRNA training 数据重叠
- spatial 数据是否来自同一患者或同一 cohort
- response label 是否由后续治疗影响
- therapy_context 是否混杂

H. 更新 output_manifest
为本步骤记录：
- 输入文件
- 输出文件
- 路径
- 主要字段
- 已知缺失
- 下游依赖

本步骤不要做：
- 不做 cell annotation
- 不做 pseudobulk
- 不做 differential expression
- 不做 baseline model
- 不做 module discovery
- 不做任何药物或机制推荐

请最终输出：
1. cohort_registry_v6_1.csv
2. sample_metadata_master_v6_1.csv
3. patient_metadata_master_v6_1.csv
4. treatment_context_flags.csv
5. dataset_role_assignment.csv
6. data_leakage_risk_log.md
7. missing_data_request.md
8. step1_summary.md
```

---

# Milestone A：Scope/Data 合练

## 合练时机

完成 Step 0 和 Step 1 后。

## 合练目的

确认项目范围、数据角色、治疗上下文、metadata 结构可以支撑后续所有步骤。这个合练不产生生物结论，只检查“地基有没有歪”。考虑到人类项目最爱在地基歪掉后开始讨论屋顶颜色，这一步非常必要。

## 合练 prompt

```text
请执行 v6.1 Milestone A：Scope/Data 合练。

输入：
- scope_lock_v6_1.md
- analysis_contract.md
- cohort_registry_v6_1.csv
- sample_metadata_master_v6_1.csv
- patient_metadata_master_v6_1.csv
- treatment_context_flags.csv
- dataset_role_assignment.csv
- data_leakage_risk_log.md

任务：
1. 检查 v6.1 的科学目标是否能由当前数据结构支撑。
2. 检查 PD1_ICI_anchor、PD1X_extension、high_confounding_support、external_anchor_only 四类治疗上下文是否足够清晰。
3. 检查 HCC deep-dive、多癌种 shared module、bulk external anchor、spatial adjudication、perturb prior 是否都有对应数据角色。
4. 列出当前最可能阻断后续施工的 5 个数据风险。
5. 给出 Step 2-3 的输入清单，不做任何模型训练。

输出：
- milestone_A_scope_data_rehearsal.md
- next_step_input_checklist_for_step2_3.md
```

---

# Step 2：Immune-state measurement

## 任务目的

把单细胞数据整理成 patient/sample 层面的 TME 状态表。这里的目标不是发现机制，而是稳定测量：细胞组成、pseudobulk、signature/pathway、基础 immune state。

## 目标

输出后续强基线和模块发现的基础矩阵：

* `cell_state_annotation_v6_1.csv`
* `cell_fraction_by_sample.csv`
* `pseudobulk_by_sample.parquet`
* `immune_signature_scores.csv`
* `tf_pathway_activity_scores.csv`
* `myeloid_QC_report.md`

## Prompt

```text
执行 v6.1 Step 2：Immune-state measurement。

基于 Step 1 的 metadata，整合可用 scRNA 数据，生成 sample/patient 层面的 immune-state measurement 表。优先复用已有注释和 MVP 的 sidecar/pseudobulk 资产，但必须重新检查髓系状态，不允许把旧 TME-MOA 当主证据。

输出 cell fraction、pseudobulk、curated immune signature scores、pathway/TF activity scores，并单独生成 myeloid_QC_report.md。scVI/atlas 可作为 QC 和表征底座，不作为主模型结果。
```

---

# Step 3：强基线包

## 任务目的

建立 v6.1 的最低可防守结果。之后所有复杂模型都必须和它比。没有强基线，复杂模型就是穿西装的噪声。

## 目标

输出一套 baseline comparison：

* fraction baseline
* curated signature baseline
* pseudobulk pathway baseline
* elastic net / logistic / Cox / XGBoost baseline
* static module baseline
* scVI latent baseline

## Prompt

```text
执行 v6.1 Step 3：Strong baseline first。

基于 Step 2 的 cell fraction、signature scores、pathway/TF activity、pseudobulk 和 scVI latent，建立强基线包。任务不是追求最高分，而是建立无泄漏、可复现、可解释的 baseline reference。

至少包含：
1. cell fraction baseline
2. curated immune signature baseline
3. pseudobulk pathway baseline
4. elastic net / logistic / Cox / XGBoost baseline
5. static module baseline
6. scVI latent + simple classifier baseline

所有 baseline 使用相同 split、相同 treatment_context 定义、相同 response mapping。输出 baseline_results.csv、baseline_model_manifest.yaml、baseline_comparison_report.md。
```

---

# Milestone B：Baseline 合练

## 合练时机

完成 Step 2 和 Step 3 后。

## 合练目的

确认 immune-state measurement 与强基线能形成一条自洽链路。这个阶段要判断：现在是否已经有足够强的 baseline 支撑后续复杂分析，还是其实简单模型已经把能榨的都榨干了。

## 合练 prompt

```text
执行 v6.1 Milestone B：Baseline 合练。

输入 Step 1-3 全部输出，完成一次从 cohort registry → metadata → immune-state measurement → strong baseline 的完整串联检查。

请输出：
1. 当前数据能否支持 PD-1 anchor 分析。
2. 哪类 baseline 最强。
3. baseline 是否主要在学习 center/cohort/batch。
4. HCC 与 pan-cancer 的表现是否分离。
5. 后续 module discovery 应优先使用哪些输入矩阵。
6. 需要暂时排除或降权的 cohort/sample。

输出 milestone_B_baseline_rehearsal.md。
```

---

# Step 4：Shared/HCC-specific immune module discovery

## 任务目的

从强基线和 pseudobulk/signature/pathway 层进入机制模块发现。v6.1 在这里不急着宣称 causal GRN，而是先找到稳定、可解释、跨队列可复现的 shared module 和 HCC-specific module。

## 目标

输出：

* `shared_immune_modules.csv`
* `hcc_specific_modules.csv`
* `module_score_matrix.csv`
* `module_response_association.csv`
* `module_stability_report.md`

## Prompt

```text
执行 v6.1 Step 4：Shared/HCC-specific immune module discovery。

基于 pseudobulk、immune signature、pathway/TF activity 和 baseline 结果，发现跨癌种 shared immune modules 与 HCC-specific barrier modules。不要直接宣称 causal GRN，先输出 inferred immune modules，并记录跨队列、跨中心、LOCO 或 leave-center-out 稳定性。

重点寻找：
1. shared PD-1 sensitivity module
2. shared primary resistance module
3. HCC-specific myeloid/Treg/stromal/vascular barrier module
4. antigen presentation / IFN deficient module

输出 module_score_matrix、shared_immune_modules、hcc_specific_modules、module_stability_report。
```

---

# Step 5：PD-1 monotherapy anchor 建立

## 任务目的

用较干净的 PD-1/ICI 单药或 anchor 队列建立基础免疫敏感性和原发耐药轴。后续 PD-1+X 的所有故事都必须接在这个 anchor 上，不能凭空发明一个 X 修复宇宙。

## 目标

输出：

* `PD1_anchor_modules.csv`
* `PD1_sensitivity_axis.csv`
* `PD1_primary_resistance_axis.csv`
* `PD1_anchor_patient_scores.csv`
* `PD1_anchor_report.md`

## Prompt

```text
执行 v6.1 Step 5：PD-1 monotherapy anchor。

基于 Step 4 的模块和 Step 1 的 treatment_context，只使用 PD1_ICI_anchor 或较干净 ICI 队列建立 PD-1/ICI 基础敏感性与原发耐药模块。输出 responder-like module、primary resistance module、immune-cold/APC-poor/myeloid-suppressed 等障碍类型。

不要把 PD-1+X 队列混入 anchor 主训练；PD-1+X 只能作为后续 extension support。
```

---

# Step 6：PD-1+X repair logic 构建

## 任务目的

在 PD-1 anchor 的基础上解释：为什么还需要 X，X 可能修复什么 TME 障碍。这里不是直接推荐药，而是把 residual barrier 映射到机制类。

## 目标

输出：

* `PD1X_repair_logic_table.csv`
* `residual_barrier_modules.csv`
* `X_mechanism_class_mapping.yaml`
* `PD1_anchor_to_X_repair_report.md`

## Prompt

```text
执行 v6.1 Step 6：PD-1+X repair logic。

基于 PD-1 anchor modules 和 HCC/shared resistance modules，构建 PD-1+X repair logic。对每个 residual barrier 判断其对应的 X mechanism class，包括 myeloid reprogramming、vascular/stromal remodeling、antigen presentation/IFN restoration、T-cell co-stimulation、epigenetic immune sensitization。

PD-1+X 队列用于支持“X 修复什么”，不能变成另一个完全独立的主问题。
```

---

# Milestone C：Mechanism 合练

## 合练时机

完成 Step 4-6 后。

## 合练目的

确认项目已经从“测量状态”进入“机制链条”：shared/HCC-specific module → PD-1 anchor → residual barrier → X repair class。这里要把故事第一次完整讲一遍，不需要空间和扰动证据都齐，但主逻辑必须闭合。

## 合练 prompt

```text
执行 v6.1 Milestone C：Mechanism 合练。

输入 Step 0-6 的全部输出，形成第一次完整机制链条演练。

请回答：
1. 当前最稳定的 shared PD-1 sensitivity module 是什么。
2. 当前最有价值的 HCC-specific barrier module 是什么。
3. PD-1 anchor 是否稳定。
4. 哪些 residual barrier 可以自然连接到 PD-1+X repair logic。
5. 哪些 X mechanism class 目前只是猜想，缺少空间/扰动/外部证据。
6. 当前论文故事更像 pan-cancer shared mechanism，还是 HCC-first translational story。

输出 milestone_C_mechanism_rehearsal.md 和 candidate_core_story_v1.md。
```

---

# Step 7：Spatial/tissue adjudication

## 任务目的

用空间/组织数据裁判核心 TME 机制是否真实存在于组织生态里。空间数据在 v6.1 里不是装饰，而是判断机制能否升格的重要证据。

## 目标

输出：

* `spatial_module_scores.csv`
* `myeloid_tumor_neighborhood_scores.csv`
* `T_cell_exclusion_scores.csv`
* `APC_niche_scores.csv`
* `LRG_spatial_support.csv`
* `spatial_adjudication_report.md`

## Prompt

```text
执行 v6.1 Step 7：Spatial/tissue adjudication。

基于可用 spatial / tissue-level 数据，对 Step 4-6 的核心模块进行空间裁判。重点分析 T cell exclusion/infiltration、myeloid-tumor neighborhood、APC niche、LRG spatial support、vascular/stromal barrier。

不要把空间数据塞进主模型；把它作为 biological adjudicator。输出 spatial_adjudication_report，并标记哪些机制获得空间支持、哪些需要降级。
```

---

# Step 8：Perturbation support 与 mICS/MRI 计算

## 任务目的

用 perturb prior、druggability、module centrality、spatial support、external consistency 来计算简化版 mICS/MRI。v6.1 不做复杂反事实炫技，先做可解释、可执行的 intervention candidate score。

## 目标

输出：

* `mIMS_scores.csv`
* `MRI_scores.csv`
* `mICS_scores.csv`
* `perturbation_support_table.csv`
* `candidate_target_priority.csv`

## Prompt

```text
执行 v6.1 Step 8：Perturbation support 与 mICS/MRI。

基于 Step 4-7 的模块、PD-1 anchor、repair logic 和 spatial support，整合 perturb prior、TF/pathway support、druggability、module centrality，计算 mIMS、MRI、mICS。

mICS 定义为第一阶段可执行版本，整合 module centrality、perturb direction、druggability、spatial/tissue support、external anchor consistency。不要凭 LINCS reversal 单独提名具体药物。
```

---

# Step 9：Bulk/external anchor

## 任务目的

把 bulk 外部队列从“趋势展示”升级成外部锚点验证。IMbrave150 这类数据可以服务 HCC 转化逻辑，但不能用简化投影糊弄自己，糊弄自己一向是人类科研的低成本娱乐。

## 目标

输出：

* `bulk_module_projection_scores.csv`
* `external_anchor_results.csv`
* `HCC_external_support_report.md`
* `external_consistency_flags.csv`

## Prompt

```text
执行 v6.1 Step 9：Bulk/external anchor。

将 shared modules、HCC-specific modules、PD-1 anchor scores、mIMS/MRI/mICS 投影到可用 bulk ICI 或 HCC 外部队列。重点评估 module score 与 response/survival/clinical subgroup 的方向一致性。

bulk 结果作为 external anchor，不替代空间或单细胞机制证据。输出 external_anchor_results 和 external_consistency_flags。
```

---

# Milestone D：Evidence 合练

## 合练时机

完成 Step 7-9 后。

## 合练目的

确认每个候选机制是否具备 L1/L2/L3 证据链。这个阶段的重点不是再发现新东西，而是判断哪些结果配进入主文，哪些只能滚去补充材料。科学世界残酷但节省版面。

## 合练 prompt

```text
执行 v6.1 Milestone D：Evidence 合练。

输入 Step 0-9 全部输出，对所有 candidate modules / mechanism classes / targets 做证据链合练。

请为每个候选生成初版 evidence level：
- L1 Robust association
- L2 Mechanistic support
- L3 Translational nomination

每个候选必须检查：
1. PD-1 anchor association
2. PD-1+X repair hypothesis
3. strong baseline comparison
4. cross-cohort robustness
5. spatial/tissue support
6. perturbation support
7. bulk/external support
8. candidate target feasibility

输出：
- milestone_D_evidence_rehearsal.md
- candidate_module_evidence_matrix.csv
- candidate_drop_or_keep_decision.md
```

---

# Step 10：Evidence cards 与 nomination package

## 任务目的

把最终候选从一堆表格整理成可以汇报、写论文、设计实验的 evidence cards。这里的输出不是“临床推荐”，而是 mechanism-prioritized PD-1+X nomination。

## 目标

输出：

* `PD1X_nomination_cards.csv`
* `top_modules_summary.md`
* `top_targets_summary.md`
* `minimal_wetlab_validation_plan.md`
* `final_evidence_table.csv`

## Prompt

```text
执行 v6.1 Step 10：Evidence cards 与 nomination package。

基于 Milestone D 的证据矩阵，为保留的 top modules、top X mechanism classes 和 top candidate targets 生成 evidence cards。每张 card 必须包含：module name、biological interpretation、PD-1 anchor association、PD-1+X repair hypothesis、strong baseline comparison、spatial/tissue support、perturbation support、bulk/external support、candidate target、experimental feasibility、evidence level。

输出 final_evidence_table、PD1X_nomination_cards、minimal_wetlab_validation_plan。命名为 nomination，不要写 clinical recommendation。
```

---

# Step 11：论文主图与复现包整理

## 任务目的

把所有结果压成 v6.1 的论文/答辩主线。不是为了现在投稿，而是形成能被阶段性审查、答辩、组会、老板质询的完整证据包。人类老板的自然栖息地就是质询场，提前准备能少死几次。

## 目标

形成 6 张主图逻辑和复现包：

* Fig1：PD-1-centered study design
* Fig2：Strong baseline and immune-state measurement
* Fig3：Shared and HCC-specific immune modules
* Fig4：Spatial adjudication of TME mechanism
* Fig5：PD-1+X repair logic
* Fig6：External validation and experimental nomination

## Prompt

```text
执行 v6.1 Step 11：论文主图与复现包整理。

基于 Step 0-10 的全部输出，整理 6 张主图的内容规划、对应数据表、关键统计结果和可视化草图说明。同步整理 reproducibility package，包括 manifest、参数、输入输出表、脚本索引、随机种子、版本信息和主要分析路径。

输出：
1. figure_plan_v6_1.md
2. main_results_narrative.md
3. reproducibility_manifest.md
4. table_figure_source_mapping.csv
5. unresolved_risks_for_next_iteration.md
```

---

# Milestone E：Paper Story 合练

## 合练时机

完成 Step 10-11 后。

## 合练目的

将所有已完成步骤合成一条完整故事线，判断当前版本是：

* HCC-first translational story；
* pan-cancer shared mechanism story；
* robust module + spatial TME mechanism story；
* 还是需要回炉重做的数据治理灾难片。

## 合练 prompt

```text
执行 v6.1 Milestone E：Paper Story 合练。

输入 Step 0-11 的全部输出，进行最终阶段合练。请不要新增分析，只根据已有结果组织故事线。

请输出：
1. 当前最适合的论文主路线：
   - HCC-first translational story
   - pan-cancer shared mechanism story
   - robust immune module + spatial TME mechanism story
   - method-oriented nomination framework
2. 六张主图每张图的核心结论。
3. 哪些结果能进主文，哪些只能进补充材料。
4. 哪些 claim 需要降级。
5. 下一轮最优先补的三个分析。
6. 最小实验验证优先级。
7. 20分钟组会汇报大纲。
8. 40分钟正式汇报大纲。

输出：
- milestone_E_paper_story_rehearsal.md
- presentation_outline_20min.md
- presentation_outline_40min.md
- next_iteration_priority.md
```

---

# 最终施工顺序压缩版

1. **Step 0**：锁 scope、术语、输出契约
2. **Step 1**：建 registry / metadata / treatment flags
3. **Milestone A**：Scope/Data 合练
4. **Step 2**：immune-state measurement
5. **Step 3**：强基线包
6. **Milestone B**：Baseline 合练
7. **Step 4**：shared/HCC-specific modules
8. **Step 5**：PD-1 anchor
9. **Step 6**：PD-1+X repair logic
10. **Milestone C**：Mechanism 合练
11. **Step 7**：spatial/tissue adjudication
12. **Step 8**：perturbation support + mICS/MRI
13. **Step 9**：bulk/external anchor
14. **Milestone D**：Evidence 合练
15. **Step 10**：evidence cards + nomination package
16. **Step 11**：主图与复现包
17. **Milestone E**：Paper Story 合练

