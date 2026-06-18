# Milestone A：Scope/Data 合练（v6.1）

## 0. 输入与审计范围
本次合练基于以下输入：
- `docs/v6_1/scope_lock_v6_1.md`
- `docs/v6_1/analysis_contract.md`
- `results/v6_1/step1/cohort_registry_v6_1.csv`
- `results/v6_1/step1/sample_metadata_master_v6_1.csv`
- `results/v6_1/step1/patient_metadata_master_v6_1.csv`
- `results/v6_1/step1/treatment_context_flags.csv`
- `results/v6_1/step1/dataset_role_assignment.csv`
- `results/v6_1/step1/data_leakage_risk_log.md`

核心统计快照：
- `cohort_registry`: 59 cohorts
- `sample_metadata`: 1,634 samples
- `patient_metadata`: 874 patients
- 四类 context（cohort 级）：`PD1_ICI_anchor=5`, `PD1X_extension=23`, `high_confounding_support=12`, `external_anchor_only=19`

## 1. 科学目标可支撑性判断
结论：**部分可支撑（conditional pass）**。

可支撑部分：
- PD-1 centered continuum 在 cohort 级已成型，四类 context 已可分层。
- 多癌种 shared module 数据角色在声明层已覆盖（28 个 main_scRNA_training cohorts）。
- HCC deep-dive 声明覆盖 11 cohorts，且有本地已入库 cohort 支撑（如 TASK01/TASK02/GSE206325/lambrecht_hcc）。
- bulk external anchor、spatial adjudication、perturb prior 在角色表上均有入口。

当前不足：
- **PD1 anchor 的 response 标签不可用度过高**：`PD1_ICI_anchor` 样本 360/360 为 `response_binary=unknown`。
- Step 合同路径存在漂移：`analysis_contract` 约定 `results/v6_1/step0/*`，当前实际主表在 `results/v6_1/step1/*`。
- registry 与 sample 实体不完全对齐：59 个 registry cohorts 中仅 22 个进入 sample 主表。

判定：
- 可进入 Step2 的数据输入整理与特征构建准备。
- 不建议直接进入依赖 PD1 anchor 标签监督的强基线比较或机制主结论。

## 2. 四类治疗上下文清晰度评估

### 2.1 PD1_ICI_anchor
- cohort 数：5
- sample 数：360
- 置信度：`medium=360`（无 `high`）
- 主要问题：response 标签全 unknown，当前仅可作为免疫状态表征入口，难以支持“敏感/耐药监督锚点”任务。

### 2.2 PD1X_extension
- cohort 数：23
- sample 数：735
- 置信度：`high=267`, `medium=468`
- 主要问题：混合方案 token 和 source_ambiguity 较多，需在 Step2 前完成 sample-level regimen 精修。

### 2.3 high_confounding_support
- cohort 数：12
- sample 数：167
- 置信度：`high=12`, `medium=86`, `low=69`
- 评价：定位基本清晰，符合“敏感性/支持性用途”。

### 2.4 external_anchor_only
- cohort 数：19
- sample 数：372
- 置信度：`high=372`
- 评价：定位清晰，适合外部锚定，不应进入主训练标签。

总体判断：
- context 定义方向正确。
- **PD1_ICI_anchor 与 PD1X_extension 的 sample 级标签质量仍需增强**，否则后续比较易失真。

## 3. 数据角色覆盖检查

声明覆盖（来自 `dataset_role_assignment.csv`）：
- `HCC_specific_deep_dive`: 11 cohorts
- `main_scRNA_training`: 28 cohorts
- `bulk_external_anchor`: 11 cohorts
- `spatial_adjudication`: 1 cohort
- `perturb_prior`: 6 cohorts

已入库 sample 表覆盖（按 cohort_id 交叉）
- HCC deep-dive：8/11 已入库
- main_scRNA_training：16/28 已入库
- bulk external anchor：3/11 已入库
- spatial adjudication：1/1 已入库
- perturb prior：0/6 已入库（仅声明，无 sample 级可用主表）

判断：
- 角色“有定义”但不完全“有可执行输入”。
- 其中 perturb prior 与 bulk external anchor 的落地程度不足，会阻断 Step3+ 的稳定推进。

## 4. 最可能阻断后续施工的 5 个数据风险

1. **PD1 anchor 标签不可监督化风险（高）**
- 现状：`PD1_ICI_anchor` 样本 response 全 unknown。
- 影响：无法完成 anchor 监督基线与稳健比较，直接阻断 v6.1 核心链路。

2. **response 映射与语义归一不完整风险（高）**
- 现状：存在大量 `response_raw` 非空但 `response_binary=unknown` 的样本（例如 NR/R/High/Medium/Low 未完全归一）。
- 影响：标签噪声导致强基线对比与模块关联显著偏差。

3. **registry-sample 漂移风险（高）**
- 现状：59 个 registry cohorts 中 37 个未进入 sample 主表。
- 影响：角色声明与可执行数据脱节，shared/HCC/external 任务的样本基数被高估。

4. **sample 级治疗与时间信息 ambiguity 过高（中高）**
- 现状：`treatment_context_flags` 中 ambiguity 非空 983/1634；`unknown_timepoint=90`。
- 影响：PD1X repair 逻辑与高混杂隔离边界不稳。

5. **合同路径与执行路径不一致风险（中）**
- 现状：合同依赖 `results/v6_1/step0/*`，当前主表输出在 `results/v6_1/step1/*`。
- 影响：下游脚本与审计流程可能误读输入，导致 no-bypass 规则失效。

## 5. Milestone A 判定与准入建议
判定：**有条件通过（Pass with blockers）**。

准入建议：
- 可以进入 Step2 的输入标准化与特征构建准备。
- 进入 Step3 前必须先完成 response 映射修复、PD1 anchor 标签可用性修复、合同路径对齐。
