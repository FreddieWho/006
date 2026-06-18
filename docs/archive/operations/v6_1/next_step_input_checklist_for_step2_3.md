# Step2-3 输入清单（仅输入检查，不执行训练）

## 0. 说明
为避免序号歧义，本清单同时标注两套编号：
- 实施序号：Step2（免疫状态测量准备）-> Step3（强基线与模块发现准备）
- 契约序号：S1 -> S2/S3（见 `analysis_contract.md`）

本文件只做输入检查与门槛定义，不包含任何模型训练或机制推断。

## 1. Step2（对应 contract S1）输入检查
目标：确保 immune-state measurement 前的数据入口稳定。

必需输入文件：
1. `results/v6_1/step1/cohort_registry_v6_1.csv`
2. `results/v6_1/step1/sample_metadata_master_v6_1.csv`
3. `results/v6_1/step1/patient_metadata_master_v6_1.csv`
4. `results/v6_1/step1/treatment_context_flags.csv`
5. `results/v6_1/step1/dataset_role_assignment.csv`
6. `mvp/sidecars/*/obs.parquet` 与 `mvp/sidecars/*/var.parquet`（按 cohort 可用性）
7. `data/processed/srt/raw/*.h5ad`（本地主数据源）

字段级检查（必须通过）：
1. `sample_metadata` 必须包含并仅使用：
- `timepoint_normalized in {pre,on_treatment,post,long_post,unknown}`
- `response_binary in {responder,non_responder,unknown,not_applicable}`
2. `treatment_context_flags.confidence_level in {high,medium,low}`。
3. `sample_id` 不得为空；`unknown sample_id` 行必须可追溯并在日志中解释。
4. 主训练候选集合必须满足：
- `scRNA_available=True`
- `treatment_context in {PD1_ICI_anchor, PD1X_extension}`
- `usable_in_main_analysis=True`

本轮状态：
- 已具备：文件存在，枚举值合法。
- 待补齐：`results/v6_1/step0/*` 合同路径尚未显式对齐到当前 `step1` 输出。

Step2 施工前 gate：
1. 建立 `step0 -> step1` 路径映射或软链接策略并写入 manifest。
2. 固化 patient-group split 规则（写成可复用 split table）。
3. 固化 cohort inclusion table（哪些 cohort 进入主分析，哪些仅 sensitivity/external）。

## 2. Step3（对应 contract S2/S3）输入检查
目标：在不训练模型的前提下，检查是否具备进入强基线/模块发现的输入前提。

必需输入文件（来自 Step2 输出或前置主表）：
1. `results/v6_1/step1/sample_metadata_master_v6_1.csv`
2. `results/v6_1/step1/patient_metadata_master_v6_1.csv`
3. `results/v6_1/step1/treatment_context_flags.csv`
4. `results/v6_1/step1/data_leakage_risk_log.md`
5. `results/v6_1/step1/source_records_local.tsv`
6. `results/v6_1/step1/source_records_web_verified.tsv`
7. `results/v6_1/step1/dataset_role_assignment.csv`

标签可用性门槛（建议硬门槛）：
1. `PD1_ICI_anchor` 中需存在可监督 response 标签（responder/non_responder），否则仅允许做无监督表征，不得做 anchor 监督比较。
2. 对 `response_raw` 非空但 `response_binary=unknown` 的样本必须先完成映射修复。
3. `not_applicable` 样本禁止进入监督标签训练。
4. `high_confounding_support` 与 `external_anchor_only` 样本默认不进入主训练标签。

角色覆盖门槛：
1. HCC deep-dive、shared module、bulk external anchor、spatial adjudication、perturb prior 各角色必须给出“声明 cohort 列表 + 已入库 cohort 列表 + 缺口列表”。
2. perturb prior 若仍为“声明存在但 sample 表无入库”，只能在 Step3 标记为 optional，不可写入核心证据链。

泄漏控制门槛：
1. 同 patient 多样本必须同折。
2. paired pre/post 必须同折。
3. bulk external anchor 与 scRNA 训练样本必须隔离。
4. 不得使用 post-only cohort 定义基线敏感性标签。

## 3. 当前阻断项（Step3 前必须处理）
1. `PD1_ICI_anchor` response 标签可用性不足。
2. response 映射仍存在大量 unknown（尤其 raw 非空样本）。
3. registry 与 sample 实体覆盖差距较大（声明与入库不一致）。
4. 合同路径 (`step0`) 与实际输出路径 (`step1`) 尚未统一。

## 4. 建议的最小输入冻结包（用于下一步施工）
1. `frozen_cohort_inclusion_v6_1.csv`：明确每个 cohort 的 `main/sensitivity/external/excluded`。
2. `frozen_patient_split_v6_1.csv`：patient-group split + fold id。
3. `response_mapping_dictionary_v6_1.csv`：`response_raw -> response_binary` 的可审计字典。
4. `treatment_context_override_log_v6_1.csv`：人工修订项与证据来源。
5. `step2_input_manifest_v6_1.yaml`：对齐 contract 与实际路径的输入清单。
