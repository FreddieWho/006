# PD1_anchor_decision_note (v6.1 Task5)

## 1. 输入审计
- 已使用：
  - `results/v6_1/step1/PD1_anchor_label_availability_report.md`
  - `results/v6_1/step1/frozen_cohort_inclusion_v6_1.csv`
  - `results/v6_1/step1/frozen_patient_split_v6_1.csv`
  - `results/v6_1/step1/sample_metadata_master_v6_1.csv`
- 未找到：
  - `response_repair_report.md`
- 替代依据：
  - `results/v6_1/step1/response_mapping_dictionary_v6_1.csv`
  - `results/v6_1/step1/sample_metadata_master_v6_1.csv`
  - `results/v6_1/step1/patient_metadata_master_v6_1.csv`

## 2. 当前可监督样本量

以 `PD1_ICI_anchor` 且 `inclusion_status=main` 的 cohort 为准，当前可执行 anchor cohort 只有：
- `GSE206325`
- `TASK01`

### 2.1 Sample-level label availability
- anchor samples total: `360`
- repaired binary labels:
  - `responder=130`
  - `non_responder=210`
  - `unknown=20`

分 cohort：
- `GSE206325`
  - `responder=92`
  - `non_responder=176`
  - `unknown=20`
  - 但全部为 `post`
- `TASK01`
  - `responder=38`
  - `non_responder=34`
  - 含 `pre/on_treatment/post`

### 2.2 Patient-level label availability
- executable anchor patients total: `71`
- patient best-response:
  - `responder=21`
  - `non_responder=30`
  - `unknown=20`

分 cohort：
- `GSE206325`
  - `responder=7`
  - `non_responder=17`
  - `unknown=20`
  - 问题：不是 clean baseline，因为全部为 post-treatment / pathological-response context
- `TASK01`
  - `responder=14`
  - `non_responder=13`
  - 可用于 baseline response supervision

### 2.3 Truly baseline-supervisable anchor set
若要求“pre-treatment + known responder/non_responder + patient-level split-compatible”，当前只剩：
- cohort: `TASK01`
- patients: `25`
- label distribution:
  - `responder=13`
  - `non_responder=12`

按 frozen split 分布：
- `train=16`
- `validation=5`
- `test=4`

这个规模只来自单一 cohort，不足以作为 robust supervised PD1 anchor 主链。

## 3. 可用 cohort 判断

### 3.1 可作为 weak supervision 的 cohort
- `TASK01`
  - 优点：`mRECIST`、有 pre/on/post、pre-treatment 患者标签完整
  - 局限：单 cohort、样本量小、泛化能力弱

### 3.2 只能作为 response-associated / auxiliary anchor 的 cohort
- `GSE206325`
  - 优点：`R/NR` 已解析，标签并非缺失
  - 局限：全部 `post`，且 endpoint 是 `PATHOLOGICAL_RESPONSE`
  - 结论：可用于 anchor-associated module consistency check，不能作为 clean baseline supervised anchor

### 3.3 当前不可执行但声明为 anchor 的 cohort
- `GSE123813`
- `GSE176021`
- `GSE272734`

这些 cohort 仍是 `registry_only` / `excluded_pending_metadata`，当前不能贡献监督标签。

## 4. 推荐路线

推荐路线：`2. weakly supervised PD1/ICI anchor`

不推荐 `1. supervised PD1 anchor` 的原因：
- 真正 baseline-supervisable 的 patient 只有 `25`
- 全部来自单一 cohort `TASK01`
- 无法支撑强 baseline benchmark、跨 cohort 稳健性、或主结论级别的 supervised anchor claim

不选择 `3. unsupervised immune-state anchor` 作为主路线的原因：
- 当前并非完全无监督可用
- `TASK01` 已提供小规模但真实的 pre-treatment supervision
- 完全放弃 supervision 会损失 anchor 定向性

不选择 `4. 暂缓 PD1 anchor，优先 HCC-specific mechanism` 作为当前主路线的原因：
- anchor 不是不可用，而是“弱可用”
- 更合理的做法是保留 anchor，但降级为 weak supervision / auxiliary evaluation
- 同时允许 HCC-specific branch 在机制层面承担更高权重

## 5. 执行解释

当前最合理的 anchor 定义是：
- 用 `TASK01 pre-treatment` 作为小规模 weak supervision anchor
- 用 `GSE206325 post/pathology` 作为 anchor-associated support，而非 baseline supervision
- 用其他 `PD1X_extension` main cohorts 做 shared immune-state/module discovery

因此 Step3 不应把 PD1 anchor 设计成“先训练一个强监督响应分类器，再向外推广”的路线。
更合理的是：
- 先做 immune-state / shared-specific module discovery
- 再检查这些模块是否在 `TASK01 pre-treatment` 中区分 responder vs non_responder
- 再检查这些模块是否在 `GSE206325` 中呈现 responder-associated consistency

## 6. 对 Step3 的影响

Step3 应调整为：
- `PD1 anchor` 只作为 weak supervision gate，不作为唯一 module discovery driver
- module discovery 仍以 multi-cohort main scRNA 数据为主
- `TASK01 pre-treatment` 仅用于：
  - anchor directionality check
  - responder/non_responder ranking sanity check
  - baseline comparison support
- `GSE206325` 仅用于：
  - post-treatment / pathological-response consistency support
  - 不得与 `TASK01 pre-treatment` 混为同一监督任务

Step3 不应输出的结论：
- “跨 cohort robust supervised PD1 anchor classifier”
- “基于多个 clean PD1 monotherapy cohort 的稳定 responder module”

Step3 可以输出的结论：
- “weakly supervised PD1/ICI anchor-consistent module”
- “anchor-associated immune-state axis”

## 7. 对 Step5 的影响

Step5 的空间/组织裁判应当服务于：
- HCC-specific barrier module
- PD-1+X repair logic
- weak anchor-consistent immune-state module

Step5 不应承担的任务：
- 替代缺失的 supervised anchor 证据
- 用空间图强行证明一个本来监督基础不足的 pan-cancer PD1 anchor

Step5 更合适的使用方式：
- 优先裁判 HCC-specific mechanism 是否有组织生态意义
- 检查 weak anchor-consistent module 是否具备 spatial plausibility
- 如果空间支持强，可以提升 L2 mechanistic support；但不能把弱监督直接提升成 robust supervised anchor

## 8. 最终决策

- 当前可监督 baseline anchor patient 数：`25`
- 当前 baseline label 分布：`responder=13`, `non_responder=12`
- 当前可执行 anchor cohort：
  - weak supervision: `TASK01`
  - auxiliary support only: `GSE206325`
- 推荐路线：`weakly supervised PD1/ICI anchor`

一句话总结：
当前 PD1 anchor 不是“不可做”，但只能作为 `weak supervision + auxiliary consistency support`，不能作为 v6.1 主链中的强监督锚点。
