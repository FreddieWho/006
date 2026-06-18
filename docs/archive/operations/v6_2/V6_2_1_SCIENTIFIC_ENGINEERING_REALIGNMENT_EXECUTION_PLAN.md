# v6.2.1 Scientific & Engineering Realignment 详细执行方案

版本：v1.0-draft
日期：2026-08-20
状态：**PLAN ONLY / NOT EXECUTED**
控制输入：`docs/v6_2/refine_1.md`

---

## 1. 方案目的

本方案用于恢复以下四者的一致性：

1. scientific plan 定义的科学问题；
2. roadmap 定义的工程顺序；
3. 当前已完成资产真正能够支持的统计对象；
4. measurement、identifiability、representation、claim 和 repair 各 gate 应承担的职责。

本方案不是新的生物分析 phase，也不以“让现有 Phase8 跑起来”为目标。其近期交付是一个可审计的架构校准包；只有该校准包获得明确批准，后续才允许进行 response-blind representation、anchor feasibility 和未来 response/context 分析。

本轮不运行模型、不读取 response 做选择、不修改 frozen membership、不调整阈值、不下载数据、不生成新 biological claim。

---

## 2. 当前冻结事实

### 2.1 当前 source of truth

科学目标由以下文件定义：

- `docs/plan/v6_2/plan.md`
- `docs/plan/v6_2/roadmap.md`

当前实际状态由以下文件定义：

- `results/v6_2/V6_2_PROJECT_SCIENTIFIC_LOGIC_AND_BLOCKER_REPORT_20260731.md`
- `results/v6_2/data_interface_repair_v1/context_repair_release_manifest.yaml`
- `results/v6_2/phase7_module_measurement_and_barrier_identifiability_repair_v1/handoff/phase7a_to_phase7b_handoff.yaml`
- `results/v6_2/phase7_module_measurement_and_barrier_identifiability_repair_v1/handoff/phase7b_to_phase8_handoff.yaml`
- `results/v6_2/phase8_anchor_repair_and_statistical_strengthening_repair_v1/phase8_repair_input_contract.yaml`
- `results/v6_2/phase8_anchor_repair_and_statistical_strengthening_repair_v1/phase8_contract_and_run_manifest.yaml`

早于 `context_repair_v1_full_rerun` 的旧 Phase8 effect、route 和 sample-surface 产物只能用于追溯漂移与风险，不能作为当前科学结果。

### 2.2 当前不可变事实

- release：`context_repair_v1_full_rerun`；
- release status：`COMPLETE_REPAIRS_PHASE8_GATE_BLOCKED`；
- Phase6 frozen membership SHA256：`285c11de4b071c45b0ec43464fb963650340c29bb2de3037570578b29ee44036`；
- 8 个 FM 均由 response-blind 流程产生；
- Phase7A/7B 数值结果保留；
- 当前仅 FM07 具备进入未来独立 mechanism-attribution 审计的前置资格；
- 当前 Phase8：`ABORTED_PRE_MODEL`；
- `response_model_run=false`；
- Phase9、repair 和 X-class 均未获准进入。

这些事实不能被本次架构校准改写。

---

## 3. 独立架构判断

### 3.1 发生了什么漂移

原 plan 的统计对象是多模块 patient immune-tumor state。FM 是 patient state 的可解释坐标，SRB Stage-A 应先以 response-blind 方式学习联合表征；response 只在后续 anchor/context head 中进入；identifiability 保护的是独立 barrier attribution、dominant-barrier claim 和 repair，不是用来先删掉联合表征中的相关信息。

当前实现逐步变成：

```text
8 个 FM
  → 逐 FM measurement gate
  → 逐 FM identifiability gate
  → 只有可独立归因 FM 才可进入 response
  → 当前只剩 FM07
  → Phase8 三-FM contract 与 FM07-only 输入冲突
```

这是 **high-severity、可恢复的 scientific/engineering drift**。measurement 和 identifiability 计算本身仍有价值；错误在于它们被赋予了过宽的下游输入淘汰权。

### 3.2 plan 与 roadmap 的内部冲突

`plan.md` 的 Barrier Identifiability Gate 明确列出 `SRB Stage-A latent` 作为输入，且 SRB 训练制度要求先用 response-blind reconstruction 学习稳定表征，再接入稀缺监督。

`roadmap.md` 却把 Phase7 identifiability 放在 Phase8 SRB Stage-A 之前，并由 Controller 冻结 eligible barrier set。后续实现进一步把该 eligible set 变成 response 模型的唯一 feature set。

校准解释如下：

- pre-representation 阶段需要的是 **measurement reliability gate**；
- response-blind Stage-A 使用的是 **representation set**；
- response/context 分析后再执行 **attribution/claim identifiability gate**；
- repair 只读取通过 attribution gate 的 independent module、co-dominant set 或 coarse barrier group。

该解释更符合原 plan 的多模块、context-dependent、potentially nonlinear patient-state 目标，也保留了 identifiability audit 对强机制 claim 的保护作用。

### 3.3 当前是否应继续 FM07-only Phase8

**不应直接继续。**

原因有两层：

1. 当前三-FM contract 与 FM07-only 输入不一致，属于 governance hard fail；
2. 直接改成 FM07-only 会把架构漂移固化为单特征 response 问题，且不能解决 anchor 在 cancer、drug、endpoint 和 cohort 上同步变化的可识别性问题。

FM07 可继续保留为当前唯一 independent-claim candidate，但不应成为 representation 的唯一输入，也不能在架构校准前启动 response analysis。

---

## 4. 校准后的科学与工程架构

### 4.1 三个集合

#### A. Scoreable set

在冻结 scoring contract 下可以计算 patient-timepoint score 的 FM。它只说明分数可生成，不说明测量可靠、可独立归因或与 response 有关。

#### B. Representation set

用于表达完整 patient immune state。允许包含：

- 高可靠 FM；
- 具有 measurement uncertainty、但仍保留结构信息的 FM；
- 彼此共线、但可作为 block/coarse group 共同表示的 FM；
- 显式 missingness、reliability weight 和 uncertainty mask。

进入 representation 不等于获得独立 mechanism claim。

#### C. Claim-eligible set

仅用于独立机制归因、独立 effect、独立 repair 和单独连接 X-class。它要求：

- measurement 可靠；
- attribution 可辨识；
- response/context 方向合法；
- replication 和 negative controls 合格；
- 不与 endpoint/cohort confounding 混淆。

### 4.2 校准后的关键顺序

```text
Frozen data/interface
  → Frozen response-blind FM membership
  → Measurement reliability + uncertainty map
  → Representation-set contract
  → Response-blind representation necessity gate
      ├─ direct module backbone
      ├─ linear/hierarchical latent
      └─ SRB Stage-A candidate
  → Anchor feasibility gate
  → Joint response/context modeling on frozen representation
  → Post-model attribution/claim identifiability
  → barrier-state verdict
      ├─ single dominant
      ├─ co-dominant
      ├─ coarse barrier group
      ├─ diffuse/mixed
      └─ unresolved/abstain
  → Repair-readiness gate
  → transport / spatial / perturbation / X-class
```

### 4.3 Response head 的边界

response head 不发现 representation，不重训 FM，不选择 module membership，也不决定 measurement threshold。它只在已经冻结的 representation 上回答：

- patient state 与 response 的关系；
- shared 与 context-modulated structure；
- endpoint/cancer/treatment interaction；
- responder-compatible direction；
- uncertainty 与 abstention。

---

## 5. 当前 FM 的预注册式重新解释

下表只重新解释已有结果，不重新筛选 FM。

| FM | 当前 measurement 结果 | 当前 identifiability 结果 | representation 候选角色 | 当前 independent-claim 资格 |
|---|---|---|---|---|
| FM01 | bootstrap 0.910；top-decile Jaccard 0.492，measurement sensitivity | VIF 2.188；可分但 measurement 未过主门 | 保留，带 measurement-uncertainty tag | 否 |
| FM02 | bootstrap 0.796；Jaccard 0.426 | 与 FM03 共线；CB02 | 以 FM02/FM03 block 或 coarse group 保留 | 否 |
| FM03 | bootstrap 0.811；Jaccard 0.557 | VIF 11.005；`blocked_nonidentifiable` | 与 FM02 联合表示，不作独立归因 | 否 |
| FM04 | bootstrap 0.864；Jaccard 0.492，measurement sensitivity | VIF 1.830；可分但 measurement 未过主门 | 保留，带 measurement-uncertainty tag | 否 |
| FM05 | bootstrap 0.835；Jaccard 0.485 | 与 FM08 高相关；CB04 | 以 FM05/FM08 block 或 coarse group 保留 | 否 |
| FM06 | bootstrap 0.826；completeness 0.706；Jaccard 0.448 | VIF 2.372；跨分辨率有限 | 保留为低完整度 sensitivity coordinate | 否 |
| FM07 | bootstrap 0.870；Jaccard 0.646 | VIF 2.188；`eligible_independent` | 独立 coordinate，仍属于联合 state | **仅具备进入未来 claim 审计的前置资格** |
| FM08 | bootstrap 0.854；Jaccard 0.469 | 与 FM05 高相关；CB04 | 以 FM05/FM08 block 或 coarse group 保留 | 否 |

所有 FM 的 biological name 继续保持 provisional；没有任何 FM 可在当前阶段称为 response/failure/resistance barrier。

---

## 6. Blocker 分解与处置

| blocker | 立即阻止代码运行 | 立即阻止强 claim | 合同 patch 可解决 | 需要新数据 | 需要重构 Phase7–9 | 当前性质 |
|---|---:|---:|---:|---:|---:|---|
| A. 三-FM contract 与 FM07-only 输入冲突 | 是 | 是 | 是，仅解决程序性冲突 | 否 | 是 | governance hard fail |
| B. identifiability 被用作 representation feature selection | 应阻止现架构继续 | 是 | 部分 | 否 | 是 | high-severity architecture drift |
| C. cancer/drug/endpoint/cohort 同步变化 | 不阻止 response-blind Stage-A | 是 | 否 | 通常是 | response/context 设计需重构 | scientific identifiability blocker |
| D. 缺独立可交换 anchor 与 HCC baseline replication | 不阻止 response-blind Stage-A | 是 | 否 | 是 | 不必返工冻结模块 | evidence/replication blocker |
| E. completeness 0.734 与书面 0.80/实现 0.60 不一致 | 阻止新 freeze | 是 | 是 | 否 | 否 | governance uncertainty |

---

## 7. 执行总览

执行分成两个 release：

### Release R0：Realignment control package

只做阅读、合同、schema、资产处置和验证器设计，不运行任何科学计算。

### Release R1：Realigned scientific execution

只有 R0 获批后启动。先进行 response-blind representation 与 anchor feasibility；response 解锁需要单独 gate。R1 不在本轮实施。

依赖关系：

```text
WP0 Snapshot freeze
  → WP1 Realignment audit
  → WP2 Governance decisions
  → WP3 Contract/schema patch
  → WP4 Asset disposition and interface validation
  → G0 REALIGNMENT_FREEZE
  → WP5 Response-blind representation necessity
  → WP6 Anchor feasibility
  → G1 RESPONSE_UNLOCK or BLOCKED
  → WP7 Joint response/context analysis
  → WP8 Attribution/claim gate
  → G2 REPAIR_READINESS
```

---

## 8. Work Package 0 — 当前快照冻结

### 目标

确保 realignment 不覆盖现有 repair lineage，也不把旧 Phase8 结果重新纳入当前结论。

### 输入

- 第 2.1 节列出的 source-of-truth；
- 当前 git status；
- Phase6 membership hash；
- Phase7A/7B handoff；
- Phase8 abort manifest。

### 任务

1. 建立 `realignment_parent_snapshot.yaml`；
2. 记录所有输入路径、SHA256、mtime、状态与优先级；
3. 标记历史 Phase8 目录为 `historical_risk_evidence_only`；
4. 记录当前 `response_model_run=false`；
5. 记录禁止修改的资产和允许新建的版本化目录；
6. 为 R0 分配新 `run_id`，不得复用历史 run_id。

### 输出

建议目录：`results/v6_2/scientific_engineering_realignment_v1/`

- `realignment_parent_snapshot.yaml`
- `canonical_input_hashes.tsv`
- `historical_asset_quarantine.tsv`
- `scope_and_prohibitions.yaml`

### Gate G0.0

- 所有 canonical 文件存在且 hash 已记录；
- membership hash 与当前 release 一致；
- 当前 Phase8 为 NOT RUN；
- 未写入任何旧结果目录；
- 未读取 response 数据值，只允许读取 response schema/availability 元数据。

失败时：停止，不进入 WP1。

---

## 9. Work Package 1 — Scientific & Engineering Realignment Audit

### 目标

执行 `refine_1.md` 要求的只读审计，形成一个唯一解释，不实施修复。

### 任务

1. 恢复 intended scientific architecture；
2. 重建 actual implemented architecture；
3. 生成 component-level drift matrix；
4. 对 8 个 FM 进行 representation/claim 双层解释；
5. 分解 governance、architecture、data identifiability、replication blockers；
6. 给出 corrected scientific target；
7. 给出 corrected engineering target；
8. 将所有现有资产分类；
9. 给出 corrected phase map；
10. 列出需要人为批准的决定。

### 唯一主输出

- `docs/v6_2/V6_2_1_SCIENTIFIC_ENGINEERING_REALIGNMENT_REPORT.md`

### 辅助机器输出

- `results/v6_2/scientific_engineering_realignment_v1/drift_matrix.tsv`
- `results/v6_2/scientific_engineering_realignment_v1/fm_role_reinterpretation.tsv`
- `results/v6_2/scientific_engineering_realignment_v1/blocker_matrix.tsv`
- `results/v6_2/scientific_engineering_realignment_v1/asset_disposition.tsv`

### 验收

- Executive verdict 不超过 10 句话；
- 明确 `do not continue FM07-only Phase8`；
- 未将旧 Phase8 effect 当成当前结果；
- 未将 identifiability audit 废弃；
- 未改变任何阈值或 FM membership；
- 每个结论可追溯到 plan、roadmap、blocker report 或机器 handoff。

---

## 10. Work Package 2 — 需要冻结的治理决定

### 目标

把不能由执行 agent 自行决定的科学解释形成显式 decision records。

### 必须批准的决定

| decision_id | 决定 | 推荐选项 | 未批准时状态 |
|---|---|---|---|
| D01 | representation set 与 claim-eligible set 是否分离 | **分离** | 保持 BLOCKED |
| D02 | identifiability gate 的职责 | **改为 attribution/claim gate；不作 representation hard filter** | 保持 BLOCKED |
| D03 | dominant barrier 允许的结果 | **single/co-dominant/coarse/diffuse/unresolved 全部允许** | 不得训练 transport |
| D04 | 原三-FM Phase8 contract | **废止并保留历史；不直接改成 FM07-only response contract** | Phase8 NOT RUN |
| D05 | response-blind Stage-A 顺序 | **置于 joint response/context 之前** | 不得继续 Stage-A/response |
| D06 | completeness 规则 | 在 response 前冻结单一规则与 conditional zone | 不得发布新 representation contract |
| D07 | response lock | **R0 和 WP5/6 全程保持锁定** | 任一读取即本轮失败 |
| D08 | SRB 地位 | 通过 necessity gate 才采用；否则 direct module backbone | 不允许默认保留 SRB |

### 输出

- `docs/v6_2/V6_2_1_REALIGNMENT_DECISION_RECORD.md`
- `results/v6_2/scientific_engineering_realignment_v1/governance_decisions.yaml`

每项记录：决策人、日期、选项、理由、被替代合同、适用 run_id 和是否可逆。

### Gate G0.1

所有 D01–D08 均有明确值；不得使用“默认同意”或由执行 agent 猜测。

---

## 11. Work Package 3 — 合同与 schema patch

### 目标

只修订职责、输入/输出和 phase 顺序；不重算数值。

### 需要新建或版本化的合同

1. `representation_set_contract_v1.yaml`
2. `measurement_reliability_contract_v1.yaml`
3. `claim_attribution_contract_v1.yaml`
4. `barrier_state_verdict_schema_v1.yaml`
5. `response_lock_contract_v1.yaml`
6. `anchor_feasibility_contract_v1.yaml`
7. `representation_necessity_eval_contract_v1.yaml`
8. `repair_readiness_contract_v1.yaml`
9. `realigned_phase_map_v1.yaml`

### representation contract 最小字段

- `module_id`
- `scoreable`
- `measurement_status`
- `reliability_weight_policy`
- `uncertainty_mask_policy`
- `coarse_group_id`
- `representation_candidate`
- `independent_claim_candidate`
- `response_used_in_assignment=false`
- `source_hash`

### barrier verdict schema

只允许：

- `single_dominant`
- `co_dominant`
- `coarse_group`
- `diffuse_mixed`
- `unresolved_abstain`

不得将“没有单一 dominant FM”编码为 pipeline failure。

### 旧合同处理

- 旧 Phase8 三-FM contract：`SUPERSEDED_BY_REALIGNMENT_V1`；
- 不删除历史文件；
- 新 manifest 必须记录 `supersedes` 和 parent hash；
- 旧 Phase7 数值继续有效，但 `allowed_for_phase8_primary` 被解释为旧 claim/input 混合字段，不直接沿用为新 representation eligibility。

### 验收

- schema validator 全部 PASS；
- representation 与 claim 字段不可互相推导；
- response lock 可检测 response/outcome/split 字段泄漏；
- contracts 不包含新生物结果；
- 没有文件覆盖旧 release。

---

## 12. Work Package 4 — 资产处置与接口验证

### 目标

最大限度复用可靠资产，只返工被漂移影响的接口。

### KEEP

- data interface 与身份修复；
- Phase4A coarse/mid annotation 及 reliability；
- Phase4B fraction-only repair 主交接；
- Phase5 confounding/negative-control 结果；
- Phase6 frozen membership 与算法比较；
- Phase6 fixed-membership projection；
- Phase7 measurement 数值、variance decomposition 和 scoring contract；
- Phase7 identifiability statistics、coarse merge map 和 strict audit。

### KEEP_BUT_REINTERPRET

- FM01/FM04 `sensitivity_only`：measurement/claim 限制，不等于不能进入联合 representation；
- FM03 `blocked_nonidentifiable`：禁止独立归因，允许与 FM02 作为 CB02 block；
- FM05/FM08：允许作为 CB04 block；
- FM07 `eligible_independent`：只表示具备未来独立 claim 审计前置资格；
- `eligible_barrier_set_v0`：保留为旧合同产物，不作为新 representation set。

### CONTRACT_PATCH_REQUIRED

- Phase7A→7B handoff 的 representation/claim 字段；
- Phase7B→Phase8 handoff；
- Phase8 expected primary set；
- completeness 规则；
- response unlock 和 repair readiness。

### LOGIC_REORDER_REQUIRED

- response-blind Stage-A 移至 joint response/context 之前；
- identifiability 拆成 pre-model measurement/coarse grouping 与 post-model attribution；
- dominant barrier 从前置单一 feature 恢复为模型后可 abstain 的 verdict。

### FUTURE_RERUN_REQUIRED

- response-blind representation necessity comparison；
- frozen-representation anchor/context analysis；
- post-model attribution；
- 仅在新证据改变 geometry 时重跑相关 identifiability；
- repair readiness 通过后才运行 transport。

### INVALID/OBSOLETE

- 旧 Phase8 response effect/route；
- 任何把旧三-FM结果写成当前 release 的摘要；
- 任何把 FM07 写成已确认 failure barrier 的资产；
- 任何把 registry/sample rows 当作独立 response patients 的口径。

### 输出

- `asset_disposition.tsv`
- `canonical_interface_map.tsv`
- `stale_reference_audit.tsv`
- `rerun_dependency_graph.yaml`

### 验收

- 每个 retained asset 有唯一 canonical path 和 hash；
- 所有 obsolete 引用从新主入口移除；
- 没有删除 provenance-rich 历史资产；
- 没有把未重建 Phase4B activity families升级为当前主输入。

---

## 13. Gate G0 — Realignment freeze

R0 只有同时满足以下条件才完成：

1. realignment report 通过 independent scientific review；
2. D01–D08 获得显式批准；
3. 所有合同/schema 通过 validator；
4. asset disposition 无悬空 canonical 引用；
5. frozen membership hash 未变化；
6. response data 未被读取用于选择或设计；
7. 当前 scientific verdict 仍是 Phase8 NOT RUN；
8. 生成 `REALIGNMENT_FREEZE_APPROVED` manifest。

输出：`realignment_freeze_manifest.yaml`。

若任一项失败，状态为 `REALIGNMENT_INCOMPLETE_RESPONSE_LOCKED`。

---

## 14. Work Package 5 — Response-blind representation necessity

状态：**FUTURE / G0 后才能运行**。

### 科学问题

联合 patient-state representation 是否比直接 frozen module vector 提供稳定、可解释且低泄漏的增益？

### 输入

- frozen module membership；
- frozen patient-timepoint module score matrix；
- representation set contract；
- reliability weights、missingness masks 和 coarse group map；
- cohort/cancer/platform/QC 仅用于分层、审计或 adversarial control；
- **禁止 response label、outcome、split-derived response information**。

### 强制比较路线

1. Direct module backbone：原始/标准化 module vector + masks；
2. Linear latent：PCA/NMF 或可回投线性 factor；
3. Hierarchical module model：处理 cohort/context 与不平衡覆盖；
4. SRB Stage-A candidate：module-factorized encoder + reconstruction；
5. 可选 graph regularization，只能作为 interaction proxy。

### split 与选择规则

- patient-grouped；
- cohort-held-out；
- cancer-held-out sensitivity；
- 所有模型使用相同 frozen folds；
- 不用 response AUC 选择 representation；
- 不用下游 response direction 调参；
- 失败 seed 全部登记。

### 主要指标

- held-out reconstruction；
- module-wise rank preservation；
- module/block interpretability；
- seed/bootstrap stability；
- cohort/cancer/platform leakage；
- missingness/coverage sensitivity；
- parameter-to-effective-unit ratio；
- coarse-group preservation；
- abstention/uncertainty availability。

### SRB adoption gate

SRB Stage-A 只有在至少一项预注册价值上稳定优于 direct backbone，且不增加 leakage、不损失 module 可解释性时采用。否则正式选择 direct module 或 linear/hierarchical route，SRB 降级为 optional implementation。

### 输出

- `representation_baseline_comparison.csv`
- `representation_leakage_audit.csv`
- `representation_interpretability_audit.csv`
- `representation_stability_audit.csv`
- `representation_route_verdict.yaml`
- `frozen_representation_manifest.yaml`

### 停止条件

- response 泄漏；
- 同一 patient 跨 split；
- membership 被修改；
- SRB 只在 training reconstruction 上更好；
- latent 主要编码 cohort/cancer/platform；
- 结果无法回投到 FM/block。

---

## 15. Work Package 6 — Anchor feasibility gate

状态：**FUTURE / 可与 WP5 并行，但仍 response-locked**。

### 科学问题

现有或新增数据是否在设计上足以区分 response effect、cancer、drug、endpoint、cohort 和 resolution？

### 本阶段允许读取

- response 是否存在；
- endpoint 类型；
- R/NR 计数；
- treatment、timepoint、cancer、cohort、resolution；
- expression/representation availability；
- patient duplication 与 paired structure。

本阶段不估计 FM/latent 与 response 的关联。

### 结构审计

1. 冻结独立统计单位；
2. 区分 sample row 与 unique patient；
3. 分开 PD-1、PD-L1、dual ICI；
4. 分开 baseline、on/post-treatment；
5. 分开 RECIST、mRECIST、pathologic response；
6. 分开 direct-mid、coarse 和 projection bridge；
7. 识别 cohort=cancer=drug=endpoint 的完全耦合；
8. 为 shared、context-specific、HCC R2 分别判断 estimability；
9. 预注册最小精度/样本要求，不以现有样本反推阈值；
10. 给出 data gap 和 acquisition priority。

### Gate G1

可能输出：

- `RESPONSE_UNLOCK_SHARED_CONTEXT`
- `RESPONSE_UNLOCK_CONTEXT_SPECIFIC_ONLY`
- `RESPONSE_UNLOCK_HCC_R2_ONLY`
- `ANCHOR_FEASIBILITY_BLOCKED`

shared/context 解锁至少需要：

- 关键 endpoint 内存在独立可比较环境；
- 每个进入估计的环境有真实 R 与 NR；
- treatment cleanliness 明确；
- patient unit 与 split 合法；
- representation 可在相同合同下生成；
- 设计矩阵不存在无法拆分的完全混杂；
- 预注册精度门通过。

若不满足，不能用复杂模型、IPW 或 prior 代替设计信息。

### 输出

- `anchor_feasibility_matrix.csv`
- `estimand_registry.yaml`
- `endpoint_exchangeability_table.csv`
- `design_rank_and_confounding_audit.csv`
- `response_unlock_verdict.yaml`
- `data_gap_priority.tsv`

---

## 16. Work Package 7 — Joint response/context analysis

状态：**FUTURE / 仅 G1 解锁后运行**。

### 科学问题

冻结的 multi-module patient state 是否与 anti-PD1/PD-L1 response/failure 有稳定关系；哪些部分共享，哪些由 context 改写？

### 约束

- representation 在读取 response 前冻结；
- 不删除共线 FM，只允许 group/block parameterization；
- within-cohort/endpoint effect 为主；
- pooled cross-endpoint 只能 sensitivity；
- 同一 patient 只贡献符合 estimand 的统计单位；
- Phase5 metadata/QC negative controls 必须同场比较；
- HCC 不得凭细胞数支配 pan-cancer effect；
- 允许无 shared effect、context-specific、mixed 或 abstain。

### 强制路线

1. direct module/block baseline；
2. frozen representation joint model；
3. hierarchical/mixed-effect cross-check；
4. leave-one-cohort/cancer/endpoint；
5. cohort/endpoint-block permutation；
6. metadata-only and QC-only controls；
7. calibration 与 abstention。

### 输出

- `joint_state_response_association_by_environment.csv`
- `shared_context_effect_decomposition.csv`
- `responder_compatible_direction.csv`
- `hcc_context_residuals.csv`
- `negative_control_comparison.csv`
- `joint_response_context_verdict.yaml`

本阶段输出的是 state-level association 与 direction，不直接产生 independent barrier claim。

---

## 17. Work Package 8 — Post-model attribution/claim gate

状态：**FUTURE / WP7 后运行**。

### 目标

将“联合状态与 response 相关”与“某个 FM 可独立归因”分开。

### 输入

- frozen representation；
- WP7 effect/decomposition；
- Phase7 measurement 与 identifiability statistics；
- coarse groups；
- bootstrap/posterior covariance；
- negative-control 结果；
- replication 结果。

### 审计

- conditional attribution stability；
- leave-one-module/block stability；
- assignment flip；
- collinearity/posterior dependence；
- cross-resolution consistency；
- carrier-state consistency；
- environment replication；
- claim sensitivity to representation route。

### 合法 verdict

- single independent barrier；
- multiple co-dominant barriers；
- coarse barrier group；
- diffuse/mixed barrier state；
- unresolved/abstain。

### 输出

- `claim_attribution_table.csv`
- `co_dominance_and_group_audit.csv`
- `barrier_state_verdict.yaml`
- `claim_boundary_by_module_or_group.csv`
- `phase_to_repair_handoff.yaml`

### Gate G2 — Repair readiness

只有具有合法 response direction、replication、claim boundary 和 target/mappability 接口的 independent barrier 或 coarse group 才能进入 repair。diffuse/mixed/unresolved 可以形成科学结果，但不得强行进入单-barrier transport。

---

## 18. 后续 repair、spatial、perturbation 与 X-class 边界

这些工作不属于 realignment immediate scope。

### 可以提前准备但不能形成机制 claim

- spatial module mapping；
- perturbation coverage inventory；
- X-class descriptor schema；
- external bulk projection feasibility；
- evidence-card templates。

### 必须等待 G2

- dominant/co-dominant barrier repair；
- responder-direction transport；
- X-conditioned displacement；
- barrier-specific X-class alignment；
- strong spatial failure niche claim；
- candidate evidence card 升级。

任何提前准备的资产必须标为 `preparatory_only_not_barrier_conditioned`。

---

## 19. 工程实现边界

### 19.1 目录与 lineage

- 当前 repair release 只读；
- 所有新资产写入 versioned `scientific_engineering_realignment_v1/`；
- 每个未来 run 使用新 `run_id`；
- manifest 必须记录 parent release、input hashes、contract version、code version；
- 不覆盖旧 Phase7/8 文件；
- 不删除历史 fail/blocked 证据。

### 19.2 Validator 优先

在任何 rerun 前先实现：

- contract schema validator；
- representation/claim separation validator；
- response leakage validator；
- patient split validator；
- endpoint pooling validator；
- historical-result leakage validator；
- membership immutability validator；
- phase dependency validator。

### 19.3 测试层级

1. schema unit tests；
2. fixture-based contract tests；
3. mutation tests：故意注入 response、同 patient 跨 split、旧 Phase8 effect、错误 endpoint pooling；
4. deterministic dry-run；
5. read-only hash audit；
6. independent reviewer；
7. 获批后才运行正式 computation。

---

## 20. Agent 分工与 ownership

| 角色 | ownership | 禁止事项 | 交付 |
|---|---|---|---|
| Controller | source-of-truth、run_id、freeze、最终 verdict | 不自行替用户作 D01–D08 决定 | manifests、gate verdict |
| Scientific architect | intended/actual architecture、drift、corrected phase map | 不编辑数据、不运行模型 | realignment report |
| Governance/contract agent | schemas、contracts、decision records | 不改 frozen assets | contract pack |
| Asset/interface agent | disposition、canonical path、hash、stale refs | 不移动 payload data | interface maps |
| Representation agent | WP5 各路线统一实现 | 不读 response | comparison artifacts |
| Anchor/statistics agent | WP6/7 estimand、exchangeability、joint model | 不改 representation | feasibility/response outputs |
| Attribution agent | WP8 claim gate | 不把 attribution 当 feature selection | barrier verdict |
| Tester | validators、fixtures、mutation、determinism | 不改 scientific verdict | verification report |
| Reviewer | correctness、leakage、claim boundary、false-green | 不替执行 agent补结果 | signed review |

写入型任务必须按目录划分 ownership；同一文件只允许一个 owner。

---

## 21. 风险与停止规则

| 风险 | 早期信号 | 强制处置 |
|---|---|---|
| 用 realignment 为 FM07-only 找理由 | 先写 FM07 response 模型再补合同 | 立即停止，保持 NOT RUN |
| response 泄漏到 representation | feature/selection 含 response、outcome、split proxy | run invalid |
| 取消质量控制 | 因多模块目标而忽略 measurement uncertainty | run invalid |
| SRB 仪式化 | SRB 未优于 direct route 仍被保留 | 降级 direct route |
| 共线被误当无信息 | 相关 FM 被删除而非 block/group | 回退 contract |
| endpoint 假同质 | mRECIST/pathologic/RECIST 直接池化 | 主分析 invalid |
| sample rows 冒充患者 | 重复 timepoint 被计为独立患者 | 统计结果 invalid |
| 旧 Phase8 泄漏 | 旧 effect/route 进入当前 verdict | release invalid |
| threshold chasing | 根据 0.492 或 response 改阈值 | realignment fail |
| 新数据只增加同一混杂 | 患者增加但 design rank 不变 | 不解锁 response |

---

## 22. 里程碑与退出标准

| milestone | 完成条件 | 允许进入 |
|---|---|---|
| M0 Snapshot frozen | hashes、prohibitions、lineage 完整 | WP1 |
| M1 Audit complete | realignment report + 4 machine tables + review PASS | WP2 |
| M2 Decisions frozen | D01–D08 全部显式批准 | WP3 |
| M3 Contracts valid | schemas、validators、mutation tests PASS | WP4 |
| M4 Realignment frozen | asset map clean；response lock；freeze manifest | WP5/6 |
| M5 Representation frozen | necessity gate 选定 direct/linear/SRB route | G1 |
| M6 Anchor feasible | response unlock verdict 非 blocked | WP7 |
| M7 Joint direction valid | shared/context verdict + controls + abstention | WP8 |
| M8 Claim gate complete | independent/coarse/mixed/unresolved verdict | G2 |
| M9 Repair ready | 至少一个合法 repair target/group | transport |

任一 milestone 允许输出 `BLOCKED` 或 `UNRESOLVED`；不得为推进 phase 而降低 gate。

---

## 23. 本次只应立即执行的范围

在用户批准本方案后，本轮建议只执行 R0：

1. WP0 当前快照冻结；
2. WP1 生成 `V6_2_1_SCIENTIFIC_ENGINEERING_REALIGNMENT_REPORT.md`；
3. WP2 提交 D01–D08 决策包；
4. 在决策获批后执行 WP3/4；
5. 生成 `REALIGNMENT_FREEZE_APPROVED` 或 `REALIGNMENT_INCOMPLETE_RESPONSE_LOCKED`。

本轮明确不执行 WP5–8，不运行 response model，不训练 SRB，不重新筛选 FM，不下载数据。

---

## 24. 最终验收清单

- [ ] 当前 release 与 parent hashes 已冻结；
- [ ] 旧 Phase8 已明确降为 historical risk evidence；
- [ ] intended 与 actual architecture 均有可审计描述；
- [ ] drift matrix 覆盖 refine_1 要求的 10 个 component；
- [ ] 8 个 FM 均有 representation 与 claim 双层解释；
- [ ] 四类 blockers 独立记录；
- [ ] D01–D08 获得显式批准；
- [ ] representation、claim、response unlock、repair readiness 合同分离；
- [ ] completeness 规则在 response 前冻结；
- [ ] old three-FM Phase8 contract 已版本化 supersede；
- [ ] frozen membership hash 未改变；
- [ ] response lock 通过；
- [ ] patient、endpoint、historical leakage mutation tests 通过；
- [ ] 无任何新 response、repair、spatial 或 X-class claim；
- [ ] independent reviewer 签署 PASS 或明确 BLOCKED；
- [ ] 所有未执行工作保持 `NOT_RUN`。

---

## 25. 最终执行判断

项目不应从当前状态直接进入 FM07-only response analysis。下一真正阶段应命名为：

> **Pre-anchor architecture realignment and response-blind representation freeze**

该阶段的目标不是选出最显著 FM，而是恢复 multi-module patient-state 作为统计对象，冻结 representation/claim 的职责边界，并证明 direct module、linear latent 或 SRB Stage-A 中哪一种 response-blind 表征值得进入未来 anchor/context 分析。

科学解锁仍依赖可交换 response anchors；工程解锁依赖 realignment contract。二者必须分别过门，不能互相替代。
# 执行模式说明（2026-08-27）

本文件记录的是早期 realignment 方案。当前项目已切换为探索性生信研究；其中预注册式 scope、冻结 estimand、response unlock 和固定 gate 不再作为新分析的硬性约束。历史版本保留用于溯源，后续执行以 `EXPLORATORY_MODE_DIRECTIVE_20260827.md` 为准。
