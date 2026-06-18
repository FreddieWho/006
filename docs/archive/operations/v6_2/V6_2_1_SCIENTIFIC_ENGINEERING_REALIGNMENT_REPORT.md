# v6.2.1 Scientific & Engineering Realignment Audit

状态：**G0 FREEZE APPROVED / WP5 COMPLETE / WP6a COMPLETE-METADATA-ONLY / G1 BLOCKED / RESPONSE LOCKED**
run_id：`v6_2_realignment_v1_20260820T212638Z`
父 release：`context_repair_v1_full_rerun`
父 release 状态：`COMPLETE_REPAIRS_PHASE8_GATE_BLOCKED`

本文是 R0 的只读审计结果，不是新的生物分析结果。它只校准路线、合同和证据职责，不重算 Phase6/7 数值，不启动 Phase8，不读取 response/outcome 数值。

## Executive verdict

1. v6.2 原本要建模的是多模块的 patient immune-tumor state，而不是单个已通过 identifiability 的 FM。
2. 当前实现把 claim-level 的 identifiability 结果过早当成 representation-level 的 feature filter，形成了 high-severity、但可恢复的架构漂移。
3. Phase7 的 measurement、identifiability、coarse merge 和 strict audit 数值本身保留有效；需要纠正的是它们在 pipeline 中承担的角色。
4. 当前 Phase8 的三-FM contract 与 Phase7 实际 FM07-only primary set 不一致，因此 Phase8 的正确状态仍是 **NOT RUN / gate blocked**。
5. 不应为了消除工程阻塞而直接把输入改成 FM07-only；这会把架构漂移固定为单特征 response 问题。
6. FM07 目前只是进入未来独立归因审计的候选，不是已确认的 response/failure/resistance barrier。
7. 下一步应把 scoreable、representation 和 claim-eligible 三个集合拆开，并在 response-blind 条件下比较 direct module、linear latent 与 SRB Stage-A。
8. response/context、repair、空间验证、perturbation 和 X-class 外推都必须等待 realignment freeze 及后续专门 gate。

## R1 执行结果（截至 2026-08-21）

- G0 已完成：D01–D08 均获导师用户批准，9 个合同/模式文件解除 draft-only，Phase8 仍明确为 `NOT RUN`，response lock 仍有效。
- WP5 已完成且全程 response-blind。透明的 `direct_module` 保留为主干；线性 latent、hierarchical module 和线性 SRB Stage-A 候选均未通过预注册的 SRB necessity gate。该结果只支持“表征路线选择”，不支持任何疗效方向或 barrier 结论。
- WP6a 已完成但只读队列/设计元数据：68 个 cohort registry 行与 68 个环境元数据行一一对应，未读取 response/outcome 数值，也未运行数值对比。元数据显示治疗情境、癌种、endpoint 和 timepoint 并非天然可交换，因此当前不能定义一个无需额外约束的共享 response 方向。
- 当前阻塞转为 G1：需先选择明确 estimand，并单独批准 WP6b 的 response 解锁。未完成前不得读取数值 response，不得启动 joint response/context model。

详细证据：`results/v6_2/scientific_engineering_realignment_v1/realignment_freeze_manifest.yaml`、`wp5_response_blind_representation/representation_route_verdict.yaml`、`wp6a_anchor_metadata_audit/anchor_metadata_summary.yaml`、`wp6a_anchor_metadata_audit/wp6a_blocker_record.yaml`。

## 1. 审计边界与输入

### 1.1 读取的控制面

- intended architecture：`docs/plan/v6_2/plan.md`、`docs/plan/v6_2/roadmap.md`；
- 当前状态：`results/v6_2/V6_2_PROJECT_SCIENTIFIC_LOGIC_AND_BLOCKER_REPORT_20260731.md`；
- release 与 handoff：context repair、Phase6 entry、Phase7A/7B handoff、Phase8 input/run manifest；
- 已生成的 `refine_1.md` 作为本次审计的直接约束。

### 1.2 不可改写的事实

- Phase6 frozen membership SHA256：`285c11de4b071c45b0ec43464fb963650340c29bb2de3037570578b29ee44036`；
- 8 个 FM 均来自 response-blind 流程；
- Phase7 measurement 与 identifiability 结果保留；
- 当前只有 FM07 具备进入未来 independent-claim audit 的前置资格；
- Phase8：`ABORTED_PRE_MODEL` / `response_model_run=false`；
- 没有当前 response direction、failure barrier、repair transport 或 X-class 结论。

完整 hash、mtime、父 release 和写入禁区见 `results/v6_2/scientific_engineering_realignment_v1/realignment_parent_snapshot.yaml`。

## 2. Intended architecture 与 actual implementation

### 2.1 Intended architecture

`plan.md` 定义的主线是：固定模块入口 → response-blind measurement → response-blind SRB Stage-A latent → anchor/context head → dominant-barrier attribution → repair transport → spatial/perturbation/external validation → X-class alignment。

其中，FM 是 patient state 的可解释坐标；identifiability gate 的作用是保护独立归因、允许 coarse merge，并在无法辨识时 abstain，而不是先把所有相关坐标从联合 representation 中删除。

`roadmap.md` 同时要求先保证数据、模块和评测可靠，再训练 SRB；但其 Phase7/Phase8 排序与 `plan.md` 中“Barrier Identifiability 输入包含 SRB Stage-A latent”的设计存在内部张力。这是需要在 realignment 中显式解决的路线问题。

### 2.2 Actual implementation

当前实际链路可概括为：

```text
8 个 response-blind FM
  → Phase7A measurement
  → Phase7B identifiability
  → eligible_barrier_set_v0
  → Phase8 contract 期望 FM01+FM04+FM07
  → 实际 primary 只剩 FM07
  → Phase8 ABORTED_PRE_MODEL
```

这条链路同时混合了三件不同的事情：能否被稳定测量、能否作为联合 representation、能否在 response/context 分析后支持独立机制 claim。混合的直接技术表现是旧 Phase8 contract 与实际 primary set 不一致。

## 3. Drift matrix 摘要

详细矩阵见 `results/v6_2/scientific_engineering_realignment_v1/drift_matrix.tsv`。

| 部件 | intended | actual | 判定 |
|---|---|---|---|
| 统计对象 | 多模块 patient state | eligible barrier set 作为唯一输入 | 高严重度架构漂移 |
| measurement gate | 给出可靠性、完整性和不确定性 | 与 downstream eligibility 混用 | 治理漂移 |
| identifiability gate | 保护独立归因并支持 coarse merge | 过早变成 feature hard filter | 高严重度架构漂移 |
| SRB Stage-A | response-blind reconstruction 后再接监督 | 尚未比较，Phase8 在 pre-model 中止 | 阶段顺序未完成 |
| Phase8 contract | 与已批准 representation/claim 语义一致 | 期望 FM01+FM04+FM07，实际 FM07 | governance hard fail |
| dominant barrier | 模型后 verdict，可 abstain/合并 | 容易被误写成 eligible FM | claim boundary 风险 |

## 4. 八个 FM 的双层解释

详细表见 `fm_role_reinterpretation.tsv`。当前解释只使用 response-blind measurement 与 identifiability 证据。

- FM01、FM04：measurement 可用但属于 sensitivity-only；不应独立进入 response claim，也不应因此被排除出联合 representation。
- FM02、FM05、FM08：存在 CB02/CB04 相关性，适合保留为联合坐标或 coarse block，不适合直接作独立 barrier。
- FM03：measurement 指标尚可，但 Phase7B 为 `blocked_nonidentifiable`，应进入 CB02 coarse representation，禁止独立归因。
- FM06：measurement completeness 更低，representation 中必须带 missingness/uncertainty mask；不应强行解释为阴性生物学。
- FM07：当前唯一 `eligible_independent` 坐标，适合保留为联合 representation 中的一维，并作为未来独立 attribution audit 的候选；当前仍是 response-blind program。

这些状态不是 response 结果，也不能用来命名 failure、resistance 或 treatment-response barrier。

## 5. 三个集合的校准

### 5.1 Scoreable set

满足输入可计算、来源可追溯、measurement contract 可审计的 FM。它回答“这个坐标能否被稳定算出”。

### 5.2 Representation set

在 scoreable 基础上，允许带可靠性权重、缺失掩码和 coarse-group 标记的 response-blind 联合表示。它回答“这个坐标是否值得作为 patient-state 的信息保留”。进入 representation 不等于拥有独立 mechanism claim。

### 5.3 Claim-eligible set

在未来冻结 representation、明确 response estimand、完成 anchor/context 分析后，才判断某一 FM 或 coarse group 是否能支持独立 attribution、dominant-barrier 或 repair claim。它回答“这个坐标能否承担强科学主张”。

校准后，`eligible_barrier_set_v0.csv` 只能作为旧 claim/input 混合合同的审计产物，不能直接作为新 representation set。

## 6. Blocker 分层

详细矩阵见 `blocker_matrix.tsv`。

### 6.1 直接工程阻塞

Phase8 三-FM contract 与实际 FM07-only primary set 不一致。修复方式不是默默删掉两个 FM，而是版本化旧 contract、记录 parent hash，并重新定义 representation 与 claim 字段。

### 6.2 架构阻塞

如果仍把 identifiability 当作 representation hard filter，就无法回答“多个共线模块是否共同组成 patient state、context 如何调制该状态、SRB 是否提供额外信息”这些原始科学问题。

### 6.3 科学识别阻塞

现有 anchors 在 cancer、drug、endpoint 和 cohort 上同时变化，无法直接支持一个可交换的 shared anti-PD-1 direction。即使 Phase8 contract 修好，也不能把 response 方向自动解释成跨环境的机制 barrier。

### 6.4 上游依赖阻塞

该问题牵涉 Phase4B 的 patient-timepoint aggregation、Phase5 的 confounding/contrast 审计、Phase6 frozen module entry、Phase7 measurement/identifiability 以及 anchor metadata。需要调整的是接口与 estimand，不是静默重跑或降低阈值。

## 7. Corrected scientific target

当前应回答的科学问题按顺序改写为：

1. 在不使用 response 的前提下，8 个 FM 中哪些坐标可可靠测量、哪些只能合并为 coarse group？
2. 在相同 frozen membership 上，direct module、linear latent 和 SRB Stage-A 哪种表征最能保留 patient-timepoint state，并且不牺牲可解释性？
3. 在 representation 冻结后，哪些 response anchors 具有可交换的 response direction；哪些只能作为 sensitivity/support-only？
4. 在通过 attribution gate 后，dominant barrier 的结果是 single、co-dominant、coarse、diffuse/mixed，还是 unresolved/abstain？
5. 只有上述问题回答后，才有资格讨论 repair、空间 niche、perturbation 或 X-class。

## 8. Corrected engineering target

工程上需要完成的不是“让旧 Phase8 运行”，而是恢复以下合同链：

```text
frozen parent snapshot
  → representation/claim split
  → response-blind necessity gate
  → anchor feasibility metadata gate (WP6a)
  → separately unlocked response-contrast audit (WP6b)
  → joint response/context model
  → post-model attribution/claim gate
  → repair readiness
```

每个合同都必须记录 source hash、response lock、适用 run_id、supersedes parent 和可逆性。任何“没有 single dominant FM”的结果都必须能编码为 `co_dominant`、`coarse_group`、`diffuse_mixed` 或 `unresolved_abstain`，不能被编码成 pipeline failure。

## 9. 资产处置

完整清单见 `asset_disposition.tsv`；历史 Phase8 表面见 `historical_asset_quarantine.tsv`。

- KEEP：data/interface repair、Phase4A/4B、Phase5 审计、Phase6 frozen membership、Phase7 measurement 和 identifiability 数值。
- KEEP_BUT_REINTERPRET：FM01/FM04 sensitivity-only、FM03 blocked_nonidentifiable、FM05/FM08 coarse block、FM07 eligible_independent、`eligible_barrier_set_v0`。
- CONTRACT_PATCH_REQUIRED：Phase7A→7B、Phase7B→Phase8、Phase8 expected set、completeness、response unlock 和 repair readiness。
- QUARANTINE：旧 Phase8 effect、route、sample-surface、bounded/readjudication 和 context-adjudication 目录；它们只用于追溯漂移。

## 10. 需要导师明确批准的决定

以下决定不能由执行 agent 默认代替。当前全部记录为 `PENDING_HUMAN_APPROVAL`，未批准前保持 `REALIGNMENT_INCOMPLETE_RESPONSE_LOCKED`：

| ID | 建议冻结内容 | 不批准的后果 |
|---|---|---|
| D01 | representation set 与 claim-eligible set 分离 | 保持 blocked |
| D02 | identifiability 改为 attribution/claim gate，不作 representation hard filter | 保持 blocked |
| D03 | dominant barrier 允许 single/co-dominant/coarse/diffuse/unresolved | 不得训练 repair transport |
| D04 | 废止旧三-FM Phase8 contract，不直接改成 FM07-only response contract | Phase8 NOT RUN |
| D05 | response-blind Stage-A 置于 joint response/context 之前 | 不得继续 response 分析 |
| D06 | response 前冻结唯一 completeness 规则及 conditional zone | 不得发布 representation contract |
| D07 | R0、WP5、WP6a 全程保持 response lock；WP6b 单独解锁 | 本轮失败并停止 |
| D08 | SRB 通过 necessity gate 后才采用，否则保留 direct-module backbone | 不得默认保留 SRB |

决策包见 `docs/v6_2/V6_2_1_REALIGNMENT_DECISION_RECORD.md` 与 `results/v6_2/scientific_engineering_realignment_v1/governance_decisions.yaml`。

## 11. 当前停止点

WP0 快照、WP1 审计、WP2 决策包准备、WP3 合同草案和 WP4 资产处置均可在 response lock 下完成。未来 WP6a 仅做 metadata-only anchor feasibility；WP6b 的 response contrast 数值审计必须另行解锁。R0 的最终 `REALIGNMENT_FREEZE_APPROVED` 不能由执行 agent 自行产生，必须等 D01–D08 的明确批准及独立 scientific review。

在此之前，禁止：启动 Phase8 response model、读取 response/outcome 数值、修改 frozen membership、降低 completeness/identifiability 阈值、把旧 Phase8 effect 当成当前结果，或将 FM07 写成已确认的 failure/resistance barrier。WP6b 仍未解锁。

审计依据包括：`refine_1.md`、`docs/plan/v6_2/plan.md`、`docs/plan/v6_2/roadmap.md`、当前 repair release/handoff、Phase7 measurement/identifiability tables 及 Phase8 abort manifest。
