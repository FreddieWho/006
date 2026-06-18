# v6.2.1 Realignment Decision Record

状态：**APPROVED**
run_id：`v6_2_realignment_v1_20260820T212638Z`
适用范围：R0 realignment control package

批准人：导师用户
批准时间：2026-08-21 23:23:34 CST
批准依据：用户明确同意本记录中的推荐选项

本文件只记录需要导师明确批准的研究设计与治理决定。`recommended_option` 是审计后的建议，不是已生效的配置；任何一项未批准，R0 都保持 `REALIGNMENT_INCOMPLETE_RESPONSE_LOCKED`。

## 决策总览

| ID | 建议选项 | 为什么需要明确批准 | 未批准时 |
|---|---|---|---|
| D01 | representation set 与 claim-eligible set 分离 | 改变下游输入资格的定义 | BLOCKED |
| D02 | identifiability 改为 attribution/claim gate，不作 representation hard filter | 改变 Phase7→SRB/Phase8 的职责 | BLOCKED |
| D03 | dominant barrier 允许 single、co-dominant、coarse、diffuse/mixed、unresolved/abstain | 改变“没有单一 barrier”是否算失败 | BLOCKED |
| D04 | 旧三-FM Phase8 contract 版本化废止；不直接改成 FM07-only | 改变当前 Phase8 是否可启动 | Phase8 NOT RUN |
| D05 | response-blind Stage-A 先于 joint response/context | 改变 SRB 与监督 head 的顺序 | BLOCKED |
| D06 | response 前冻结唯一 completeness 规则和 conditional zone | 解决书面 0.80 与实际 0.734 的冲突 | 不发布 representation contract |
| D07 | R0/WP5/WP6a 全程 response lock；WP6b 需独立解锁 | 防止 response 泄漏进表征和设计 | 本轮失败并停止 |
| D08 | SRB 只有通过 necessity gate 才采用，否则保留 direct-module backbone | 改变模型主干与方法学主张 | 不得默认保留 SRB |

## 逐项记录

### D01 — representation 与 claim 分离

- 当前状态：`APPROVED`
- 推荐值：`separate_sets=true`
- 证据：`refine_1.md` 的 representation/claim 分离要求；`plan.md` 的多模块 patient-state 目标；当前 Phase7B 只产生 claim/input 混合的 `eligible_barrier_set_v0`。
- 影响：所有 scoreable、measurement-reliable 或 coarse-group 坐标可先进入 response-blind representation；独立 claim 另过 attribution gate。
- 可逆性：可逆，但会影响后续 handoff schema。

### D02 — identifiability gate 职责

- 当前状态：`APPROVED`
- 推荐值：`post_model_attribution_gate=true; pre_model_hard_filter=false`
- 证据：`plan.md` 的 Barrier Identifiability 输入包含 SRB Stage-A latent；Phase7B 的共线性和 coarse merge 结果。
- 影响：保留相关信号用于联合表征，同时禁止其承担独立机制归因。
- 可逆性：可逆，但需要版本化 Phase7B→Phase8 handoff。

### D03 — dominant barrier verdict

- 当前状态：`APPROVED`
- 推荐值：`single_dominant | co_dominant | coarse_group | diffuse_mixed | unresolved_abstain`
- 证据：`refine_1.md` 明确不强制 single dominant barrier；`plan.md` 要求不可辨识模块合并 coarse class。
- 影响：无单一 barrier 时仍可形成诚实科学结果，但不得强行进入单-barrier repair transport。
- 可逆性：可逆，属于 verdict schema。

### D04 — 旧 Phase8 contract

- 当前状态：`APPROVED`
- 推荐值：`supersede_old_three_fm_contract=true; do_not_convert_to_fm07_only=true`
- 证据：当前 Phase8 contract 期望 `FM01+FM04+FM07`，Phase7B primary 只有 `FM07`，且 `response_model_run=false`。
- 影响：保留 NOT RUN 语义；新 contract 必须从 representation/claim 分层重新定义。
- 可逆性：旧文件保留，新版本可回溯 parent hash。

### D05 — response-blind Stage-A 顺序

- 当前状态：`APPROVED`
- 推荐值：`response_blind_stage_a_before_joint_response_context=true`
- 证据：`plan.md` 的 SRB Stage-A reconstruction 与 context modulation 设计；`roadmap.md` 的响应盲模块/表征原则。
- 影响：先验证表征是否有必要，再引入 response；避免稀缺监督反向定义 representation。
- 可逆性：可逆，但会改变 phase map。

### D06 — completeness 规则

- 当前状态：`APPROVED`
- 已选值：`tiered_explicit_completeness_v1`
- 证据：当前文字主线 `completeness_mainline_min=0.80` 与 FM07 实际可进入旧 eligibility 的 `0.734` 存在冲突。
- 互斥选项：
  - `strict_claim_v1`：`mainline_min=0.80`；`0.60<=x<0.80` 只能作为 scoreable/conditional support，不得进入独立 claim；`x<0.60` 仅保留审计记录。
  - `tiered_explicit_completeness_v1`：`mainline_min=0.80`；`conditional=[0.60,0.80)`，允许进入带 mask 的联合 representation 但只能 support-only；`x<0.60` 不进入 representation。
  - `legacy_reproduction_v1`：复现当前实现的 `<0.60` 阻断与 `0.60–<0.80` conditional，但必须明确标记为 legacy，不得冒充书面 mainline 规则。
- 导师必须写入 `selected_option`、`mainline_min`、`conditional_lower_inclusive`、`conditional_upper_exclusive`、`below_conditional_behavior` 和生效版本；在此之前不发布新的 representation eligibility。
- 影响：不得通过降低阈值消除阻塞。
- 可逆性：需新版本合同，不改写历史数值。

### D07 — response lock

- 当前状态：`APPROVED`
- 推荐值：`lock_R0_WP5_WP6a=true; WP6b_requires_separate_unlock=true`
- 证据：`refine_1.md` 要求 response-blind realignment；当前 Phase8 尚未运行。
- 影响：WP6a 只允许读取 schema/availability/标签 provenance 元数据；WP6b 才能在独立解锁后读取 contrast 数值。R0/WP5/WP6a 禁止读取 response/outcome 值用于选择、调参、分组或关联检验。
- 可逆性：只有通过单独 gate 才能解锁。

### D08 — SRB 地位

- 当前状态：`APPROVED`
- 推荐值：`necessity_gate_then_adopt_or_fallback_to_direct_module`
- 证据：当前尚未完成 direct module、linear latent、SRB Stage-A 的 response-blind necessity comparison。
- 影响：SRB 不能因为原计划偏好而默认保留；也不能在没有比较时宣称方法学必要性。导师还需批准最小优势/不劣界、held-out unit、bootstrap 规则和失败判据。
- 可逆性：可逆，属于模型选择 gate。

## 批准栏

| decision_id | mentor_decision | decision_date | rationale | signature_or_trace |
|---|---|---|---|---|
| D01 | APPROVED | 2026-08-21 | 用户同意推荐架构 | mentor_user |
| D02 | APPROVED | 2026-08-21 | 用户同意推荐架构 | mentor_user |
| D03 | APPROVED | 2026-08-21 | 用户同意推荐架构 | mentor_user |
| D04 | APPROVED | 2026-08-21 | 用户同意推荐架构 | mentor_user |
| D05 | APPROVED | 2026-08-21 | 用户同意推荐架构 | mentor_user |
| D06 | APPROVED | 2026-08-21 | 用户选择 tiered_explicit_completeness_v1 | mentor_user |
| D07 | APPROVED | 2026-08-21 | 用户同意推荐 response lock | mentor_user |
| D08 | APPROVED | 2026-08-21 | 用户同意 necessity gate | mentor_user |
