# scope_lock_v6_1

## 1) Project title
- English: **PD-1-centered TME Mechanism Nomination Framework with Shared/Specific Immune Modules and Spatially Anchored Validation (v6.1)**
- 中文：**v6.1：以 PD-1 为中心的 TME 机制提名框架：共享/特异免疫模块与空间锚定验证**

## 2) Main scientific question
在多癌种 ICI 数据中，先以 PD-1/ICI monotherapy 或较干净 ICI 队列学习基础免疫敏感性与原发耐药模块，再在 PD-1+X 场景解释 X 修复了哪些残余 TME 障碍，最终输出可被空间/组织证据、扰动先验、外部 bulk 锚点与最小实验方案支持的机制假设与候选靶点。

## 3) Primary analysis
- P1. Cohort registry + treatment context lock（先锁样本与治疗上下文，后做建模）
- P2. Immune-state measurement（cell fraction / pseudobulk / patient-level immune state features）
- P3. Strong baseline first（fraction/signature/pathway + elastic net/logistic/Cox/XGBoost + static module + scVI latent baseline）
- P4. Shared immune module + HCC-specific barrier module discovery（先稳定，再谈因果）
- P5. PD-1 anchor -> PD-1+X repair logic（仅解释 X 修复何种残余障碍）
- P6. Spatial/tissue biological adjudication（关键生物学裁判层，不是装饰层）
- P7. Translational nomination（L1/L2/L3 证据卡 + 外部锚点 + minimal wet-lab feasibility）

## 4) Secondary analysis
- S1. 高混杂联合治疗（TACE/HAIC/radiotherapy/multi-line）仅用于敏感性分析与叙事支持，不进入主训练标签核心。
- S2. 动态 pre/post 轨迹分析只在 timepoint 与 center 质量足够时作为 secondary，不能覆盖主结论。
- S3. MVP archetype/trajectory/axis/TME-MOA/virtual drug shift 仅允许作为 proxy 展示、baseline 或 sanity check。
- S4. CTLA-4 仅可作为 ICI 参考轴，不能替代 PD-1 主轴。

## 5) Excluded analysis
- E1. 将 PD-1 单药与 PD-1+X 拆成两个独立课题进行平行主结论。
- E2. 不经过强基线对照直接上复杂模型并给出主结论。
- E3. 仅凭 marker 或 marker 均值/简单组合推断机制。
- E4. 仅凭 LINCS reversal 独立提名药物或机制。
- E5. 将空间图作为装饰性验证而非机制裁判证据。
- E6. 将高混杂联合治疗标签作为主训练标签核心。
- E7. 将手工 MOA、heuristic counterfactual、pseudo-archetype bulk assignment、demo ranking 作为主结果。
- E8. 搭建网页作为主交付。

## 6) Treatment continuum definition
- Context A (Anchor): PD-1/PD-L1 monotherapy + 较干净 ICI 队列。用途：学习基础敏感性/耐药轴。
- Context B (Extension): PD-1+X（anti-VEGF/TKI/other immune or targeted with clear metadata）。用途：解释 X 的 repair 作用。
- Context C (Sensitivity only): 高混杂治疗场景。用途：敏感性分析、外部支持、反证；不进入主标签训练。

## 7) Evidence levels (3-tier only)
- **L1 Robust association**：跨队列/跨中心稳健、强基线可比、关键敏感性分析通过。
- **L2 Mechanistic support**：pathway/TF、perturb prior、ligand-receptor、spatial/tissue 等机制证据支持。
- **L3 Translational nomination**：外部 bulk/clinical anchor、druggability、spatial adjudication、最小实验可行性支持。

升级规则：L3 必须建立在 L1+L2 之上。任何层级若出现关键冲突（如空间反证、方向反证），自动降级并转入附录或敏感性结果。

## 8) MVP inheritance policy
- 直接继承（可用于正式主线输入）：cohort registry 思路、sidecar-first 数据组织、patient/sample metadata、pseudobulk 工作流、scVI atlas（QC/表征底座）、patient-level 表示经验、scope lock 与输出契约。
- 降级继承（只能 proxy/sanity）：archetype、trajectory、axis、TME-MOA、virtual drug shift。
- 禁止主结果继承：手工 MOA 直推机制、heuristic counterfactual、pseudo-archetype bulk assignment、demo-oriented ranking。

## 9) Main output list
- Core tables:
  - `cohort_registry_v6_1.csv`
  - `treatment_context_flags.csv`
  - `strong_baseline_results.csv`
  - `immune_module_scores.csv`
  - `shared_specific_modules.csv`
  - `spatial_adjudication_results.csv`
  - `M_signals_v6_1.csv`
  - `PD1X_nomination_cards.csv`
  - `external_validation_results.csv`
  - `wetlab_priority_list.csv`
- Core documents:
  - `scope_lock_v6_1.md`
  - `analysis_contract.md`
  - `decision_log.md`
  - `spatial_adjudication_report.md`
  - `PD1_anchor_to_X_repair_logic.md`
  - `minimal_validation_plan.md`
- Core claim package:
  - 至少 1 个 shared PD-1 sensitivity module
  - 至少 1 个 HCC-specific barrier module
  - 至少 1 个经空间/组织裁判支持的机制链
  - 至少 1 个可实验落地的候选靶点

## 10) Closure / fallback logic
- C1. 若复杂模型未超过强基线：关闭“复杂预测主叙事”，保留“稳健模块 + 机制解释”路径。
- C2. 若空间/组织不支持核心机制：机制主结论降级为 association，不进入强 PD-1+X nomination。
- C3. 若 PD-1 anchor 不稳定：停止 Step4 及之后，回滚到 cohort/label/treatment flag 重锁。
- C4. 若 perturb prior 与 patient+spatial 方向冲突：perturb 仅作反证，不用于提名推进。
- C5. 若仅有 Context C 支持而无 Context A/B 核心证据：结果只能归类 sensitivity，不得写入主结论。

## 11) Execution boundary (hard gate)
- Gate-G0：未完成 registry/context lock，不得进入任何建模步骤。
- Gate-G1：未完成 strong baseline report，不得提交复杂模型主结果。
- Gate-G2：未形成 PD-1 anchor，不得进入 PD-1+X repair nomination。
- Gate-G3：无 spatial adjudication，不得发布 L2/L3 机制主结论。
- Gate-G4：所有主结论必须可追溯到 `output_manifest` 记录与中间表链路。
