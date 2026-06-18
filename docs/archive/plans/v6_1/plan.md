# 科研方案 v6.1

**PD-1-centered TME mechanism nomination framework with shared-specific immune modules and spatially anchored validation**

中文名：

**v6.1：以 PD-1 为中心的 TME 机制提名框架：共享/特异免疫模块与空间锚定验证**

---

## 0. v6.1 相比 v6 的核心变化

v6.1 不推翻 v6，而是做三处简化和一处强化。

第一，治疗上下文不再拆成很多互相隔离的任务，而是改成 **PD-1-centered continuum**：

* PD-1 / ICI monotherapy：学习基础 ICI 敏感性与耐药轴；
* PD-1+X：解释在 PD-1 基础上还需要修复什么 TME 障碍；
* 高混杂联合治疗：作为敏感性分析或外部支持，不作为主训练标签的核心来源。

这保留了你说的“PD-1 solo 是 monotherapy effect 的学习锚点，然后再往上堆 PD-1+X 更有说服力”。我认为这比硬拆治疗层更合理。

第二，证据等级从复杂分层压缩成三档：

* **Level 1：Robust association**
* **Level 2：Mechanistic support**
* **Level 3：Translational nomination**

第三，强基线前置。v6.1 不允许高级模型先行。先用 cell fraction、curated signatures、pseudobulk pathway、elastic net / XGBoost、静态模块等强基线建立“最低可防守结果”。

第四，空间/组织数据从“增强项”提升为 **key biological adjudicator**。它不一定进主模型，但必须参与判断核心 TME 机制是否真实存在于组织生态中。

这个方向仍然继承 v5.3 的主轴：多癌种免疫治疗、共享+特异免疫机制、M-signals、ICSmini/e_X-mini、外部锚点和阶段门控。v5.3 原本就把共享+特异因果 GRN、PD-1+X / X+X 排序、体外/体内验证和机制读数作为主目标。

---

# 1. 总目标

v6.1 的总目标是：

**在多癌种 ICI 数据中学习可复现的共享/特异 TME 免疫模块，以 PD-1 单药/ICI 单药反应作为基础免疫敏感性锚点，进一步解释 PD-1+X 中 X 所修复的残余免疫障碍，并输出可被空间组织证据、扰动先验和最小实验支持的机制类 X 与候选靶点。**

注意，这里不再把输出叫“临床推荐系统”，而叫：

**PD-1+X hypothesis nomination framework**

也就是 **机制假设提名系统**。
这不是怂，是不装。没有前瞻临床验证之前，直接说“推荐系统”属于科研版开盲盒，审稿人不会感动，只会磨刀。

---

# 2. 核心科学问题

v6.1 聚焦三个问题。

## Q1：什么是跨癌种可复现的 PD-1 敏感性基础模块？

这主要用 PD-1 / ICI monotherapy 或较干净的 ICI 队列学习。

核心不是“谁 responder”，而是：

* 哪些 T cell / APC / IFN / antigen presentation / myeloid 状态与 ICI 敏感性稳定相关；
* 哪些模块跨癌种共享；
* 哪些模块在 HCC 中被肝脏免疫耐受背景改写。

## Q2：PD-1 不够时，主要残余障碍是什么？

这一步承接 PD-1 solo 学到的 monotherapy response biology。

核心问题是：

* non-responder 是因为 T cell 没来？
* 来了但功能被髓系/Treg/基质压住？
* 抗原呈递不够？
* IFN/炎症通路有但空间上不可接触？
* VEGF/血管/基质屏障阻断了免疫进入？

这就是 PD-1+X 里 X 的生物学位置。

## Q3：哪些 X 机制类最可能修复这些障碍？

v6.1 不优先输出一串具体药名，而是先输出机制类：

* myeloid reprogramming
* vascular / stromal remodeling
* antigen presentation / IFN restoration
* T cell activation / co-stimulation support
* epigenetic immune sensitization

然后才映射到候选靶点或药物家族。

---

# 3. 治疗结构：从“分层”改为“PD-1-centered continuum”

这是 v6.1 最关键的结构修改。

## 3.1 主分析层：PD-1 / ICI monotherapy anchor

这层用于学习基础 ICI 敏感性。

纳入：

* PD-1 单药；
* PD-L1 单药；
* CTLA-4 可作为辅助 ICI 参考，但不作为 PD-1 主轴；
* 较干净的 ICI 队列。

输出：

* PD-1 sensitivity module；
* primary resistance module；
* responder-like immune state；
* non-responder immune barriers。

## 3.2 扩展解释层：PD-1+X combination extension

这层不是独立重开一个问题，而是在 monotherapy anchor 上解释：

> 这个患者或这个亚型，在 PD-1 基础上还缺什么？

纳入：

* PD-1 + anti-VEGF / TKI；
* PD-1 + locoregional therapy，谨慎；
* PD-1 + other immune or targeted agents，若 metadata 清楚。

输出：

* X mechanism class；
* combination-compatible module；
* PD-1 anchor 与 X repair module 的关系。

## 3.3 敏感性分析层：high-confounding context

这层只作为支持或反证。

包括：

* TACE / HAIC / radiotherapy 强干预混合队列；
* 多线治疗后样本；
* timepoint 不清楚的 post-treatment；
* therapy regimen 过度混杂的 public cohort。

用途：

* 不作为主训练标签；
* 可用于验证模块方向；
* 可用于 HCC 转化叙事；
* 不可用于直接声称 PD-1+X 推荐。

这比之前把治疗拆成五六类更简单，也更符合你“PD-1 solo 是底座，PD-1+X 是增量解释”的逻辑。

---

# 4. 证据等级：三档，不搞迷宫

v6.1 使用三档证据。

## Level 1：Robust association

回答：

> 这个模块是否稳定关联 ICI response / resistance？

证据来源：

* 多队列复现；
* LOCO / leave-center-out；
* strong baseline 对比；
* batch / center / purity 敏感性分析；
* label permutation / negative control。

产物：

* shared module；
* HCC-specific module；
* module-response association；
* strong baseline report。

## Level 2：Mechanistic support

回答：

> 这个模块是否真的像一个可解释的 TME 机制，而不是统计噪声？

证据来源：

* TF-target / pathway 支持；
* perturb prior 支持；
* ligand-receptor / cell-cell interaction 支持；
* spatial / tissue co-localization；
* 与已知 PD-1 或 PD-1+X 生物学一致。

产物：

* mechanism-supported module；
* candidate bottleneck target；
* X mechanism class。

## Level 3：Translational nomination

回答：

> 这个模块或靶点是否值得进入实验或后续转化验证？

证据来源：

* 外部 bulk / clinical anchor；
* spatial / tissue adjudication；
* perturbation direction consistency；
* druggability；
* 最小体外实验可行性。

产物：

* Top module；
* Top-1 / Top-3 target；
* Top PD-1+X mechanism class；
* minimal wet-lab validation plan。

这三层足够了。再多就不是严谨，是人类层级崇拜。

---

# 5. 数据架构

v6.1 使用四类数据，不再把它们全都塞进一个平等的大锅。

## 5.1 主训练数据：ICI scRNA / single-cell-level data

用途：

* 学习 TME 状态；
* 建立 PD-1 sensitivity anchor；
* 发现共享/特异免疫模块；
* 构建 patient-level module scores。

要求：

* response label 尽量明确；
* timepoint 尽量明确；
* metadata 必须记录治疗方案；
* HCC 为第一深描癌种；
* 多癌种用于共享模块学习。

## 5.2 基础外部数据：bulk ICI cohorts

用途：

* 验证 module score；
* 验证 patient-level M-signals；
* 验证 HCC 转化相关性。

v6.1 仍然保留 IMbrave150 这类 bulk anchor，但它不能再只是“趋势展示图”。v6 已经明确要求把 bulk projection 升级为外部验证，而不是继续停留在 demo 层。

## 5.3 关键裁判数据：spatial / tissue-level data

v6.1 上调空间数据地位。

用途不是炫模型，而是裁判核心 TME 机制：

* 髓系抑制是否真的贴近肿瘤巢？
* T cell exclusion 是否空间成立？
* APC niche 是否与 responder-like module 共定位？
* LRG 是否具有空间相邻基础？
* VEGF / stromal remodeling 是否解释 PD-1+X 增益？

空间数据不一定进入主模型，但必须进入主证据链。

## 5.4 干预先验数据：LINCS / perturb-seq / scPerturb

用途：

* 作为方向性软约束；
* 支持 target prioritization；
* 支持 X mechanism class；
* 不直接替代 patient TME 证据。

v5.3 已经强调扰动类先验可作为因果软约束但不替代数据驱动，v6.1 继续保留这个限制。

---

# 6. 模型与分析主线

v6.1 不再强调 retrieval / memory bank。它们只保留为解释层。

主线改成六步。

## Step 1：Immune-state measurement

目标：

把每个患者样本表示为可解释的 TME 状态。

输入：

* scRNA；
* metadata；
* cell-state annotations；
* pseudobulk；
* patient-level aggregation。

输出：

* cell fraction；
* pseudobulk module matrix；
* immune state scores；
* basic TME feature table。

重点：

* 髓系必须单独 QC；
* 不允许直接把 MVP 的手工 TME-MOA 当主证据；
* atlas / scVI 只作为表征底座，不作为主创新。

MVP 的价值主要在数据组织、patient-level 表示和中间表，而不是最终机制层；gap 分析也明确指出，MVP 还没有进入 causal GRN、M-signals、ICSmini/e_X-mini 和投稿级验证闭环。

## Step 2：Strong baseline first

这是 v6.1 新增的硬门槛。

在任何复杂模型前，必须建立强基线。

基线至少包括：

1. **Cell fraction baseline**
   CD8、Treg、myeloid、DC、NK、B/plasma 等比例。

2. **Curated immune signature baseline**
   IFN、cytotoxicity、exhaustion、antigen presentation、Treg、myeloid suppression。

3. **Pseudobulk pathway baseline**
   GSVA / AUCell / PROGENy / TF activity。

4. **Simple ML baseline**
   elastic net、random forest / XGBoost、logistic / Cox。

5. **Static module baseline**
   WGCNA-like module 或 graph-free module。

6. **scVI latent baseline**
   latent + simple classifier，不允许直接拿复杂模型欺负空气。

验收：

* baseline 必须经过同样 split；
* 同样 LOCO / leave-center-out；
* 同样校准和置信区间；
* 同样 negative control。

复杂模型只有在以下任一方面超过强基线，才有资格进入主文：

* 预测性能更好；
* 泛化更稳；
* 机制解释明显更强；
* 外部锚点更一致；
* wet-lab candidate 更清晰。

否则主线要诚实转为“机制解释优先”，别硬说自己是 predictive model。v5.3 原本也设定了“预测不够亮眼时转向机制修复与可达性”的降级路径。

## Step 3：Shared/specific immune module discovery

目标：

发现跨癌种共享和 HCC 特异的免疫模块。

与 v5.3 的不同：

* v5.3 更强调 shared-specific causal GRN；
* v6.1 更强调 **先稳定模块，再谨慎升级 causal claim**。

也就是说，不再一上来宣称 causal GRN，而是分阶段命名：

* inferred immune module；
* perturbation-constrained module；
* causal-candidate module；
* experimentally prioritized module。

输出：

* shared immune modules；
* HCC-specific modules；
* module-response association；
* module-treatment-context association；
* module stability report。

v5.3 的主干仍然保留：共享层在免疫核心基因集，特异层保留器官性残差，扰动方向作为软约束。

## Step 4：PD-1 anchor → PD-1+X repair logic

这是 v6.1 的生物学核心。

先用 PD-1 / ICI monotherapy anchor 学：

* responder module；
* primary resistance module；
* immune-cold module；
* myeloid/Treg suppressed module；
* antigen-presentation deficient module。

然后在 PD-1+X 中问：

> X 修复的是哪个 PD-1 单药不足以解决的障碍？

例子：

* PD-1 anchor 显示 T cell 已存在但功能受抑 → X 可能是 myeloid/Treg/stromal remodeling；
* PD-1 anchor 显示 APC/antigen presentation 不足 → X 可能是 IFN/APC restoration；
* PD-1 anchor 显示 vascular/stromal exclusion → X 可能是 anti-VEGF/TKI-like remodeling；
* PD-1 anchor 显示 T cell priming 不足 → X 可能是 co-stimulation / vaccine-like / antigen release support。

这样 monotherapy 和 combination 是一条逻辑链，不是两个课题。

## Step 5：Spatial/tissue adjudication

v6.1 把空间数据作为关键裁判层。

空间层要回答四个问题：

1. 模块是否存在于正确组织位置？
2. 关键细胞是否空间相邻？
3. ligand-receptor 是否有空间共定位基础？
4. HCC-specific barrier 是否真的具有组织生态意义？

输出：

* spatial module map；
* T cell exclusion / infiltration score；
* myeloid-tumor neighborhood score；
* APC niche score；
* LRG spatial support；
* spatial adjudication report。

空间不一定进主模型，但没有空间/组织证据支撑的 TME claim 不应过度上升为主机制结论。

## Step 6：Mechanism nomination, not clinical recommendation

最终输出不是“给病人推荐药”，而是：

* Top shared module；
* Top HCC-specific module；
* Top repair mechanism class；
* Top-1 / Top-3 candidate target；
* PD-1+X hypothesis list；
* minimal wet-lab validation plan。

组合排序仍可保留，但命名为：

**mechanism-prioritized PD-1+X nomination**

不是 clinical recommendation。人类医学伦理和监管不会因为我们图画得漂亮就自动让路，烦，但合理。

---

# 7. v6.1 的简化 M-signals

v5.3 的 M-signals 三主读数 FG / mICS / MRI 是有价值的，但 v6.1 要压缩复杂度。v5.3 原本还列了九个辅读数，这些东西很容易膨胀成指标动物园。

v6.1 只保留三主两辅。

## 主读数 1：mIMS

**module immune mechanism score**

替代原先分散的 module score / mICS 初级版本。

含义：

* 某患者在某个免疫机制模块上的激活/抑制程度；
* 是所有后续分析的基础。

## 主读数 2：MRI

**mechanism restoration index**

保留。

含义：

* 某个 X 机制类是否把患者状态推向 PD-1 responder-like module；
* 重点服务 PD-1+X nomination。

## 主读数 3：mICS

**module intervention candidate score**

保留，但从复杂反事实降级为第一阶段可执行版。

整合：

* module centrality；
* perturb direction；
* druggability；
* spatial/tissue support；
* external anchor consistency。

## 辅读数 1：ESR

edge / evidence support rate

用于评估模块结构证据。

## 辅读数 2：LRG-spatial

ligand-receptor gain with spatial support

用于 TME 互作和空间证据。

暂时下线或作为补充：

* FG；
* RT；
* VA；
* FΔH；
* BNS；
* SPP；
* TAM-ΔM2→M1。

如果 pre/post 质量足够，再把 FG/RT/FΔH 作为 secondary dynamic analysis。v5.3 已经允许在前后样本不足或跨中心差异大时关闭 SS-CGF，只保留静态因果 GRN + 先验；v6.1 更坚定地采用这个关闭逻辑。

---

# 8. 简化后的证据输出

每个候选机制模块输出一张证据卡。

## Evidence Card 字段

1. **Module name**
2. **Biological interpretation**
3. **PD-1 anchor association**
4. **PD-1+X repair hypothesis**
5. **Cross-cohort robustness**
6. **Strong baseline comparison**
7. **Perturbation support**
8. **Spatial/tissue support**
9. **Bulk/external support**
10. **Candidate target**
11. **Experimental feasibility**
12. **Evidence level：L1 / L2 / L3**

这比复杂证据等级表更容易进入论文和答辩。

---

# 9. 验证体系

v6.1 的验证分四类。

## 9.1 强基线验证

必须回答：

> 新模块是否比 cell fraction / signature / pathway / simple ML 更好？

指标：

* AUC / AUCPR；
* calibration；
* decision curve，可选；
* effect size；
* confidence interval；
* leave-center-out；
* LOCO。

## 9.2 稳健性验证

包括：

* label permutation；
* center prediction audit；
* tumor purity sensitivity；
* cell number sensitivity；
* max-center removal；
* response-label mapping sensitivity；
* timepoint sensitivity。

## 9.3 机制验证

包括：

* pathway / TF support；
* perturb prior support；
* spatial co-localization；
* LRG spatial support；
* module-target consistency。

## 9.4 转化验证

包括：

* bulk external anchor；
* HCC-specific external cohort；
* Top target 最小体外验证；
* 若可行，再进入 PD-1/PD-L1 联用验证。

v5.3 原本要求 LOCO、留中心、消融、结构阈值、ZSL、跨折稳定性和外部锚点；v6.1 保留这些精神，但压缩成更可执行的验证包。

---

# 10. MVP 继承原则

v6.1 继续采用审慎继承。

## 直接继承

* cohort registry 思路；
* sidecar-first 数据组织；
* patient/sample metadata；
* pseudobulk；
* scVI atlas 作为 QC / 表征底座；
* patient-level 表示的工程经验；
* scope lock 与输出契约。

## 降级继承

* archetype：展示层；
* patient trajectory：展示层 / secondary analysis；
* axis：baseline / sanity check；
* TME-MOA：语义标签层；
* virtual drug shift：只作为演示灵感，不作为正式证据。

## 不继承为主结果

* 手工 MOA 直接当机制；
* heuristic counterfactual；
* pseudo-archetype bulk assignment；
* demo-oriented ranking。

当前 gap 已经明确指出，MVP 完成的是故事骨架和 patient-level 原型，但核心算法、M-signals、排序和验证闭环仍然差距很大。

---

# 11. 分阶段施工

## Phase 1：可信数据底座

输出：

* cohort registry；
* metadata master；
* treatment context flag；
* response label dictionary；
* pseudobulk；
* cell fraction；
* state QC，尤其是髓系 QC。

验收：

* metadata 不乱；
* timepoint 不乱；
* treatment flag 可用；
* 髓系状态不靠玄学。

## Phase 2：强基线包

输出：

* fraction baseline；
* signature baseline；
* pathway baseline；
* elastic net / XGBoost baseline；
* static module baseline；
* baseline comparison report。

验收：

* 所有高级模型必须和它比较；
* 如果强基线已经足够好，高级模型必须证明自己不是装饰品。

## Phase 3：shared/HCC-specific module discovery

输出：

* shared modules；
* HCC-specific modules；
* module stability；
* module-response association；
* PD-1 anchor module。

验收：

* 至少一个 shared PD-1 sensitivity module；
* 至少一个 HCC-specific resistance/barrier module；
* 跨队列稳定性可报告。

## Phase 4：空间/组织裁判层

输出：

* spatial support report；
* myeloid-tumor neighborhood；
* T cell exclusion/infiltration；
* APC niche；
* LRG-spatial support。

验收：

* 至少一个核心模块获得空间或组织层支持；
* 如果空间不支持，主机制结论降级。

## Phase 5：PD-1+X mechanism nomination

输出：

* repair mechanism class；
* mICS；
* MRI；
* candidate target；
* evidence card；
* top PD-1+X hypotheses。

验收：

* 每个 nomination 必须连接到 PD-1 anchor；
* 不能凭 LINCS reversal 单独提名；
* 不能凭 marker 均值单独提名。

## Phase 6：外部锚点与实验计划

输出：

* bulk external validation；
* Top module / target；
* minimal wet-lab plan；
* final evidence table。

验收：

* 至少一个外部队列支持；
* 至少一个实验可测 candidate；
* 形成可投稿主图逻辑。

---

# 12. 论文故事线

v6.1 的故事线可以压成 6 张主图。

## Fig. 1：PD-1-centered study design

展示：

* PD-1 monotherapy anchor；
* PD-1+X extension；
* scRNA / bulk / spatial / perturb / experiment；
* 三档证据等级。

## Fig. 2：Strong baseline and immune-state measurement

展示：

* cell states；
* strong baseline comparison；
* PD-1 response association；
* HCC vs pan-cancer。

## Fig. 3：Shared and HCC-specific immune modules

展示：

* shared PD-1 sensitivity module；
* HCC-specific resistance module；
* module heatmap；
* module stability。

## Fig. 4：Spatial adjudication of TME mechanism

展示：

* T cell exclusion；
* myeloid-tumor neighborhood；
* APC niche；
* LRG spatial support。

这是 v6.1 相比 v6 更重要的一张图。

## Fig. 5：PD-1+X repair logic

展示：

* PD-1 anchor；
* residual barrier；
* X mechanism class；
* MRI / mICS；
* candidate target。

## Fig. 6：External validation and experimental nomination

展示：

* bulk external anchor；
* Top target；
* minimal validation design；
* evidence card summary。

---

# 13. 关闭条件

v6.1 的关闭条件更简单。

## 关闭条件 1：复杂模型不超过强基线

结果处理：

* 取消复杂模型主叙事；
* 论文转向 robust immune module + spatial TME mechanism。

## 关闭条件 2：空间不支持核心 TME 机制

结果处理：

* 机制结论降级；
* 不做强 PD-1+X nomination；
* 保留为 association module。

## 关闭条件 3：PD-1 anchor 不稳定

结果处理：

* 不进入 PD-1+X 逻辑；
* 先重做 cohort / label / treatment flag；
* 必要时改成 HCC-only mechanism paper。

## 关闭条件 4：perturb prior 与 patient TME 方向冲突

结果处理：

* perturb 只作为反证；
* 不输出该 X nomination；
* 优先相信 patient-level 与 spatial/tissue 证据。

这比 v6 的多关闭条件更短，也更能执行。别把失败路线写成百科全书，失败本身已经很烦了。

---

# 14. 最终产出

## 表格

* `cohort_registry_v6_1.csv`
* `treatment_context_flags.csv`
* `strong_baseline_results.csv`
* `immune_module_scores.csv`
* `shared_specific_modules.csv`
* `spatial_adjudication_results.csv`
* `M_signals_v6_1.csv`
* `PD1X_nomination_cards.csv`
* `external_validation_results.csv`
* `wetlab_priority_list.csv`

## 文档

* `scope_lock_v6_1.md`
* `baseline_report.md`
* `module_discovery_report.md`
* `spatial_adjudication_report.md`
* `PD1_anchor_to_X_repair_logic.md`
* `minimal_validation_plan.md`

## 主结果

* PD-1 sensitivity anchor module；
* HCC-specific resistance/barrier module；
* spatially supported TME mechanism；
* Top PD-1+X repair mechanism class；
* Top-1 / Top-3 candidate target；
* 最小实验验证计划。

---

# 15. v6.1 的一句话版本

**v6.1 以 PD-1 单药/ICI 单药反应为基础锚点，先用强基线和共享/特异免疫模块定义 ICI 敏感性与残余耐药障碍，再用空间/组织证据和扰动先验判断哪些 X 机制最可能修复这些障碍，最终输出 PD-1+X 机制假设与可实验验证靶点。**

这比 v6 更简单，更像能做完，也更像能被评审读懂。科研方案不是越多层越高级，有时候只是越多层越像人类组织结构图，一眼望去全是部门，没人负责结果。
