你现在执行 v6.2.1 项目的 **Scientific & Engineering Realignment Audit**。

这不是一个新的分析 phase，也不是修复算法、补数据或继续 Phase8 的执行任务。

你的唯一目标是：

**重新阅读当前 scientific plan、roadmap 和最新项目状态报告，从第一性原理恢复项目原本要回答的科学问题；判断实际执行过程在哪里发生了科学目标、统计对象、gate 含义或工程阶段的漂移；然后给出一个明确、可执行但尚不实施的目标校准方案。**

本轮禁止任何实质性计算、模型训练、response analysis、重新筛选 FM、阈值优化或数据修复。

---

# 0. 必须优先阅读的 source of truth

至少完整阅读：

1. `plan.md`
2. `roadmap.md`
3. `V6_2_PROJECT_SCIENTIFIC_LOGIC_AND_BLOCKER_REPORT_20260731.md`

如果项目目录中存在这些文档所指向的最新 Phase6/Phase7/Phase8 decision manifest、final report 或 handoff，可以读取它们以确认“实际执行了什么”，但：

* scientific plan 定义科学目标；
* roadmap 定义预期工程路径；
* 20260731 blocker report 定义当前实际状态；
* 更早版本只能用于追溯漂移来源，不能覆盖上述 source of truth。

不要引用外部文献。
不要重新设计一个新项目。
不要因为当前结果困难而主动降低原始科学问题。

---

# 1. 首先恢复项目真正的科学问题

先不要看当前哪个 FM 通过了 gate。

从 plan 本身回答：

## 1.1 项目的研究对象是什么？

明确区分：

* patient immune-tumor state；
* response-blind immune programs/modules；
* shared vs context-modulated failure structure；
* residual barriers；
* responder-compatible direction；
* X-class repair hypothesis。

说明 FM/module 在整个体系中的正确角色究竟是什么。

重点检查：

**FM 是否本来应该作为 patient state 的多维表示坐标，而不是预先筛选后的 response predictors。**

## 1.2 项目的核心统计问题是什么？

明确项目是否原本试图研究：

* 多模块联合状态；
* 模块间 interaction；
* shared/context modulation；
* 非线性或 module-factorized representation；
* multi-barrier profile；
* repair trajectory；
* X-class alignment；

还是本来就只是：

`筛出若干 FM → 用线性模型预测 response → 找最显著 FM`

如果后者不是 plan 的本意，必须明确指出。

## 1.3 dominant barrier 的原始含义是什么？

判断它应该是：

* 模型和数据支持后得到的可能结果；

还是：

* 在进入模型之前必须先选出的单一 feature。

允许的最终结果是否本来应该包括：

* single dominant barrier；
* multiple co-dominant barriers；
* coarse barrier group；
* diffuse/mixed barrier state；
* abstention / unresolved。

---

# 2. 重建当前“实际完成了什么”

不要根据 phase 名称判断进度。

依据真实输出，把当前成果拆成以下层级：

### A. Data/interface

哪些数据、metadata、patient-timepoint、expression layer 已经可靠？

### B. Immune-state measurement

哪些 cell-state / patient-timepoint features 已经可靠？

### C. Response-blind module discovery

目前 frozen FM/module 是怎样得到的？

是否读取了 response？

### D. Module measurement reliability

每个 FM 当前真正证明了什么？

区分：

* 可以稳定评分；
* measurement uncertainty 较高；
* route disagreement；
* coverage/completeness 问题。

### E. Module identifiability

目前真正证明的是：

* 模块有没有信息；

还是：

* 模块的独立效应是否可以从相关模块中被单独归因。

必须严格区分这两个问题。

### F. Response / PD-1 anchor

到底有没有合法完成正式 response modeling？

不要把旧 Phase8 或被废弃结果当作当前结果。

### G. SRB / context / transport / repair / ZSL / spatial evidence

分别判断：

* 已完成；
* 部分完成；
* 仅有前置资产；
* 尚未开始；
* 因上游阻隔而 blocked。

最终给出一个**基于科学证据链，而不是 phase 编号的真实项目进度**。

---

# 3. 专门审计 FM gate 是否发生概念漂移

这是本次审计的核心。

检查当前 pipeline 是否发生了下面这种变化：

原始目标：

`多个 response-blind modules`
→ `共同描述 patient immune state`
→ `联合 representation / context model`
→ `response association`
→ `再判断哪些 module 或 module group 可以获得独立机制归因`
→ `repair`

实际执行是否逐渐变成：

`多个 FM`
→ `逐个 measurement gate`
→ `逐个 identifiability gate`
→ `只有单独可辨识 FM 才允许进入后续 response model`
→ `最后只剩 FM07`
→ `FM07 与 response`

如果是，必须判断这是不是 scientific/engineering drift。

特别回答：

1. measurement reliability gate 应决定什么？
2. identifiability gate 应决定什么？
3. 哪些 gate 应该影响 **representation eligibility**？
4. 哪些 gate 只应该影响 **claim eligibility**？
5. 共线是否意味着 module 应从模型中删除？
6. top-decile Jaccard 略低于阈值是否意味着 module 不再具有联合建模价值？
7. FM07 是“唯一有生物信息的 module”，还是“目前唯一允许独立机制归因的 module”？

不得通过修改阈值来解决这个问题。

---

# 4. 检查 plan 与 roadmap 本身是否存在内部矛盾

重点检查：

* Barrier Identifiability Gate 在 plan 中需要哪些输入；
* roadmap 中它实际被放在什么位置；
* SRB Stage-A 与 identifiability 的先后顺序是否一致；
* gate 原本是为了保护 dominant-barrier claim，还是为了筛选 SRB 输入；
* Phase8 contract 是否把一个 claim-level gate 错误转化成了 feature-level hard gate。

如果 plan 和 roadmap 本身存在矛盾：

不要偷偷选择一个版本。

明确写出：

1. 冲突是什么；
2. 哪个解释更符合项目核心科学问题；
3. 为什么；
4. 哪些合同需要版本化 patch。

---

# 5. 判断当前到底存在几类 blocker

不要把所有问题混成一个“Phase8 blocked”。

至少分别判断：

## Blocker A — Governance / contract blocker

例如：

Phase8 frozen expected barrier set 与 Phase7 eligible set 不一致。

## Blocker B — Architecture / scientific drift

例如：

identifiability 被错误地用成 representation feature selection，导致多模块问题退化为单 FM 问题。

## Blocker C — Data / identifiability blocker

例如：

现有 response anchors 在 cancer、drug、endpoint、cohort 上同时变化，shared anti-PD1 effect 难以识别。

## Blocker D — Evidence / replication blocker

例如：

缺乏独立、可交换的 anchor replication 或 HCC baseline replication。

对每个 blocker 标记：

* 是否立即阻止代码运行；
* 是否立即阻止科学 claim；
* 是否可以通过合同 patch 解决；
* 是否需要新增数据；
* 是否必须重构 Phase7–9；
* 是否只是 uncertainty，而不是 failure。

---

# 6. 校准科学目标

在不降低原项目科学野心的前提下，重新定义当前阶段应该追求什么。

应特别考虑以下原则：

### 6.1 Representation 与 Claim 分离

考虑是否应建立两个不同集合：

**Representation set**

用于描述完整 patient immune state。

原则上允许保留：

* 高可靠 module；
* 存在 measurement uncertainty 但仍有信息的 module；
* 共线但可作为 module block/coarse barrier 使用的 module。

它们可以进入：

* SRB Stage-A；
* joint patient-state representation；
* reconstruction；
* context modeling；
* spatial/perturbation mapping。

**Claim-eligible set**

只有达到严格 measurement + identifiability 要求的 module/module group 才能：

* 单独命名为 independent barrier；
* 单独估计独立 effect；
* 做独立 `do(repair)`；
* 单独连接 X-class；
* 支撑强机制 claim。

判断这种分层是否更符合 plan。

### 6.2 不强制 single dominant barrier

判断后续目标是否应从：

`找到一个 dominant FM`

恢复为：

`估计 multi-barrier patient state，并在证据允许时识别 dominant/co-dominant/coarse/unresolved barrier`

### 6.3 Response head 不承担 representation discovery

判断是否应：

* 在 response-blind 数据上先学习稳定联合 representation；
* 再让稀缺 response 数据只用于 anchor/context head；
* 而不是先用 response 可归因性筛 representation。

---

# 7. 校准工程目标

本轮不要实施，只给出应修改的工程边界。

明确判断下面每部分：

### Keep

哪些当前产物仍然完全有效，无需返工？

例如：

* data interface；
* Phase4A/4B；
* frozen response-blind module membership；
* measurement audits；
* identifiability statistics。

### Reinterpret

哪些结果数值本身有效，但含义需要修改？

例如：

* FM01/FM04 sensitivity-only；
* FM03 non-identifiable；
* FM07 eligible-independent。

判断它们是否应该从：

“能否进入模型”

重新解释为：

“能否获得独立 mechanism claim”。

### Patch

哪些 governance/contract 文件需要改，但不需要重跑计算？

### Reorder

哪些 phase 的逻辑顺序需要调整？

尤其检查：

`module discovery`
→ `measurement audit`
→ `SRB Stage-A representation`
→ `joint response/context modeling`
→ `identifiability / attribution`
→ `repair`

是否比当前顺序更合理。

### Rerun later

哪些分析只有在合同校准完成后才值得重新执行？

本轮禁止实际 rerun。

---

# 8. 当前 phase 应重新如何命名

不要机械说“Phase8”。

根据真实证据链，给当前阶段一个更准确的定位，例如：

* module measurement complete;
* representation/claim gate conflated;
* pre-anchor architecture realignment required;

或你认为更准确的描述。

同时回答：

1. 当前已经完成科学链条的哪一步；
2. 哪一步实际上没有完成；
3. 哪一步因为设计漂移而需要重新定义；
4. 校准后下一真正 phase 应该是什么。

---

# 9. 本轮禁止事项

本轮绝对禁止：

* 不运行 response model；
* 不训练 SRB；
* 不重做 module discovery；
* 不改变 frozen FM membership；
* 不降低或提高任何已有阈值；
* 不重新选择 FM；
* 不看 response 后调整 module；
* 不下载新数据；
* 不补空间数据；
* 不做 perturbation；
* 不做 X-class ranking；
* 不生成新的 biological claim；
* 不以“让 Phase8 能跑起来”为优化目标。

**本轮优化目标只有一个：恢复科学问题、统计对象、gate 含义和工程阶段之间的一致性。**

---

# 10. 最终只输出一个报告

输出：

`V6_2_1_SCIENTIFIC_ENGINEERING_REALIGNMENT_REPORT.md`

报告必须使用以下结构：

## 1. Executive verdict

用不超过 10 句话回答：

* 当前真正进度；
* 是否发生 drift；
* drift 的严重程度；
* 当前 immediate blocker；
* 更深层 blocker；
* 是否应该继续现有 FM07-only Phase8。

## 2. Intended scientific architecture

用简洁语言恢复 plan 原始科学逻辑。

## 3. Actual implemented architecture

描述当前实际 pipeline。

## 4. Drift matrix

表格字段：

* component
* intended role
* actual role
* drift
* severity
* consequence

至少覆盖：

* FM/module
* measurement gate
* identifiability gate
* dominant barrier
* Phase8
* SRB Stage-A
* response anchor
* shared/context modeling
* repair
* X-class

## 5. FM-by-FM reinterpretation

对每个 frozen FM 输出：

* measurement status
* identifiability status
* representation eligibility
* independent-claim eligibility
* recommended current role

注意：

这里不是重新筛选 FM，只是重新解释已有结果。

## 6. Blocker decomposition

分别列：

* governance blocker
* architectural drift
* data identifiability blocker
* replication/evidence blocker

## 7. Corrected scientific target

明确下一阶段真正应该回答的问题。

## 8. Corrected engineering target

只写需要修改什么合同、phase 顺序和输入/输出关系。

禁止执行。

## 9. Asset disposition

分为：

* KEEP
* KEEP_BUT_REINTERPRET
* CONTRACT_PATCH_REQUIRED
* LOGIC_REORDER_REQUIRED
* FUTURE_RERUN_REQUIRED
* INVALID/OBSOLETE

## 10. Corrected phase map

给出校准后的：

`current state → next phase → later phases`

不要展开成长 roadmap，只需要把关键依赖关系纠正。

## 11. Decisions requiring explicit approval

只列真正需要人为冻结的科学治理决定。

例如：

* representation set 与 claim-eligible set 是否正式分离；
* identifiability gate 是否从 input gate 改为 attribution/claim gate；
* Phase8 原三-FM contract 是否废止并版本化；
* 是否允许 multi/coarse/unresolved barrier；
* 是否在架构校准前禁止读取 response。

## 12. Final one-paragraph project status

用普通研究者能直接理解的语言重新描述：

“项目现在到底做到哪里，为什么停，校准以后下一步究竟是什么。”

---

# 最重要的判断原则

你不是来“帮现有 pipeline 找理由继续跑”的。

你需要独立判断：

**当前实现是否仍然忠实回答 scientific plan 定义的问题。**

如果发现当前工程设计已经把：

**多模块、context-dependent、potentially nonlinear immune failure state**

退化为：

**筛一个或几个独立 FM 再做 response association**

必须明确指出，并追踪漂移发生在哪一个 gate、contract 或 phase transition。

同时也不要反向过度修正：

* identifiability audit 本身仍然有价值；
* measurement reliability audit 本身仍然有价值；
* 当前 Phase7 数值结果不应因为架构漂移而被废弃；
* 需要纠正的是这些结果在 pipeline 中承担的角色，而不是为了恢复多维模型而取消质量控制。

最终目标是：

**保留已经获得的可靠测量和可辨识性证据，同时恢复 v6.2 原本的多模块 patient-state、shared/context modulation、repair 和 X-class 科学结构。**
