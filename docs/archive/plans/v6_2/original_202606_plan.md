# 科研方案 v6.2

**PD-1-anchored causal immune-mechanism framework with counterfactual repair scoring and spatial-niche adjudication**

中文名：

**v6.2：以 PD-1 为锚的因果免疫机制框架——模块级反事实修复评分与空间 niche 裁判**

---

## 0. v6.2 相比 v6.1 的核心变化

v6.2 不推翻 v6.1，而是把 v6.1 在执行中暴露的"主链断裂 + 方法学只到关联层"两大问题，结合学界近半年方向，做四处升级、一处收紧、一处明确退路。

### 0.1 四处升级（提升发表竞争力）

1. **从"关联模块"升级为"因果模块图 + 反事实修复"**。v6.1 止步于 inferred module 的关联证据；v6.2 把核心方法学定义为 **PASCAR 三层栈**：因果结构层（共享/特异因果模块图）→ 反事实干预层（模块级个体处理效应/ITE 的修复评分）→ 多模态证据裁判层。这正对应近半年成熟的 CINEMA-OT / scDRP（反事实 ITE，2026）/ GRIT（最优传输因果 GRN，2025）路线，是可直接 adapt 而非自造的方法学硬通货。

2. **空间从"裁判层"升级为"共主轴"**。近半年 ICI 抗性研究的前沿已明确转向**空间 niche 与细胞互作**（ISCHIA niche 构建、GraphTME 图预测、responder 中 C1QC⁺巨噬-CD4 T niche、CAF 介导的互作破坏）。v6.2 把空间 niche 分析提升为与模块发现并列的一级发现+证据轴，做 niche 级（而非仅 spot 级）的定量裁判。

3. **引入"虚拟细胞/in-silico 扰动"作为可选交叉验证引擎**。近半年单细胞基础模型（STATE@Arc、Tahoe-x1、X-Cell、scOTM；Virtual Cell Challenge@Cell 2025）使 in-silico perturbation（ISP）成为热点。但同时有强批判证据（大规模 CRISPRi 中仅约 41% 扰动有可测全转录组效应；scFM 未必稳超简单基线）。因此 v6.2 把 scFM **明确定位为可选的表征底座 + 第二套 ISP 引擎**，与最优传输反事实、扰动先验做**三角验证（triangulation）**，而**不押注为主创新**。

4. **新增"证据三角化"为命名方法学贡献**。一条机制提名只有在 {因果方向、反事实 ITE、扰动先验、空间 niche、外部 bulk} 中 ≥2 条正交证据一致时才能升级——把 robustness 做成算法的一等输出，直接迎合顶刊对可重复性的审查口味。

### 0.2 一处收紧（提升工程可行性）

5. **一切在模块级（module-level）而非基因级展开，且优先 adapt 已发表方法**。因果图、反事实、niche、ISP 全部在数十–数百个免疫机制模块的维度上运行；因果用不变性筛选（ICP/IRM 风格）+ 扰动软约束，反事实直接 adapt CINEMA-OT/scDRP 的条件最优传输，niche 直接用 ISCHIA 风格构建。这把 v5.3/v6.1 一直没落地的"causal GRN 野心"压缩为"可在时间窗内做完"的工程量，同时保留方法学新颖性（新颖性在于**整合架构 + PD-1 锚生物学 + HCC 转化**，而非单算法）。

### 0.3 一处明确退路（双轨 PD-1）

6. **PD-1 anchor 双轨制度化**。v6.1 因 anchor 混杂被封锁导致主链断裂。v6.2 把 anchor 救援与 HCC-first 退路写成制度：用三个 GEO 开放的干净 anti-PD1 单药队列（BCC/NSCLC/HNSCC，均带 R/NR + scTCR）救援 anchor；在 Milestone 处依救援门一次性裁决主叙事（R1 全栈 PD-1 中心 / R2 HCC-first 因果转化 / R3 证据图方法学头条），三轨共用 PASCAR 栈。

v6.2 仍继承 v6.1/v5.3 主轴：多癌种免疫治疗、共享+特异免疫机制、强基线前置、阶段门控、机制提名而非临床推荐。

---

## 1. 总目标

v6.2 的总目标是：

**在多癌种 ICI 数据中学习可复现的共享/特异 TME 免疫机制模块，将其升级为以 PD-1 单药反应为锚的因果模块图；在图上用模块级反事实（个体处理效应）量化"修复某机制把患者推向 PD-1 responder-like 多远"，并由空间 niche、扰动/虚拟细胞 in-silico 扰动与外部 bulk 的三角证据裁定，最终输出可被空间组织证据与最小湿实验支持的 PD-1+X 机制假设与候选靶点。**

输出名称仍为 **PD-1+X mechanism nomination framework**（机制假设提名系统），明确排除"临床推荐系统"。

---

## 2. 核心科学问题

v6.2 聚焦四个问题（前三沿用 v6.1，第四为升级新增）。

### Q1：什么是跨癌种可复现的 PD-1 敏感性**因果**基础模块？

不再只问"哪些模块关联敏感性"，而是问"哪些模块构成跨癌种**环境不变（invariant）的因果骨架**"——把"共享"定义为因果不变性，而非聚类约定。

### Q2：HCC 中主要的**特异因果屏障**是什么？

承接近半年 HCC 抗性生物学，重点检验髓系（S100A9⁺CD14⁺单核、TREM2⁺巨噬、巨噬源血管生成）、基质（CAF）、血管（VEGF）是否构成 HCC 特异的因果残差边。

### Q3：哪些 X 机制类的**反事实修复效应**最大？

把"X 修复了什么"量化为模块级个体处理效应：对模块施 `do(repair)`，患者沿 responder-like 流形的可达位移（MRI）。X 机制类按反事实修复效应 + 可干预性排序。

### Q4（新增）：这些因果屏障是否组织为**免疫抑制 niche**，且 X 修复是否预测 niche 解构？

把空间 niche 升为科学问题：核心因果屏障是否在组织中以 niche 形式存在（如髓系-肿瘤巢、T 细胞排斥带、CAF 屏障）；反事实修复方向是否对应 niche 的空间解构。

---

## 3. 治疗结构：PD-1-centered continuum（双轨制度化）

沿用 v6.1 三层连续体，新增双轨判决协议。

* **主分析层 PD1_ICI_anchor**：学习基础 ICI 敏感性的因果锚。救援队列见 §5。
* **扩展解释层 PD1X_extension**：解释 X 修复的残余障碍（含 anti-VEGF/TKI/lenvatinib）。
* **敏感性分析层 high_confounding_support**：TACE/HAIC/radiotherapy 等强混杂，仅支持/反证。
* **外部锚定层 external_anchor_only**：bulk/spatial/clinical 外部验证。

**双轨判决**：anchor 救援门（cohort-response 混杂降至 MEDIUM 以下且 within-cohort 方向一致）通过 → 锁 R1 PD-1 中心主叙事；不过 → 转 R2 HCC-first，PASCAR 栈不变，仅把反事实的 responder 流形从 anchor 切换为 HCC responder-like / 稳健模块流形。

---

## 4. 证据等级：三档 + 三角化

沿用 L1/L2/L3，新增三角化为升级前置。

* **L1 Robust association**：多队列复现、LOCO/leave-center-out、强基线对比、置换/负对照。
* **L2 Mechanistic support**：**因果方向**（Layer A）、**反事实 ITE 方向**（Layer B）、TF/pathway、扰动/ISP 先验、**空间 niche 共定位**、ligand-receptor。
* **L3 Translational nomination**：外部 bulk/clinical 锚、空间裁判、druggability、最小湿实验可行性。

**三角化规则（新增、命名贡献）**：任一 candidate 升 L2/L3 前，必须在 {因果方向、反事实 ITE、扰动/ISP 先验、空间 niche、外部 bulk} 中至少 2 条正交证据方向一致；任一关键冲突（如反事实方向与空间 niche 相反）自动降级并记录。

---

## 5. 数据架构

四类数据沿用 v6.1，新增精确外部队列清单与可选 scFM 表征。

### 5.1 主训练数据：ICI scRNA（含 anchor 救援）

* HCC 深描（第一癌种）+ 多癌种共享。
* **anchor 救援（GEO 开放，决定 R1）**：`GSE123813`（Yost BCC/SCC，anti-PD1 单药，pre/post，R/NR，+scTCR）、`GSE176021`（Caushi NSCLC，新辅助 nivolumab 单药，+scTCR，560K T）、`GSE200996`（Luoma HNSCC，新辅助 nivolumab，含联合臂须剥离，+scTCR）。
* HCC 深描补强：`GSE149614`（10 例 HCC scRNA）。
* 可选受控上行：dbGaP `phs002065`（Bi RCC）、EGA `EGAS00001004290-92`（Braun RCC）。

### 5.2 基础外部数据：bulk ICI cohorts（升级为外部锚，非趋势图）

* HCC：`GSE215011`（nivolumab 单药）、`GSE235863`（anti-PD1+lenvatinib，配对）、`GSE140901`（nivo/pembro±ipi NanoString）、`GSE279750`（anti-PD-L1 联合）、IMbrave150（已有，升级为模块级投影+CI）。
* 泛癌：`GSE78220`（Hugo）、`GSE91061`（Riaz）、`PRJEB23709`（Gide）。
* 投影方式：模块级 score 投影 + 置信区间 + 失败案例报告。

### 5.3 关键裁判+发现数据：spatial / niche（升级为共主轴）

* 已暂存 `GSE238264`（HCC 空间，support_only）；需定向补 ≥1 套 ICI 相关 HCC 10x Visium（含核/侵袭前沿/基质分区）。
* 可选受控对照：Pozniak 黑色素瘤 `EGAS00001006488`（scRNA+Visium+GRN，作 A 层方法对照）。
* 分析：niche 构建（ISCHIA 风格）+ 细胞互作图 + T 细胞排斥/浸润 + 髓系-肿瘤邻近 + APC niche + LRG 共定位 + 置换零控。

### 5.4 干预先验数据：扰动 + 虚拟细胞 ISP（升级为可三角验证）

* 本地扰动：`GSE133344`、`GSE90063`、`GSE193736`、`GSE306429`、LINCS-L1000（构建 perturb-prior 矩阵）。
* 可选 scFM ISP 引擎：在算力允许下用公开单细胞基础模型对核心模块做 in-silico 扰动，与最优传输反事实、扰动先验做三角验证；**须配 leave-cancer 混杂审计**，不作主创新。

### 5.5 可选表征底座：单细胞基础模型 embedding

仅用于跨癌种/跨队列 batch 协调（提升 shared 因果层对齐），不作主表征；须审计"学到的是癌种还是免疫态"。

---

## 6. 模型与分析主线：PASCAR 三层栈

v6.2 的方法学核心。每层设 v0（快速可防守）/v1（顶刊创新）/v2（野心），层间固定数据契约。

### Step A — 结构层：Shared-Specific Causal module Graph（SS-CGF）

* 输入：module_score_matrix（Phase4 已产）+ 治疗上下文 + 扰动先验。
* v0：模块级 graphical lasso（GGM），shared/HCC 池差集 = specific 边，within-cohort meta + 置换零控（关联图，可发）。
* **v1（主推）**：模块级**不变因果发现**（跨癌种作环境，ICP/IRM 风格筛 shared invariant 父集；HCC 特异边用环境交互识别；扰动方向作软约束，方向冲突降权）。
* v2：模块级可微因果发现（NOTEARS/DCDI 变体）+ 扰动作软干预数据。
* 产物：shared causal backbone + HCC-specific residual DAG。

### Step B — 干预层：Counterfactual Repair / module-ITE

* 输入：因果图 + PD-1 anchor（或 HCC responder-like）标签。
* v0：mIMS = 模块激活/抑制标准化；MRI = responder-like 判别方向上的模块扰动梯度一致性；mICS = 规则整合（centrality+perturb/ISP+druggability+niche+external）。产出正式 `mIMS/MRI/mICS` 文件。
* **v1（主推）**：**条件最优传输反事实**（adapt CINEMA-OT/scDRP）——在 anchor responder/non-responder 流形间求位移场；对模块 M 施 `do(repair)` 后计算个体处理效应 ITE = 沿位移到 responder-like 的可达增益 = MRI；mICS = ITE × 可干预性 × 三角证据。
* v2：基于 Layer A 因果图的 structural counterfactual（do-calculus 沿因果路径传播）。
* 产物：module-ITE/MRI、X mechanism class 排序、PD-1+X repair logic。

### Step C — 裁判层：Multimodal Evidence-Graph Triangulation

* 输入：candidate（module, X-class, target）+ 各模态证据。
* v0：异质证据加权聚合 + 置换标定 → L1/L2/L3 + 显式降级规则 → `candidate_module_evidence_matrix.csv`。
* **v1（主推）**：正交证据贝叶斯融合（各模态条件独立证据更新后验，置换/null 校准似然，输出后验+可信区间）+ 三角化规则。
* v2：异质图神经网络排序（仅作敏感性对照，不作主结论）。
* 产物：概率化证据卡、L1/L2/L3 矩阵、冲突清单。

> 重点：髓系必须单独 QC；atlas/scVI/scFM 只作表征底座；不允许把 MVP 手工 TME-MOA 当主证据。

---

## 7. v6.2 的 M-signals（重定义为模块级 ITE）

保留三主两辅，但 MRI/mICS 重定义为反事实量。

* **主读数 1：mIMS**（module immune mechanism score）——模块激活/抑制程度，基础量。
* **主读数 2：MRI**（mechanism restoration index）——模块级**个体处理效应**：反事实修复把患者推向 responder-like 的可达位移。
* **主读数 3：mICS**（module intervention candidate score）——整合 ITE、module centrality、扰动/ISP 方向、druggability、空间 niche 支持、外部一致性。
* **辅读数 1：ESR**（evidence support rate，含三角化命中数）。
* **辅读数 2：LRG-niche**（ligand-receptor gain with spatial niche support）。

暂时下线：FG/RT/VA/FΔH/BNS/SPP/TAM-ΔM2→M1（pre/post 充足时再作 secondary dynamic）。

---

## 8. 证据卡字段（升级）

每个候选机制输出一张证据卡：

1. Module name
2. Biological interpretation
3. **Causal role（shared invariant / HCC-specific residual）**（新增）
4. PD-1 anchor association
5. **Counterfactual repair effect（module-ITE / MRI）**（新增）
6. PD-1+X repair hypothesis（X mechanism class）
7. Cross-cohort robustness（LOCO/leave-center-out）
8. Strong baseline comparison
9. Perturbation / in-silico perturbation support（含三角化）（升级）
10. **Spatial niche support（niche co-localization）**（升级）
11. Bulk/external support
12. Candidate target & druggability
13. Experimental feasibility
14. **Evidence level L1/L2/L3 + triangulation count**（升级）

---

## 9. 验证体系

四类验证沿用 v6.1，新增方法学基准与三角化验证。

* **9.1 强基线验证**：fraction/signature/pathway/elastic net·XGBoost/static module/scVI latent；同 split、同 LOCO、同校准。复杂方法须在预测/泛化/机制解释/外部一致/湿验候选清晰度上至少一项超基线。
* **9.2 稳健性验证**：label permutation、center prediction audit、purity/细胞数敏感性、max-center removal、response-mapping/timepoint 敏感性。
* **9.3 机制验证**：**因果方向稳定性、反事实 ITE 一致性、ISP 引擎间收敛（OT vs scFM vs prior 三角）**、TF/pathway、空间 niche 共定位、module-target 一致。
* **9.4 转化验证**：bulk 外部锚、HCC 外部队列、Top target 最小湿验（IHC/mIF → 共培养 → CRISPR/小分子）。

---

## 10. MVP 继承原则

沿用 v6.1：直接继承数据底座（registry、sidecar、metadata、pseudobulk、scVI/可选 scFM 作 QC/表征）；降级继承 archetype/trajectory/axis/TME-MOA/virtual drug shift 为展示/baseline；不继承手工 MOA、heuristic counterfactual、demo ranking 为主结果。

---

## 11. 分阶段施工

（详见 `roadmap.md`；此处为概览）

* Phase 1 可信数据底座 + anchor 救援队列激活
* Phase 2 强基线包（streamlined）+ 可选 scFM 表征
* Phase 3 shared/HCC-specific 模块 → **Layer A 因果模块图**
* Phase 4 PD-1 anchor（双轨门）+ **Layer B 反事实修复/ITE**
* Phase 5 **空间 niche 定量裁判**（共主轴）
* Phase 6 扰动 + ISP 三角验证 → **Layer C 证据图三角化 + L1/L2/L3 矩阵**
* Phase 7 外部锚激活 + evidence cards + nomination + 最小湿验
* Phase 8 论文主图 + 复现包

---

## 12. 论文故事线（6 主图，升级）

* **Fig. 1**：PD-1-anchored causal framework 设计（continuum + PASCAR 三层 + 三档证据 + 三角化）。
* **Fig. 2**：强基线 + immune-state measurement（cell states、baseline 对比、HCC vs pan-cancer）。
* **Fig. 3**：**Shared invariant 因果骨架 + HCC-specific 因果屏障**（因果图、稳定性、髓系/基质/血管屏障）。
* **Fig. 4**：**空间 niche 裁判**（免疫抑制 niche、髓系-肿瘤邻近、T 细胞排斥、APC niche、CAF 屏障、LRG-niche）——升级为重头图。
* **Fig. 5**：**PD-1+X 反事实修复逻辑**（anchor → residual barrier → module-ITE/MRI → X mechanism class → 三角证据）。
* **Fig. 6**：外部验证 + 实验提名（bulk 外锚、Top target、最小湿验设计、evidence card 汇总）。

---

## 13. 关闭条件

* **关闭 1：因果/反事实不超过关联强基线** → 取消因果主叙事，转 robust module + spatial niche mechanism。
* **关闭 2：空间 niche 不支持核心屏障** → 机制降级，不做强 PD-1+X nomination。
* **关闭 3：PD-1 anchor 不稳定** → 转 R2 HCC-first（双轨退路），PASCAR 栈不变。
* **关闭 4：ISP/扰动与患者 TME 反事实方向冲突** → 该 X 不提名，优先信 patient-level + spatial niche。
* **关闭 5：scFM ISP 不收敛或不稳** → scFM 退为可选佐证，三角化以 OT 反事实 + 扰动先验为准（直接采纳 scFM 批判性证据的防御性设计）。

---

## 14. 最终产出

### 表格
* `cohort_registry_v6_2.csv`、`treatment_context_flags.csv`
* `strong_baseline_results.csv`
* `module_score_matrix.csv`、`shared_specific_causal_edges.csv`（新增）
* `mIMS_scores.csv`、`MRI_scores.csv`、`mICS_scores.csv`（正式 M-signals 文件）
* `spatial_niche_scores.csv`（升级）、`LRG_niche_support.csv`
* `perturbation_isp_triangulation.csv`（新增）
* `candidate_module_evidence_matrix.csv`（L1/L2/L3 + 三角化）
* `PD1X_nomination_cards.csv`、`tcr_clone_master_v6_2.csv`、`external_validation_results.csv`、`wetlab_priority_list.csv`

### 文档
* `scope_lock_v6_2.md`、`baseline_report.md`、`causal_module_graph_report.md`（新增）、`counterfactual_repair_report.md`（新增）、`spatial_niche_adjudication_report.md`、`PD1_anchor_to_X_repair_logic.md`、`minimal_validation_plan.md`

### 主结果
* shared invariant PD-1 敏感性因果骨架
* HCC-specific 因果屏障模块（髓系/基质/血管）
* 空间 niche 支持的 TME 机制
* Top PD-1+X 反事实修复机制类 + Top-1/Top-3 候选靶点
* 最小湿实验验证计划

---

## 15. v6.2 的一句话版本

**v6.2 把共享/特异免疫模块升级为以 PD-1 单药反应为锚的因果模块图，用模块级反事实（个体处理效应）量化每个 X 机制把患者推向 responder-like 的修复力，再以空间 niche、扰动/虚拟细胞 in-silico 扰动与外部 bulk 的三角证据裁定，最终输出可被空间组织与最小湿实验支持的 PD-1+X 机制假设与可验证靶点。**

比 v6.1 更有方法学锋芒（因果+反事实+niche+三角化，全部对齐近半年前沿），又因"模块级 + adapt 已发表方法 + scFM 仅可选"而更可做完。

---

## 附录：v6.2 调整的学界依据（近半年，2025末–2026中）

| 趋势 | 代表工作 | v6.2 的采纳方式 |
|---|---|---|
| 虚拟细胞 / scFM in-silico 扰动 | STATE(Arc)、Tahoe-x1(3B)、X-Cell、scOTM、Virtual Cell Challenge(Cell 2025) | 可选 ISP 引擎 + 表征底座，三角验证，不押主线（含 41% 可测效应的批判性防御） |
| 因果 GRN + 最优传输反事实 | CINEMA-OT(Nat Methods 2023)、GRIT(Bioinformatics 2025)、scDRP(2026 反事实 ITE) | Layer A 因果 + Layer B 反事实直接 adapt |
| 空间 niche / 细胞互作 | ISCHIA、GraphTME(2025)、C1QC⁺巨噬-CD4 niche(Cell Discovery 2025) | 空间升为共主轴，niche 级定量裁判 |
| HCC 抗性=髓系+基质+血管 | S100A9⁺CD14⁺单核、TREM2⁺巨噬(post-TACE)、巨噬源血管生成(2026)、SIRPα/PD-L1⁺髓系(2025) | HCC 特异因果屏障的具体 marker 假设与 X-class 映射 |
