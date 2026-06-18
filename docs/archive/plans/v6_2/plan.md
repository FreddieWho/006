# 科研方案 v6.2

> 当前执行模式（2026-08-27）：本项目按探索性生物信息学研究推进。本文及其引用的早期“预先规定/冻结/gate”表述仅作历史方案记录，不是新分析的硬约束。执行入口见 `docs/v6_2/EXPLORATORY_MODE_DIRECTIVE_20260827.md`。

# PD-1-anchored immune failure mechanism and X-class repair framework

中文名：

# 以 PD-1 为锚的免疫治疗失败机制与 X-class 修复框架

---

## 0. 文档定位

> Revision note, 2026-06-24（v6.2 生物学补强；依据 `results/v6_2/PLAN_ROADMAP_BIOLOGICAL_REVIEW_20260624.md`，用户批准）：
> 已并入 **B1** 与 **B2**；**B3–B6** 落档为补强审阅挂钩，到达对应 phase 时唤起（见 roadmap 顶部挂钩表）。
> - **B1（HCC 病因分层）**：etiology 作为**可插拔、低权重、HCC-context 限定**的 sensitivity 分层；**默认不进入 pan-cancer 主线 module discovery 与 responder-direction，不得主导 IO 重点**；覆盖不足则关闭并记录缺口。目的是避免在"有/无明确病因癌种"之间引入偏倚。
> - **B2（模块本体补强）**：在模块清单与 cell-state 本体中纳入 TLS/B 细胞体液模块、neutrophil/TAN/NET、tumor-intrinsic WNT/β-catenin immune exclusion、CD8 progenitor(TCF1⁺)/terminal(TIM3⁺) 拆分（仍 response-blind 发现，仅进入注释与稳定性雷达）。

本方案定义一个面向 anti-PD1 / ICI 失败机制与 PD1+X 机制修复假设的肿瘤免疫计算框架。

本研究不以“预测某个患者是否应该接受某种治疗”为目标，也不建立临床用药推荐系统。本研究的目标是识别 anti-PD1 失败背后的可复用免疫屏障机制，学习这些机制如何在不同癌种和器官上下文中被改写，并判断哪些 X-class 机制最可能修复这些屏障。

本方案的核心模型不是泛用大模型，也不是单纯 response classifier，而是一个以模块为结构单位的机器学习框架：

**Structured Repair Backbone, SRB：结构化修复主干。**

SRB 学习 patient-timepoint 的 module-factorized latent representation，在同一主干上完成 anti-PD1 failure module 表征、context modulation、PD-1 anchor、counterfactual repair transport 和 X-class zero-shot alignment。空间数据由于模态独立，不并入 SRB 主干，而由独立的 Spatial Niche Graph Model, SNGM 进行组织生态验证。

---

## 1. 核心科学问题

### 1.1 anti-PD1 失败是多模块免疫生态失败，而不是单一 checkpoint 失败

PD-1/PD-L1 阻断直接作用于 T cell checkpoint 抑制，但临床和组织层面的 anti-PD1 失败常常来自更上游或更外周的免疫生态障碍。

这些障碍包括：

1. tumor antigen / HLA / antigen presentation 不足；
2. IFN / inflammatory signaling 不支持免疫激活；
3. CD8 / NK cytotoxic state 不充分；
4. DC / APC priming 失效；
5. suppressive myeloid、Treg、CAF、vascular、stromal barrier 形成残余屏障；
6. T cell 在空间上无法进入或接触肿瘤巢；
7. 治疗后状态迁移没有进入 responder-compatible immune state；
8. B 细胞 / TLS / 生发中心 体液免疫支持不足（B2）；
9. 中性粒 / TAN / NET 介导的 CD8 排斥与髓系抑制（B2，与 HCC 病因联动）；
10. 肿瘤内在 WNT/β-catenin 等驱动的 immune exclusion（B2，antigen-independent）；
11. CD8 progenitor(TCF1⁺) 与 terminal(TIM3⁺) exhaustion 的失衡（B2）。

因此，本研究把 anti-PD1 response / failure 定义为 patient immune-tumor state 与 responder-compatible module state 之间的关系，而不是一个孤立的二分类标签。

### 1.2 PD1+X 的 X 是 residual barrier repair mechanism

PD1+X 中的 X 不应首先被定义为具体药名，而应定义为机制修复类。X 的作用是修复 PD-1 checkpoint release 之外仍然存在的 residual immune barrier。

例如：

* myeloid suppression 强，X 可属于 myeloid reprogramming；
* antigen presentation / IFN failure 强，X 可属于 APC/IFN restoration；
* T cell exclusion 或 stromal barrier 强，X 可属于 vascular/stromal remodeling；
* tumor-intrinsic immune silence 强，X 可属于 epigenetic immune priming；
* regulatory suppression 强，X 可属于 Treg suppression relief。

所以本研究的核心转化逻辑是：

**anti-PD1 failure module → dominant residual barrier → X-class repair hypothesis。**

### 1.3 HCC 是预指定 index cancer，而不是泛癌代表

HCC 在本研究中是预指定 index cancer 和 organ-context deep-dive。它不是“典型癌种”，也不代表所有癌种。HCC 的作用是检验共享 anti-PD1 failure modules 如何被肝脏免疫耐受、炎症、髓系、血管、CAF、内皮和局部治疗背景改写。

HCC 的研究价值在于：

1. 肝脏是免疫耐受器官；
2. HCC 常伴随慢性炎症、肝硬化、病毒或代谢背景；
3. HCC 中 myeloid、macrophage、endothelial、CAF、VEGF 轴具有明确免疫治疗相关性；
4. PD1+anti-VEGF / TKI-like 组合在 HCC 中具有现实转化意义；
5. HCC 空间组织生态适合检验 immune exclusion、myeloid-tumor niche、vascular/stromal barrier 等机制。

因此，本研究主问题不是“HCC 特异耐药机制是什么”，而是：

**泛 anti-PD1 failure modules 是否存在；这些 modules 如何被 HCC 这样的 organ context 改写；这种改写是否能产生可验证的 X-class repair hypothesis。**

---

## 2. 可检验假设

### H1：跨癌种存在可复用的 anti-PD1 failure / sensitivity modules

不同癌种中 anti-PD1 response 与 failure 的差异应部分落在共享免疫模块上。这些模块可能涉及 cytotoxic T/NK、IFN signaling、antigen presentation、APC function、T cell dysfunction、myeloid suppression、Treg regulation、vascular/stromal exclusion 等。

可观测证据：

* 模块在多个 clean / expanded ICI anchor 中方向稳定；
* leave-one-cancer 和 leave-one-cohort 后仍可复现；
* 不是由 cohort、batch、cancer type、cell fraction 或 response label definition 驱动；
* 比单基因、单 signature、粗 cell fraction 更稳健。

### H2：shared anti-PD1 failure modules 会被 cancer / organ context 改写

共享模块并不要求在每个癌种中完全相同。器官背景、肿瘤类型、治疗史、空间结构和局部免疫生态会改变这些模块的表达、互作和转化意义。

HCC 是预指定 context modulation 检验场。HCC 分析的科学目标是判断共享 anti-PD1 failure modules 在肝脏肿瘤生态中如何被改写，而不是把 HCC 特异性上升为全项目唯一主轴。

可观测证据：

* shared modules 在 HCC 中出现可解释 residual effect；
* residual effect 与 myeloid、vascular、CAF、APC/IFN、Treg、组织空间结构等因素相关；
* HCC 外部数据或空间数据支持这种 context modulation；
* 结果不能被平台、纯度、组织来源或样本质量解释。

### H3：X-class 的合理性取决于它是否修复 dominant residual barrier

一个 X-class 是否值得被提名，不取决于它是否流行，也不取决于它是否在文献中常与 PD1 联用，而取决于它是否能修复 dominant residual barrier。

可观测证据：

* dominant residual barrier 可辨识；
* X-class descriptor 与该 barrier 的 module attributes 匹配；
* perturbation / target evidence 支持该 X-class 的修复方向；
* SRB transport head 预测 X-class-conditioned repair 能把状态推向 responder-compatible direction；
* spatial niche 和 external evidence 不冲突。

### H4：强 TME 机制 claim 需要 spatial niche 支持

TME 机制不是 dissociated scRNA 的表达均值。强机制 claim 应该能在组织生态中找到 carrier cells、空间邻近、niche structure 或 ligand-receptor proximity。

可观测证据：

* candidate module carrier cells 可定位；
* module activity 与 tumor-stroma boundary、T cell exclusion、myeloid-tumor adjacency、APC-T cell niche、vascular/stromal barrier 等空间结构相关；
* ligand-receptor 互作有空间邻近基础；
* spatial niche 与 SRB 得到的 dominant barrier / repair direction 一致。

### H5：结构化 mechanism space 应支持 zero-shot / compositional extrapolation

如果 anti-PD1 failure modules 和 X-class repair ontology 是有效的，则模型应能对未见 X-class、未见 PD1+X 组合或未见机制描述进行外推排序或机制归类。

可观测证据：

* leave-one-X-class / leave-one-combination 下，SRB alignment head 超过 nearest-prototype、attribute matching、additive compositional baseline；
* 模型能 abstain，而不是强行分类；
* ZSL 输出与 repair、perturbation、spatial、external evidence 一致；
* 无 semantic leakage。

---

## 3. 数据角色

### 3.1 Clean anti-PD1 / ICI anchor

用于学习基础 anti-PD1 sensitivity / failure direction，并定义 responder-compatible state。

要求：

* PD1 / PDL1 单药或相对干净 ICI；
* response label 清楚；
* treatment context 清楚；
* timepoint 可解释；
* patient/sample identity 可追踪。

### 3.2 Expanded / support anchor

用于增强 anti-PD1 方向、扩展环境多样性和做 sensitivity analysis。它不能替代 clean anchor，但能增强或削弱 anchor 可信度。

### 3.3 Multi-cancer ICI layer

用于发现 shared anti-PD1 failure modules，并检验跨癌种环境不变性。

### 3.4 HCC index cancer layer

用于检验 shared modules 在 HCC 中的 context modulation，并连接 HCC 中 PD1+anti-VEGF / TKI-like 机制修复假设。

### 3.5 Spatial / tissue-niche layer

用于判断 dominant barrier 是否真实组织为空间 niche。空间层不是 SRB 的训练标签，而是独立组织生态验证。

### 3.6 Perturbation / target prior layer

用于提供 X-class repair direction 的先验证据。perturbation prior 只作为 soft condition 和 coverage gate，不作为真值。

### 3.7 External bulk / clinical anchor layer

用于检验 module direction 和 candidate mechanism 是否在外部队列中具有一致趋势。bulk 外锚不能证明细胞互作或空间机制。

---

## 4. 模型总览

本方案的模型系统由两个核心模型和一个证据决策层构成。

### 4.1 模型一：Structured Repair Backbone, SRB

SRB 是主模型。它把 patient-timepoint features 编码成 module-factorized latent representation，并在同一主干上完成以下任务：

1. 表征 anti-PD1 immune state；
2. 学习 shared vs context-modulated failure modules；
3. 建立 anti-PD1 responder-compatible direction；
4. 估计 dominant barrier 的 counterfactual repair trajectory；
5. 对 X-class / PD1+X 组合进行 zero-shot / compositional alignment。

SRB 不是黑箱 patient embedding。它的 latent 被强制按 module 分块：

每个 module m 对应一个 latent sub-block `z_m`。
患者状态表示为 `Z = {z_1, z_2, ..., z_M}`。

这样做的生物学意义是：

* 每个 latent block 对应一个可解释 immune module；
* dominant barrier 可以落到具体 module block；
* repair transport 可以说明主要移动了哪些 module；
* X-class alignment 可以直接对齐 barrier module 与 mechanism descriptor；
* graph attention 可以作为 module interaction proxy，但不能直接等同因果边。

### 4.2 模型二：Spatial Niche Graph Model, SNGM

SNGM 是独立空间验证模型。它不并入 SRB，因为空间数据与 scRNA / patient-level module 数据通常不完全匹配，硬合并会牺牲可信度。

SNGM 的任务是判断 SRB 识别的 dominant barrier 是否形成真实组织 niche，例如：

* myeloid-tumor niche；
* T cell exclusion niche；
* APC-T cell niche；
* CAF / ECM barrier；
* vascular / endothelial barrier；
* ligand-receptor supported suppressive neighborhood。

### 4.3 证据决策层：Evidence Card

Evidence Card 不是神经网络主创新，而是将 SRB、SNGM、perturbation、external validation 和 tissue feasibility 合并为可审计 claim 的决策层。

它负责：

* 合并正交证据；
* 记录冲突；
* 赋予 evidence tier；
* 决定 candidate 是否进入主文、补充、边界分析或剔除。

---

## 5. 输入表征：Immune-state construction

### 5.1 目标

将原始 scRNA、metadata、pathway、TF、cell fraction、pseudobulk 和 response context 组织成 patient-timepoint feature matrix，为 SRB 提供可学习输入。

### 5.2 输入

* scRNA expression；
* cell annotation；
* patient/sample metadata；
* cancer/cohort/batch；
* treatment context；
* response label；
* timepoint；
* cell-state-specific pseudobulk；
* pathway / TF / signature activity；
* optional TCR features。

### 5.3 方法

第一步是 cell-state harmonization。将不同队列的 annotation 映射到统一层级，包括 T/NK、CD8 effector、CD8 dysfunctional、Treg、NK、DC/APC、inflammatory myeloid、suppressive myeloid、macrophage、B/plasma、CAF、endothelial、tumor-like 等。

第二步是 index-context refinement。对 HCC 中高度相关的 myeloid、macrophage、CAF、endothelial 状态进行可靠性审计。若这些细分状态不稳定，后续不允许写强 HCC context modulation claim。

第三步是 patient-timepoint aggregation。对每个 patient/timepoint 生成 cell-state fraction、cell-state-specific pseudobulk、pathway/TF/signature activity、candidate module features 和 QC covariates。

第四步是 response label environment registration。RECIST/irRECIST、timepoint、药物、cohort、response definition 作为 environment 或 stratification variable，不作为普通 biological feature 直接供模型利用，防止 label-definition leakage。

（B1）对 HCC，etiology（HBV / HCV / alcohol / NASH-MAFLD / 其他）作为**可插拔、低权重**的 sensitivity 分层变量登记，**默认仅在 HCC index-context（见 H2 与 §8.4.2 context head）启用，不并入 pan-cancer 主线 module discovery，亦不得主导 responder-direction**；启用前须评估病因元数据覆盖率，覆盖不足则保持关闭并记录缺口，避免在有/无明确病因癌种间引入偏倚。

### 5.4 输出

* patient-timepoint feature matrix；
* cell-state fraction matrix；
* cell-state-specific pseudobulk matrix；
* pathway / TF / signature matrix；
* cell-state reliability report；
* response-label environment table；
* etiology-stratified sensitivity layer（B1；仅 HCC，可插拔，低权重，默认关闭，覆盖不足不启用）。

### 5.5 可靠降级

如果细分 annotation 不稳定，则降级为 coarse cell-state ontology。
如果 myeloid / stromal 细分不可靠，则合并成 broader myeloid-stromal barrier，不做细粒度机制 claim。

### 5.6 并行比较

并行比较 marker-based harmonization 与 reference-mapping harmonization。比较指标包括 label consistency、marker enrichment、module stability 和 spatial mappability。

---

## 6. Response-blind module discovery

### 6.1 目标

发现稳定、可解释、可评分的 immune modules。module 是 SRB latent 的生物学分块单位，也是后续 repair、ZSL 和 evidence card 的基本语义单位。

### 6.2 输入

* immune-core genes；
* cell-state-specific pseudobulk；
* pathway / TF / signature features；
* patient-timepoint feature matrix；
* 不直接使用 response label 做模块发现。

### 6.3 方法

模块发现采用 response-blind strategy。推荐使用两条并行路线：

1. co-expression / co-activity module route
   使用 WGCNA-like、NMF、consensus clustering、graph community detection 或类似方法发现稳定模块。

2. sparse dictionary refinement route
   使用 sparse autoencoder 或 non-negative dictionary learning 学习可解释模块字典，但不让 response label 直接决定模块结构。

模块需要通过 bootstrap、leave-cohort、leave-cancer、resolution sensitivity 和 module membership stability 评估。

### 6.4 输出

* module membership；
* module score matrix；
* module attributes；
* carrier cell-state；
* pathway / TF / LR annotation；
* module stability；
* module mappability to spatial and perturbation data。

### 6.5 可靠降级

如果 neural dictionary 不稳定，则使用 consensus WGCNA/NMF/graph community modules。
如果高分辨率 modules 不稳定，则合并为 broader immune programs。

### 6.6 并行比较

比较 WGCNA-like modules、NMF/topic modules、sparse autoencoder modules。主要比较模块稳定性、生物学富集、跨队列复现、空间可映射性和 downstream SRB utility。

---

## 7. Barrier identifiability gate

### 7.1 目标

在进行 dominant barrier 判定和 repair transport 前，确认模块之间是否可区分。myeloid、Treg、stromal、IFN、tumor purity 和 inflammation 等信号容易共线，若不先判断可辨识性，dominant barrier 会变成不稳定归因。

### 7.2 输入

* module score matrix；
* SRB Stage-A latent；
* cell-state fractions；
* deconvolution outputs；
* bootstrap resamples；
* multiple module resolutions。

### 7.3 方法

对每个候选 barrier module 计算：

1. bootstrap stability；
2. module assignment flip rate；
3. condition number / collinearity index；
4. cross-resolution consistency；
5. deconvolution consistency；
6. carrier cell-state consistency；
7. perturbation / spatial mappability。

通过 gate 的 modules 才能进入 dominant barrier 判定。未通过 gate 的 modules 不硬拆，而是合并为 coarse barrier class。

### 7.4 输出

* barrier identifiability table；
* dominant barrier eligibility；
* merged coarse barrier definitions；
* unstable barrier log。

### 7.5 生物学作用

该 gate 把“我们假设 myeloid/Treg/stromal 可分”改为“我们量化哪些 barrier 可分，哪些只能合并”。它直接保护 H2 和 H3 的可信度。

---

## 8. Structured Repair Backbone, SRB

### 8.1 输入

SRB 输入 patient-timepoint features：

* module score vector；
* cell-state fractions；
* pathway/TF/signature activities；
* treatment context；
* cancer/cohort environment；
* response-label environment；
* timepoint；
* optional TCR features；
* X-class descriptors, for alignment and transport condition；
* perturbation coverage signals。

### 8.2 Module-factorized encoder

SRB encoder 将输入映射成按 module 分块的 latent：

`Z = {z_1, z_2, ..., z_M}`

每个 `z_m` 对应一个 immune module。每个 module 可以拥有自己的 local encoder，也可以使用共享 MLP 加 module embedding 实现半共享编码。

生物学目标：

* 将患者状态拆解为多个可解释 immune barrier / sensitivity modules；
* 避免黑箱 global embedding；
* 为 dominant barrier、repair transport 和 X-class alignment 提供明确接口。

### 8.3 Cross-module graph attention

SRB 在 module latent blocks 之间加入 cross-module graph attention。注意力邻接由两部分组成：

1. learnable adjacency；
2. prior adjacency，包括 perturbation prior、TF/pathway prior、LR prior 和 module co-activity prior。

该图注意力用于学习 module interaction proxy。它不是自动的因果边。只有当某条边同时通过 environment stability、ablation necessity、perturbation consistency、paired displacement direction 或 spatial support，才可被描述为 causal-supporting edge。

生物学目标：

* 识别 shared anti-PD1 failure modules 之间的组织方式；
* 识别 context-modulated residual interactions；
* 为 repair transport 限定局部传播路径；
* 为 evidence card 提供可解释结构支持。

### 8.4 SRB heads

SRB 包含五个主要 heads。

#### 8.4.1 Reconstruction head

输入：`Z`
输出：重建 module/pathway/TF/cell-state features。

目的：

* 保证 latent 保留 immune-state 信息；
* 防止模型只为 response label 过拟合；
* 支持自监督训练。

#### 8.4.2 Context head

输入：`Z` + environment variables
输出：shared effect、context modulation、HCC index-context residual effect。

目的：

* 检验 H1：哪些 modules 是 shared anti-PD1 failure/sensitivity backbone；
* 检验 H2：哪些 shared modules 在 HCC 或其他 context 中被改写；
* 将 response label source、cohort、cancer、timepoint 作为 environment 处理，而不是当作普通生物特征；
* （B1）HCC etiology 可作为 HCC-context 限定的**低权重、可插拔** environment 协变量，仅用于 HCC residual 的 sensitivity 分析，**不进入 pan-cancer shared module 判定**。

#### 8.4.3 Anchor head

输入：`Z` from clean / expanded / support anchor
输出：responder-compatible direction、anti-PD1 failure direction、anchor reliability。

目的：

* 定义 anti-PD1 responder-compatible state；
* 判断 response/failure modules；
* 为 transport head 提供目标方向。

#### 8.4.4 Transport head

输入：`Z` + dominant barrier + X-class descriptor + graph attention + perturbation coverage
输出：repair displacement `Δ` 和 counterfactual state `Z' = Z + Δ`。

目的：

* 检验 H3；
* 估计修复某 residual barrier 是否将 patient state 推向 responder-compatible direction；
* 形成 X-class-conditioned repair hypothesis。

#### 8.4.5 Alignment head

输入：barrier latent blocks + X-class descriptors + perturbation/spatial/external attributes
输出：X-class compatibility、zero-shot ranking、abstention probability。

目的：

* 检验 H5；
* 支持未见 X-class 和未见 PD1+X 组合的 compositional extrapolation；
* 不替代 repair、spatial 或 perturbation 证据。

### 8.5 Perturbation coverage gate

perturbation prior 不作为真值，只作为 soft condition。SRB 中设置 perturbation coverage gate：

* perturbation coverage 高时，perturbation-consistency loss 权重上升；
* coverage 低时，该 loss 权重下降；
* 无覆盖时不允许模型编造强 perturbation support。

该 gate 保护 X-class repair hypothesis 不被 LINCS / perturb-seq 的 coverage bias 误导。

### 8.6 训练制度

SRB 采用分阶段训练制度，而不是所有 heads 同时乱训。

第一阶段训练 module-factorized encoder、reconstruction head 和 graph attention regularization。目标是得到稳定的 module-structured latent。

第二阶段在低学习率或部分冻结 encoder 的基础上训练 context head、anchor head、transport head 和 alignment head。目标是把稀缺监督信号接入已经稳定的表征，而不是让监督噪声破坏 latent。

多任务 loss 使用 uncertainty weighting 或 GradNorm 一类自动加权策略，避免人工手调多个 loss 系数导致验证集过拟合。

### 8.7 损失函数

SRB 总损失包括：

* reconstruction loss；
* context invariance / modulation loss；
* anchor response loss；
* transport distribution alignment loss；
* paired displacement loss；
* responder-direction loss；
* graph locality regularization；
* perturbation-consistency loss gated by coverage；
* X-class contrastive / ranking loss；
* abstention / calibration loss；
* sparsity and stability regularization。

### 8.8 SRB 输出

* module-factorized patient latent；
* shared anti-PD1 module evidence；
* context-modulated module evidence；
* HCC index-context residual barrier evidence；
* responder-compatible direction；
* dominant barrier candidates；
* repair displacement；
* X-class compatibility ranking；
* zero-shot extrapolation results；
* uncertainty and abstention outputs。

### 8.9 可靠降级

如果 SRB 全模型不稳定，则降级为：

* module score matrix + mixed-effect / meta-analysis；
* responder centroid direction；
* entropic OT / linear displacement repair；
* nearest-prototype X-class matching。

降级后不写 SRB full model claim，但仍可完成 robust immune barrier mechanism analysis。

### 8.10 并行比较

SRB 必须与以下模型比较：

* simple module matrix + logistic / elastic net；
* scVI/patient latent + classifier；
* PCA/NMF patient embedding；
* non-neural OT repair；
* nearest-prototype / attribute matching ZSL。

比较指标包括 response direction stability、paired displacement prediction、context modulation reproducibility、zero-shot ranking、calibration、abstention、confounding leakage 和 external validation consistency。

---

## 9. Counterfactual Repair Transport

### 9.1 生物学目标

Transport head 回答：

**如果修复某个 dominant residual barrier，patient immune state 是否朝 anti-PD1 responder-compatible state 移动？**

这不是临床疗效预测，也不是证明药物真实效果，而是 model-based counterfactual repair hypothesis。

### 9.2 三档实现

#### 高档：conditional flow matching on module latent

适用于 paired pre/post 数据充分、timepoint 清楚、held-out displacement 稳定的情形。模型学习连续位移场：

`dZ/dt = fθ(Z, condition, t)`

该方案可产生 repair trajectory，并自然连接 pre/post 状态迁移。

#### 中档：single-step residual MLP transport

适用于 paired 数据中等但方向可验证的情形。模型学习：

`Tθ(Z, c) = Z + Δθ(Z, c)`

这是默认稳健实现。

#### 低档：entropic OT + barycentric projection

适用于 paired 数据稀缺或 neural transport 不稳定的情形。该方案作为可靠降级和关键 baseline。

### 9.3 CRTN 可信度阶梯

#### L1：方向可信

transport direction 在 held-out patient/module displacement 中显著优于 permutation null、linear displacement、linear-OT 或 signature-reversal baseline。

#### L2：组合迁移可信

在真实 PD1+X 队列或 pre-specified held-out combination validation 中，模型只利用 seen mechanism / single-agent / partial-combination evidence 学到的 repair direction，能够预测组合治疗后的 observed module shift。

若数据时间线不能支持真正前瞻，则只能称 retrospective held-out validation，不能写 prospective。

#### L3：实验因果可信

只有湿实验或功能验证支持时，才能从 model-based counterfactual hypothesis 升级为 stronger causal therapeutic hypothesis。

### 9.4 输出

* repair direction；
* repair displacement；
* barrier-specific repair effect；
* X-class-conditioned repair effect；
* repair uncertainty；
* comparison with OT / linear / permutation baseline；
* CRTN credibility level。

---

## 10. X-class alignment and zero-shot extrapolation

### 10.1 生物学目标

Alignment head 回答：

**一个未见 X-class 或未见 PD1+X 组合，是否与某个 residual barrier 的修复方向兼容？**

### 10.2 X-class descriptor

每个 X-class 由结构化属性定义：

* target family；
* pathway effect；
* expected repair direction；
* carrier cell-state；
* perturbation signature；
* spatial niche cue；
* literature / ontology descriptor；
* known PD1-combination context；
* coverage confidence。

### 10.3 模型

Alignment head 使用 barrier latent block 与 X-class descriptor 做双塔或三塔对齐：

* barrier encoder 取 SRB 中通过 gate 的 dominant barrier `z_m`；
* X-class encoder 编码机制描述；
* evidence encoder 编码 perturbation / spatial / external evidence；
* scoring head 输出 compatibility；
* abstention head 判断是否超出已知机制空间。

训练目标：

* pairwise ranking loss；
* contrastive barrier-X matching；
* leave-one-X-class objective；
* leave-one-combination objective；
* calibration loss；
* abstention loss；
* evidence conflict penalty。

### 10.4 Zero-shot baselines

比较对象不是普通 supervised classifier，而是 zero-shot-capable baselines：

* nearest mechanism prototype；
* attribute matching；
* additive PD1 + X model；
* ontology/literature prior；
* nonparametric retrieval。

### 10.5 输出

* unseen X-class ranking；
* unseen PD1+X ranking；
* mechanism compatibility；
* calibration；
* abstention；
* explanation attributes；
* comparison against zero-shot baselines。

### 10.6 可靠降级

如果 alignment head 不超过 zero-shot baselines，或 calibration / abstention 失败，则降级为 prototype ZSL 和 evidence-card-only ranking。此时 ZSL 仍可作为外推边界分析，但不能作为强方法学主结果。

---

## 11. Spatial Niche Graph Model, SNGM

### 11.1 目标

SNGM 判断 SRB 识别的 dominant barrier 是否具有组织生态基础。

它回答：

* barrier carrier cells 是否在正确组织区域出现；
* 是否形成 myeloid-tumor、T cell exclusion、APC-T cell、CAF/vascular 等 niche；
* ligand-receptor 是否具有空间邻近基础；
* spatial niche 是否支持或反驳 SRB repair hypothesis。

### 11.2 输入

* spatial coordinates；
* spot/cell gene expression；
* cell type / cell-state deconvolution；
* module activity；
* tissue region annotation；
* tumor/stroma/interface labels；
* ligand-receptor candidates；
* optional response / treatment labels。

### 11.3 主方法

先建立 statistical niche baseline：

* neighborhood composition；
* tumor-immune distance；
* co-localization；
* ligand-receptor proximity；
* permutation null；
* region enrichment。

若空间数据质量足够，训练 hetero-GNN：

* nodes：spot / cell / region；
* edges：spatial adjacency、histology adjacency、ligand-receptor potential；
* node features：cell-state probability、module activity、region labels；
* output：niche embedding、niche type、module-niche compatibility、resistance niche probability。

训练目标：

* reconstruct neighborhood composition；
* predict tissue region；
* predict module activity from neighborhood；
* align scRNA-derived module with spatial niche representation；
* optional response association, only when labels reliable。

### 11.4 与 SRB 的关系

SNGM 不与 SRB 端到端合并。二者只在结论层对接：

* SRB 给出 dominant barrier；
* SNGM 检验该 barrier 是否形成 spatial niche；
* Evidence Card 记录支持、缺失或冲突。

### 11.5 输出

* spatial niche graph；
* niche embeddings；
* myeloid-tumor niche score；
* T cell exclusion score；
* APC-T cell niche score；
* CAF/vascular barrier score；
* ligand-receptor spatial support；
* spatial support / conflict verdict。

### 11.6 可靠降级

如果空间数据不足以训练 GNN，则使用 statistical niche adjudication。若空间 metadata 不足以支持 response/failure 判断，则空间证据只能作为组织定位支持，不升强机制 claim。

### 11.7 并行比较

比较 statistical niche model、hetero-GNN、deconvolution-only model 和 scRNA-spatial contrastive alignment。关键指标包括 niche reproducibility、permutation significance、module-niche coherence、response/failure association 和与 SRB repair direction 的一致性。

---

## 12. Evidence Card

### 12.1 目标

Evidence Card 将所有模型输出转化为可审计科学结论。

每个 candidate module / X-class / target 都必须回答：

1. 它是什么 biological barrier；
2. 它是否与 anti-PD1 failure 相关；
3. 它是否 shared 或 context-modulated；
4. 它是否通过 barrier identifiability gate；
5. SRB transport 是否支持 repair；
6. perturbation coverage 是否足够；
7. X-class alignment 是否支持；
8. spatial niche 是否支持；
9. external validation 是否同向；
10. 是否存在关键冲突；
11. 允许写到什么 claim level。

### 12.2 证据层级

#### L1：Robust immune barrier

模块与 anti-PD1 failure / sensitivity 稳定相关，且通过 baseline 和混杂审计。

#### L2：Mechanistic immune barrier

模块具有 pathway / TF / cell-state / graph attention / perturbation / spatial 中至少一种机制支持。

#### L2.5：Repair-compatible barrier

SRB transport 和 X-class alignment 支持某类 repair direction，且没有关键冲突。

#### L3：Experimentally testable PD1+X mechanism hypothesis

候选获得至少两条正交证据支持，包括 repair、spatial、perturbation、external validation 或 tissue feasibility，并且具备明确实验验证路径。

### 12.3 冲突优先规则

关键冲突优先于总分。以下情况必须降级：

* repair direction 与 spatial niche 相反；
* perturbation direction 与 SRB transport 相反；
* external validation 方向相反；
* ZSL 高分但 evidence attributes 不支持；
* signal 主要来自 cohort/cancer/batch/cell fraction；
* barrier 未通过 identifiability gate；
* response label environment 解释了主要效应。

---

## 13. 验证体系

### 13.1 SRB 验证

* reconstruction quality；
* module latent interpretability；
* cohort/cancer leakage；
* context modulation stability；
* anchor direction consistency；
* graph attention ablation；
* transport head ablation；
* alignment head ablation；
* uncertainty calibration。

### 13.2 Barrier gate 验证

* bootstrap stability；
* flip rate；
* condition number；
* deconvolution consistency；
* multi-resolution consistency；
* coarse barrier merge justification。

### 13.3 CRTN 验证

* held-out displacement；
* responder-direction movement；
* comparison with linear / OT / permutation / signature reversal；
* graph locality；
* perturbation consistency；
* spatial consistency；
* credibility level assignment。

### 13.4 ZSL 验证

* leave-one-X-class；
* leave-one-combination；
* calibration；
* abstention；
* comparison with nearest prototype / attribute matching / additive baselines；
* semantic leakage audit；
* evidence consistency。

### 13.5 SNGM 验证

* niche reproducibility；
* permutation null；
* region prediction；
* module-niche coherence；
* response/failure association where labels allow；
* consistency with SRB dominant barrier。

### 13.6 外部验证

* external bulk direction；
* HCC index-context validation；
* treatment-context sensitivity；
* response-label sensitivity；
* failure-case reporting。

---

## 14. 可独立成立的结果层级

本项目的结果可以按证据完整度形成三个层级。它们不是时间表，而是不同数据与模型强度下的科学出口。

### Result Layer A：anti-PD1 failure module atlas and context modulation

核心内容：

* shared anti-PD1 failure modules；
* context-modulated modules；
* HCC index-context residual barrier；
* barrier identifiability map；
* anchor reliability；
* baseline and robustness validation。

该层不依赖强反事实和 ZSL 结果。

### Result Layer B：counterfactual barrier repair and spatial niche support

核心内容：

* SRB transport head；
* dominant barrier repair direction；
* SNGM spatial niche support；
* perturbation coverage gate；
* CRTN credibility level；
* repair-compatible mechanism hypotheses。

该层形成机制修复主故事。

### Result Layer C：X-class zero-shot extrapolation and PD1+X hypothesis nomination

核心内容：

* X-class descriptor ontology；
* alignment head / compositional ZSL；
* unseen X-class / unseen combination ranking；
* evidence cards；
* experimentally testable PD1+X hypotheses。

该层是转化外推最强的出口。

---

## 15. 允许与禁止的 claim

### 15.1 允许的强科学 claim

在证据充分时，允许写：

* identify shared anti-PD1 failure modules；
* identify context-modulated immune barriers；
* define HCC as a pre-specified organ-context deep dive；
* learn a module-factorized structured repair backbone；
* estimate model-based counterfactual repair trajectories in module space；
* identify spatial niches supporting immune failure mechanisms；
* prioritize X-class repair hypotheses；
* demonstrate zero-shot extrapolation to unseen X-classes or combinations；
* nominate experimentally testable PD1+X mechanism hypotheses。

### 15.2 禁止的 claim

禁止写：

* clinical treatment recommendation；
* individual patient treatment prediction；
* proven therapeutic efficacy；
* drug X will improve anti-PD1 response in patients；
* graph attention equals causal edge；
* ZSL alone discovers new biology；
* LINCS / scFM / perturbation prior alone proves mechanism；
* spatial proximity alone proves functional interaction；
* HCC represents all cancers；
* model-based counterfactual equals clinical benefit duration；
* MVP-derived heuristic counterfactual is formal evidence。

---

## 16. 本方案的核心逻辑

本研究用 SRB 实现生物学目标的方式如下：

1. anti-PD1 failure 被表示为 module-factorized latent state，而不是黑箱 embedding；
2. shared failure modules 由 context head 和 environment stability 识别；
3. organ-context modulation 由 context residual effects 识别，HCC 是预指定 index context；
4. dominant barrier 先经过 identifiability gate，确保可辨识后才进入 repair；
5. transport head 学习 barrier repair 如何把状态推向 responder-compatible direction；
6. perturbation prior 通过 coverage gate 影响 repair，不作为真值；
7. alignment head 检验 X-class 和 PD1+X 组合是否能在机制空间中外推；
8. SNGM 独立验证 dominant barrier 是否形成 spatial niche；
9. Evidence Card 将所有证据转化为可审计、可降级、可实验验证的科学 hypothesis。

这条链条的科学主线是：

**anti-PD1 failure modules → context modulation → identifiable residual barrier → model-based repair transport → spatial niche validation → X-class zero-shot extrapolation → experimentally testable PD1+X mechanism hypothesis。**
