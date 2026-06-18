# 科学方案 v7

## 以 PD-1 为锚的空间屏障结构与 X-class 修复框架

版本：v7.0
日期：2026-08-31
状态：当前活跃科学方案

本文只回答“我们要知道什么、为什么重要、什么证据足以支持结论”。具体阶段、工程任务和执行顺序见同目录下的 `roadmap.md`。

v7 不降低 v6 的宏观目标。它继续研究：

- 跨癌种可复用的免疫治疗失败机制；
- 器官和肿瘤环境如何改写这些机制；
- PD1+X 是否可能修复 PD-1 释放之后仍存在的免疫屏障；
- 如何从机制空间外推到新的 X-class，并在证据不足时拒绝作答。

v7 的结构性改变是：空间信息不再是末端验证 sidecar，而是定义屏障、区分共享与情境特异机制、判断 repair 是否具有组织学意义的主线证据。v6 已形成的数据接口、response-blind 分子程序和测量审计继续复用；SRB、FM07-only、旧 Phase8 与空间末端裁判路线不再具有默认优先权。

---

## 1. 核心科学命题

### 1.1 一句话命题

**PD-1 治疗失败不是某个基因或免疫程序单独升高，而是多个分子程序和细胞状态被组织成可重复的空间屏障结构；不同癌种可以共享功能结构，但由不同细胞、尺度和组织背景实现；有效的 X 应当重排、瓦解或绕过这些结构，使组织生态向 responder-compatible 状态移动。**

### 1.2 总体证据链

```text
response-blind 分子程序与细胞状态
  → 空间屏障结构的开放发现
  → 共享功能结构与器官特异实现
  → PD-1 response / failure 临床锚定
  → PD1+X 相关的空间与纵向重排
  → 扰动、靶点、bulk 和独立队列证据
  → 可实验检验的 X-class repair hypothesis
```

这里的“开放发现”表示：不预设所有空间结构都必须是 TLS、肿瘤—基质边界、血管或坏死区。已知结构只验证部分发现；稳定但暂无名称的连续空间场可以保留为主要发现对象。

### 1.3 项目不是什么

本项目不是：

- 把几百万个细胞当成几百万名患者的疗效预测；
- 把不同癌种、药物、终点和时间点混成一个泛癌 AUC；
- 把图注意力、空间邻近或表达相关性直接称为因果关系；
- 用复杂网络生成的位移自动宣称 X 修复了 PD-1 failure；
- 仅因空间数据很多而把所有模态端到端拼接；
- 临床用药推荐或患者级决策系统。

---

## 2. 为什么重要，以及 v7 的创新落点

### 2.1 科学重要性

PD-1/PD-L1 阻断只解除一类 checkpoint 抑制。患者仍可能因抗原呈递不足、髓系抑制、CAF/ECM、血管异常、T 细胞排斥、TLS 缺失、终末耗竭或器官耐受背景而失败。这些问题不仅是“某类细胞多少”，更取决于细胞位于哪里、与谁接触、是否跨越边界、是否形成支持或排斥回路。

因此，只有把分子程序、细胞载体和组织结构连接起来，才能回答：

1. 哪些 barrier 是跨癌种可复用的；
2. 哪些只是器官或癌种的具体实现；
3. 哪些 X-class 在组织层面具有合理的修复方向；
4. 哪些看似有效的分子变化并没有改变真正的空间屏障。

### 2.2 当前领域已经做到什么

近年的工作已经证明：

- 空间细胞组成、邻域和基因特征可以与 ICI 结局相关，并能在独立 NSCLC 队列复现；
- 肿瘤—基质边界、DC/T-cell 相互作用与 CAF 结构可关联 CRC 的 anti-PD1 response；
- 跨组织空间 hub、泛癌 recurrent niche 和空间多组学预测已经可行。

代表性依据包括：

- [Aung et al., Nature Genetics, 2025：234 名 NSCLC 患者的空间多组学与免疫治疗结局](https://www.nature.com/articles/s41588-025-02351-7)
- [Feng et al., Nature Communications, 2024：CRC 肿瘤—基质边界与 anti-PD1 response](https://www.nature.com/articles/s41467-024-54710-3)
- [Dhillon et al., Nature Biotechnology, 2025：跨组织空间 hub](https://www.nature.com/articles/s41587-024-02173-8)
- [Pan-cancer spatial TME study：373 个样本、12 个癌种的 recurrent niches](https://pubmed.ncbi.nlm.nih.gov/41999752/)

因此，v7 的创新不能只是“做一个空间 atlas”“做 neighborhood clustering”或“用 GNN 预测 response”。

### 2.3 v7 的核心创新

#### 创新一：空间屏障结构是机制单位

每个候选结构同时描述：

- 分子程序是什么；
- 哪些 cell states 承载；
- 位于肿瘤、基质、边界还是免疫聚集区；
- 形成邻近、排斥、穿透、梯度、连续场还是多峰结构；
- 作用尺度和不确定性；
- 在治疗后是否被重排。

这比单纯 signature 或 cell fraction 更接近可干预的组织机制。

#### 创新二：共享功能结构与情境实现分离

跨癌共性不要求每个癌种出现完全相同的细胞标签或几何形状。可以共享的是“髓系—肿瘤抑制界面”“血管/CAF 排斥”“APC–T-cell 支持缺失”等功能结构，而 carrier、空间尺度和器官背景允许不同。

#### 创新三：开放的潜在空间场与已知结构共同进入主线

TLS、肿瘤—基质边界等已知结构是强验证锚点，但不定义发现边界。无明确名称、没有肉眼边界但跨患者可复现的连续场，可以作为候选 barrier architecture；只有得到独立证据后才命名。

#### 创新四：repair 被定义为观测到的组织重排

repair 不是模型画出的抽象 transport。它要求至少看到以下证据中的组合：

- 患者内 pre/on/post 状态移动；
- 空间屏障的瓦解、穿透、重组或支持性 niche 建立；
- X 靶点在合理 carrier cells 中可达；
- 真实扰动方向与预期一致；
- 独立组合队列或外部 bulk/scRNA 方向一致。

#### 创新五：跨模态通过共同生物学对象连接

scRNA、spatial、bulk 和 perturbation 不做原始特征硬拼接，而通过共同的 module、cell state、spatial architecture 和 patient-level evidence 连接。每个 bridge 都显式报告覆盖率、映射置信度和 abstention。

#### 创新六：复杂模型必须证明必要性

透明的距离、邻接、边界、空间统计和分层模型是正式基线。潜在场模型、图表示和 GNN 只有在整患者/整队列留出中提供稳定增量时才升级。模型复杂度不是项目创新性的来源。

---

## 3. 核心概念

### 3.1 分子程序

在不使用 response 标签的条件下，由共表达、已知机制基因集或细胞状态特异表达定义的功能坐标。v6 的 8 个 FM 保留为经验基线，但不是固定真理，也不是已经成立的 failure barriers。

### 3.2 空间屏障结构

分子程序、cell-state carrier 与空间组织共同形成的可重复对象。它可以表现为：

- 连续潜在空间场；
- 肿瘤—基质或肿瘤—免疫边界；
- 细胞邻域或生态位；
- 局部邻近/排斥回路；
- 多尺度梯度；
- TLS 等已知结构；
- 暂时无法命名但可复现的组织场。

### 3.3 空间回路

空间屏障结构中能够用 carrier、邻接/排斥关系和功能方向描述的部分。例如 myeloid–tumor suppressive interface、vascular/CAF exclusion 或 APC–T-cell support niche。空间回路是统计和机制假说，不自动等同因果通路。

### 3.4 共享结构与情境实现

- **共享结构**：在多个癌种或队列中保留相同功能拓扑和方向。
- **情境实现**：同一功能由不同 carrier、区域、尺度或器官背景实现。
- **情境特异结构**：仅在一个器官/癌种内稳定，且无法归入共享功能结构。

### 3.5 Responder-compatible ecology

在可比较的 PD-1/PD-L1 临床环境中，与响应稳定相关的多模块、细胞和空间组织状态。它不是单一 responder centroid，也不要求所有癌种完全相同。

### 3.6 Residual barrier

在 checkpoint release 之外仍与 failure 兼容，并具有分子、细胞载体和组织结构证据的障碍。允许单一、共主导、粗分组、混合、弥散或 unresolved。

### 3.7 Repair-compatible rewiring

PD1+X 后出现的观测状态变化，其方向与削弱 residual barrier、增强 responder-compatible ecology 相容。没有可交换对照时，只能称为 repair-compatible 或 combination-associated rewiring，不能称为 X 的因果增量效应。

### 3.8 X-class

按主要修复机制定义，而不是按药名定义，例如 vascular/stromal remodeling、myeloid reprogramming、APC/IFN restoration、epigenetic priming 或 regulatory suppression relief。

---

## 4. 可检验假设

### H1：跨癌种存在可复用的 response-blind 分子程序和细胞状态

这些程序在 patient-balanced、leave-cohort-out 条件下保持可测量性和相似功能方向，并优于随机基因集或纯队列特征。

**削弱条件：** 稳定性主要由单一队列、组织来源、ambient RNA 或粗细胞比例驱动。

### H2：空间组织包含超越 abundance 的增量信息

在控制 cell fraction、module abundance、区域面积、平台与 QC 后，邻接、排斥、边界穿透、空间尺度或潜在场仍能解释结构锚点、response 或跨队列重复性。

**削弱条件：** topology-aware 模型在整患者留出中不优于 abundance-only、region-only 或空间平滑基线。

### H3：未知或未命名空间场可以稳定发现

不依赖已知结构标签、定义性 marker 或形态边界，至少部分连续空间场仍在患者、组织块、切片或平台间具有可重复子空间；TLS 和边界只验证其中一部分。

**削弱条件：** 潜在场只在单 seed、单 fold、单平台或直接 marker 可见时成立，或空间置换后不再有增量。

### H4：跨癌共享的是功能结构，器官环境改变其实现

至少一类 barrier architecture 在多个癌种中保持功能拓扑，同时 carrier、空间尺度或区域分布出现可解释的 context modulation。

**削弱条件：** 所有结构只能在单癌种复现，或所谓共享完全由平台和 cell-type annotation 映射造成。

### H5：HCC 对共享结构产生可复现的肝脏情境改写

HCC 中的肝脏耐受、髓系、血管、CAF/ECM、缺氧与局部治疗背景，应改变共享 barrier architecture 的实现，而不是简单增加某类细胞。

**削弱条件：** HCC residual 被组织纯度、细胞组成、平台、取材位置或单队列效应完全解释。

### H6：PD-1 response/failure 与多结构生态状态相关

在 treatment×endpoint×cancer 环境内，responder-compatible ecology 和 residual barriers 应在患者级复制；跨环境可以共享、情境调节、特异或 unresolved。

**削弱条件：** cohort/treatment/endpoint-only 基线持续强于生物和空间特征，或独立队列方向反转。

### H7：有效 X-class 与 residual barrier 的观测重排相容

PD1+X responder 或 paired pre/on/post 样本中，应出现与 barrier 削弱方向相容的分子重排；空间重排必须由同患者纵向空间数据直接观测，或由不使用 response、并在整患者留出中验证过的 architecture surrogate 间接支持，同时得到 perturbation/target 或外部队列支持。

**削弱条件：** 位移在配对/时间置换后消失、仅由 composition/QC 解释、独立组合队列方向冲突，或靶点不在合理 carrier 中表达。

### H8：结构化 repair space 支持跨 X-class 外推

当存在至少三类机制不同的 X-class、明确正负例和独立组合数据时，基于 barrier–repair compatibility 的 class-held-out 排序应超过药名相似、nearest prototype 和简单 pathway overlap，并在证据不足时 abstain。

**削弱条件：** 性能依赖药名或标签泄漏，leave-one-X-class 不超过简单基线，或只能在单一 antiangiogenic class 内成立。

---

## 5. 证据与模型架构

### 5.1 分子与细胞状态层

输入包括：

- v6 response-blind empirical programs；
- cell-state fractions 与 cell-state-specific pseudobulk；
- VEGF/angiogenesis/endothelial、CAF/ECM、hypoxia、suppressive myeloid、APC/IFN、cytotoxic、progenitor/terminal exhaustion、TLS/B、NET、WNT exclusion 等机制轴；
- TCR、ligand–receptor 或肿瘤内在特征仅在覆盖和语义合格时进入。

该层回答“什么功能存在、由谁承载”，不直接回答“在组织中如何形成 barrier”。

### 5.2 空间主干

空间主干包含三种互补表示。

#### A. 可解释空间统计

- 局部组成；
- tumor–immune distance；
- boundary penetration；
- myeloid–tumor、APC–T-cell、CAF/vascular、TLS/B 等邻接；
- Moran/variogram、空间尺度、区域 enrichment；
- ligand–receptor proximity 与空间置换 null。

#### B. 开放潜在空间场

在 response-blind、结构标签隔离条件下发现连续场和可识别子空间。模型允许 `K_eff=0`、多个可识别场或 unresolved，不强制把场命名为已知结构。

#### C. 图表示候选

每个患者/切片独立成图，节点可以是 cell、spot、bin 或局部区域，边只表达物理邻接或版本化的功能关系。图模型必须超过可解释统计基线，否则删除，不影响科学主线。

### 5.3 跨模态桥

统一对象是：

```text
gene/module ↔ cell state ↔ spatial architecture ↔ patient-level evidence
```

每个 bridge 必须报告：

- gene coverage 与缺失；
- label-transfer/reference 来源；
- resolution 和 coordinate scale；
- mapping confidence；
- platform sensitivity；
- unmapped 与 abstain；
- 是否使用图像、分割或 ground truth 派生变量。

### 5.4 临床与 context 层

先在每个 treatment×endpoint×cancer 环境内估计，再进行 meta-analysis 或 hierarchical modeling。允许输出：

- shared；
- context-modulated；
- context-specific；
- heterogeneous/conflicting；
- unresolved/abstain。

不同 endpoint、治疗方案和时间点不原始池化。患者是临床独立单位，spot/cell 不能扩大临床样本量。

### 5.5 Repair 与 X-class 层

只有同时具有测量、空间组织和 PD-1 failure 支持的结构，才可进入 repair 分析。优先使用真实 paired displacement、低容量分层模型和透明 transport baseline；复杂反事实模型只有在可交换比较和外部验证存在时才升级。

空间 repair 必须区分三种证据：同患者纵向空间数据支持 `observed_spatial_rewiring`；经整患者留出验证的非空间代理只支持 `surrogate_inferred_architecture_rewiring`；两者都没有时，`architecture_rewiring=NOT_IDENTIFIABLE`。横断面 post-treatment 空间样本只能提供结构兼容性支持，不能填补纵向空间证据。

### 5.6 Evidence Card

每个机制候选单独记录：

- molecular support；
- carrier-state support；
- spatial topology/field support；
- PD-1 response association；
- cross-cohort/cancer replication；
- HCC context modulation；
- paired/combination displacement；
- perturbation/target support；
- bulk/external support；
- negative controls；
- conflicts；
- permitted claim 与 abstention reason。

不得用总分掩盖方向冲突。

---

## 6. 数据角色

### 6.1 `/006` 的主资产

- 44 个 count-compatible 单细胞对象：泛癌 response-blind 程序、cell-state carrier 与跨队列稳定性；
- TASK01：HCC PD-1 临床锚与患者内纵向轨迹；
- TASK02：HCC PD1+X 组合相关纵向状态；
- GSE286827、GSE301741 等：分层外部 anchor 或待修复候选；
- GSE238264：HCC post-treatment 空间支持；Mendeley 的 11 个 HCC spatial object 需先闭合 patient/response/scale provenance；GSE291246 是 BCC Xenium，不进入 HCC 证据；
- bulk 队列：患者级外部方向和临床尺度验证；
- GSE133344、GSE193736、GSE306429、GSE90063、XAtlas、L1000：真实 perturbation/target direction；
- TCR 与 ligand–receptor 资产：覆盖合格时提供辅助机制证据。

### 6.2 `~/013_spatial` 的主资产

当前已审计的高价值层包括：

- HTAN CRC：30 名患者、31 个组织块、47 张训练 section；适合 spatial foundation 与开发；
- GSE175540 ccRCC：24 个样本，含可审计 TLS 锚点；
- ST_CRC_CMS：7 名患者、14 张连续切片、12 个可审计来源，适合边界和患者内跨切片稳定性；
- TLS_VISIUM_USZ：肾癌与肺癌患者级 TLS 外部验证；
- HEST、10x、STOmicsDB 等大规模多平台资产：空间 atlas、跨平台压力测试和候选发现；临床 metadata 不足时不得承担 response 结论；
- 现有 R-04：可复用数据合同、patient-grouped split、ground-truth 隔离和优化审计；最新状态为 `RESTART_STABILITY_DIAGNOSTIC_COMPLETE_NOT_FORMAL_K_SELECTION`，不得继承为已成立空间场或有效 K。

### 6.3 数据角色原则

数据角色按“具体 claim”管理，而不是给整个项目永久贴标签。探索性 discovery/calibration 可以迭代，但一旦某数据参与 ontology、候选结构、阈值或模型选择，它对相关 claim 只能承担 development/support，不能再称独立验证。承担 claim-bearing external validation 的数据必须在该 claim 的所有选择完成前保持隔离；secondary role 只能提供不承担结论升级的支持。`~/013_spatial` 已有 role ledger 中的 GSE175540、ST_CRC_CMS、TLS_VISIUM_USZ 等隔离边界继续生效，其作用是防止选择后验证，不是预注册科学假设。

---

## 7. 独立单位、泄漏与可比性

### 7.1 独立单位

- 临床结论：patient；
- 空间发现：patient 是外层独立单位；block、section、region、spot/cell 是患者内嵌套观测；
- 配对变化：patient-pair；
- 扰动：独立 perturbation-condition/replicate；
- 跨队列结论：cohort/environment。

### 7.2 禁止泄漏

- 同一患者不同时间点、block、section 不得跨训练/验证；
- bootstrap、split、置换和外部样本计数均以 patient 为外层；block/section 只进入患者内层级模型，不增加独立样本数；
- 同一来源的镜像、重复文件、相邻切片不得伪装成独立外部验证；
- ground truth、GT distance、GT-derived mask 和定义性 marker 不得循环进入对应发现模型；
- response 参与特征选择后，必须使用未触碰队列或嵌套评估；
- 不能用 cohort、药物、endpoint 或文件名作为隐式 response 代理。

### 7.3 跨平台边界

Visium、Stereo-seq、Xenium、Visium HD 和 imaging-based 数据的分辨率、坐标和测量误差不同。共享的是患者级空间统计与功能结构，不要求原始 spot/cell embedding 完全同分布。

---

## 8. 结论等级

由弱到强：

1. 可重复分子程序；
2. 可重复空间场/结构；
3. topology 超越 abundance 的空间屏障结构；
4. 跨队列或跨癌共享结构；
5. 与 PD-1 response/failure 稳定相关的结构；
6. HCC 等器官情境改写；
7. PD1+X-associated repair-compatible rewiring；
8. spatial + perturbation + external 支持的 repair mechanism hypothesis；
9. 可实验验证的 X-class nomination；
10. X 的因果增量效应或临床推荐。

第 10 级需要随机化、可交换 treatment arms 或等价强度的因果证据，当前数据不默认支持。

---

## 9. 最低充分发表证据

### 9.1 空间屏障结构主结果

- 至少两个独立患者队列复现同类结构；
- 整患者/整队列留出；
- topology 超过 composition、module abundance、region proportion 与 spatial smoothing 基线；
- carrier、区域、距离/邻接和空间置换证据齐全；
- 至少一个已知结构锚点验证，同时允许未命名场；
- 明确失败和不可识别结构。

### 9.2 泛癌与 HCC context 主结果

- 至少一个结构在多个癌种共享功能拓扑；
- leave-one-cancer/cohort 后方向保持；
- HCC 中的实现由独立 HCC 数据支持；
- 主要结果不能由平台、纯度、区域大小、ambient 或 cell fraction 单独解释。

### 9.3 PD1+X repair 主结果

- 一个可信的 PD-1 failure-compatible barrier architecture；
- 一个真实 PD1+X paired 或 longitudinal displacement；
- 一个未参与发现的组合队列或外部患者层验证；
- spatial 与 perturbation/target 至少一类强支持，另一类不冲突；
- 至少一个机制不匹配或失败组合负例；
- 结论限定为 repair-compatible hypothesis，除非有可交换 treatment arms。

若要在标题或主结论中使用“空间结构重排”，还必须有同患者纵向空间观测，或通过整患者留出、未使用 response 的 architecture surrogate。前者可称观测空间重排，后者必须标为代理推断；两者均无时，空间 repair 保持 `NOT_IDENTIFIABLE`，只能报告分子/细胞状态纵向变化与横断面空间兼容性。

### 9.4 X-class 外推主结果

- 至少三类机制不同的 X-class；
- 正例、失败例和机制不匹配负例；
- class-held-out 评测；
- 超过简单 pathway overlap 和 prototype baseline；
- 有校准的 abstention。

---

## 10. 可接受的否定结果

v7 允许并要求报告：

- 空间 topology 不增加 abundance 之外的信息；
- 只有 context-specific、没有泛癌 shared architecture；
- HCC 差异由细胞组成或取材解释；
- 潜在场不稳定或 `K_eff=0`；
- GNN 不优于空间统计；
- PD1+X 位移不朝 responder-compatible ecology 移动；
- perturbation、spatial 和 clinical 方向冲突；
- 只有一个 X-class，无法支持 zero-shot；
- 证据不足而 abstain。

这些结果不会使项目失去科学价值；它们定义空间机制、可迁移性和 repair 推断的适用边界。

---

## 11. 最终预期叙事

若核心假设得到支持，最终主张为：

> 肿瘤免疫治疗失败由分子程序、细胞载体和组织空间共同形成的屏障结构决定。部分功能结构可跨癌种迁移，但由器官环境以不同 carrier 和尺度实现。PD1+X 的机制价值可通过这些结构是否发生与 responder-compatible ecology 相容的观测重排来评估，并由扰动、外部患者和组织证据约束为可实验检验的 X-class repair hypothesis。

若只有空间发现成立而临床/repair 证据不足，结论收缩为“泛癌空间屏障结构及其器官实现”；不强行声称 PD-1 failure 或 X repair。

---

## 12. 范围边界

- 不把横断面空间相关解释为结构因果生成；
- 不从二维切片恢复唯一真实三维组织；
- 不要求所有潜在场都有现成生物名称；
- 不把一个强癌种或结构替所有癌种背书；
- 不因某个模型失败而删除仍有生物信息的模块；
- 不把 unsupported、post-only 或 metadata 不完整数据升级为临床锚；
- 不使用 v6.1 hash 生成的 perturbation 连续分数作为真实实验效应；
- 不把 v6 历史 Phase8 结果当作 v7 现行结论；
- 不形成临床治疗建议。

v7 的价值不依赖“模型越复杂越好”，而依赖能否把跨癌免疫程序、空间组织、临床响应和可干预修复连接成一条可证伪的证据链。
