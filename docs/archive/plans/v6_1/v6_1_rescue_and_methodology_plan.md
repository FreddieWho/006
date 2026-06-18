# v6.1 修复与方法学创新方案（Rescue & Methodology Innovation Plan）

**生成时间**：2026-06-17
**依据**：`reports/v6_1_science_engineering_alignment_report.md`（二次修订版）+ `docs/plan/v6_1/{plan,roadmap}.md` + `results/v6_1/` 实际进度 + 外部数据精确检索
**性质**：在不推翻 v6.1 主线的前提下，补齐方法学创新、重排技术路线、给出外部数据精确清单与双轨救援门控。

---

## 0. 战略决策回执（用户已定）

| 维度 | 决策 | 对方案的约束 |
|---|---|---|
| 方法学创新方向 | **三方向（SS-CGF / 反事实 PD-1→残障→X修复 / 多模态证据图）整合为同一框架，且三轨可并行**；风险过高时给推荐组合 | 框架须分层解耦、接口清晰，单层失败不连坐 |
| 目标定位 | **生物医学顶刊（Nature / Cancer Cell / Immunity 级）** | 创新点=架构/组合/问题建模创新 + 机制发现 + HCC 转化闭环；不必做 benchmark 冠军，但方法须扎实新颖、可被审稿人认作"新框架"而非"标准流程拼装" |
| PD-1 中心主张 | **双轨并行，按救援结果决定** | 设 anchor 救援门；过则回归 PD-1 中心主叙事，否则降级 HCC-first |
| 资源边界 | **算力充足 / 湿实验可纳入 / 时间紧需快速可发** | 每层设 v0(快速可防守)→v1(创新)→v2(野心)分级；先交付可发版本，创新增量并行推进 |
| 排除项要求 | 任何被排除或可选未采用的方案须给**具体理由** | 见第 9 节，逐条具体论证 |

---

## 1. 设计总原则

1. **不重做地基，只补主链与方法学**。Step1/Step2 数据底座、Phase4–10 治理流水线保留为"工程脚手架"；本方案替换的是科学内核（causal/counterfactual/evidence-graph）与缺口数据。
2. **分层解耦 + 接口契约**。三大创新方向不是三套割裂算法，而是一个机制提名引擎的三层栈（结构层→干预层→裁判层），层间用固定数据契约连接，可并行开发、可单独降级。
3. **MVP 先行、创新并行**。每层都有 v0/v1/v2。v0 用现有可解释方法保证"30 天内有可发结果"；v1/v2 是冲顶刊的方法学增量，在 v0 基础上叠加，失败不影响可发性。
4. **证据驱动的 claim 门控**。所有结论必须经多模态证据图裁定到 L1/L2/L3，并带显式降级规则——这本身既是工程合规，也是方法学卖点。
5. **双轨 PD-1**。anchor 救援与 HCC-first 退路并行铺设，里程碑处一次性裁决主叙事，避免把鸡蛋放在断裂的第一环上。

---

## 2. 统一方法学框架：**PASCAR**

> **PASCAR** = **P**D-1-**A**nchored **S**hared-specific **CA**usal mechanism graph with counterfactual **R**epair scoring, adjudicated by a multimodal evidence graph.
> 一句话：把"以 PD-1 单药生物学为锚的共享/特异**因果机制图**"建出来，在图上定义"修复某机制模块会把患者推向 responder-like 多远"的**反事实修复算子**，再让**多模态证据图**对每条 claim 做 L1/L2/L3 裁判。

这正好把用户要的三个方向焊成一条链，且每层独立可发：

```
┌─────────────────────────────────────────────────────────────────┐
│  输入契约：module_score_matrix (Phase4 已产) + PD-1 anchor labels  │
│            + 治疗上下文 + 多模态层（TCR/spatial/perturb/bulk）       │
└───────────────┬─────────────────────────────────────────────────┘
                ▼
   Layer A —— 结构层：Shared-Specific Causal Graph (SS-CGF)
   "哪些免疫机制模块之间存在跨癌种不变(shared)/HCC特异(specific)的因果/方向关系？"
   产物：module-level causal DAG（shared backbone + HCC-specific residual edges）
                │  接口：因果图 G(shared), G(HCC-specific)
                ▼
   Layer B —— 干预层：Counterfactual PD-1 Repair Operator
   "在 G 上对模块 M 施加 do(repair)，患者状态向 PD-1 responder-like manifold 移动多少？"
   产物：mIMS / MRI / mICS（正式文件）、X mechanism class 排序
                │  接口：candidate (module, X-class, target) + 反事实位移分数
                ▼
   Layer C —— 裁判层：Multimodal Evidence Graph Adjudicator
   "每条 candidate 的 L1/L2/L3 证据是否成立？冲突如何降级？"
   产物：candidate_module_evidence_matrix.csv（L1/L2/L3 + 降级规则）、evidence cards
                ▼
       输出：PD-1+X mechanism nomination + 最小湿实验计划
```

### 2.1 为什么这个组合对生物医学顶刊有竞争力（创新论证）

- **不是单算法先进，而是问题建模 + 架构组合创新**：据我所知，ICI scRNA 领域大量工作止步于"差异模块/signature 关联 response"。PASCAR 的新颖性在于把**因果结构（不变性）× 反事实干预（最优传输位移）× 多模态正交证据裁定**三者首次串成一个"机制→可干预性→证据等级"的端到端提名引擎，并以 **PD-1 单药生物学为因果锚**。这是审稿人能一眼记住的"新框架"。
- **共享/特异 = 因果不变性**：把"跨癌种共享层"形式化为**跨癌种环境不变的因果机制**、"HCC 特异"为**环境特异残差边**——这给了"shared vs specific"一个有原理的定义，而非聚类约定，直接回应 plan §3.2/§4 的科学问题，是方法学硬通货。
- **MRI/mICS = 反事实位移**：用最优传输/反事实把"X 修复了什么"量化为"non-responder→responder-like 流形的位移方向与可达性"，让 plan 里一直停在概念层的 M-signals 变成有定义、可计算、可证伪的量。
- **L1/L2/L3 = 证据图的一等公民**：把证据等级做成算法输出而非事后表格，天然满足顶刊对 robustness/可重复性的审查。

---

## 3. 三层栈的分级实现（v0 / v1 / v2）

每层 v0 必做（保可发），v1 为顶刊创新增量（并行推进），v2 为野心选项（算力足时尝试）。

### Layer A — Shared-Specific Causal Graph

| 级别 | 方法 | 产物 | 风险 |
|---|---|---|---|
| **A.v0** | 模块级 **partial correlation / GGM**（graphical lasso）分别在 shared 池与 HCC 池估计，差集=specific 边；within-cohort meta + 置换零分布 | 关联性模块图（非因果），可发 | 低 |
| **A.v1（推荐创新）** | 模块级 **不变因果发现**：跨癌种作为"环境"，用 ICP / Invariant-Risk 风格筛 shared invariant 因果父集；HCC 特异边用环境交互项识别；**扰动先验作方向软约束**（perturb 方向与边方向冲突则降权） | shared causal backbone + HCC-specific residual DAG | 中 |
| **A.v2** | **可微因果发现**（NOTEARS/DCDI 变体）+ 扰动作为"软干预"数据，模块级 DAG 联合 shared/specific 分解 | 全可微因果图，benchmark 友好 | 中高 |

- **关键设计**：在**模块级**（~数十–数百模块）而非基因级做因果，降维到可解释、可识别、算力可控；这是把 SS-CGF 从 v5.3 的"野心"落地为"可做完"的核心取舍。
- **不变性 = 共享**：这是 A.v1 的方法学灵魂，建议作为主推。

### Layer B — Counterfactual PD-1 Repair Operator

| 级别 | 方法 | 产物 |
|---|---|---|
| **B.v0** | **mIMS** = 患者在模块 M 上的激活/抑制（直接由 module score 标准化）；**MRI** = 用 PD-1 anchor 训练的 responder-like 判别方向上，模块扰动的梯度方向一致性；**mICS** = 规则整合（centrality+perturb+druggability+spatial+external），即 plan §7 的第一阶段版 | 正式 `mIMS/MRI/mICS` 文件 |
| **B.v1（推荐创新）** | **反事实位移**：在 anchor 学到的 responder/non-responder 流形间用**最优传输(OT)**求位移场；对模块 M 施 `do(repair)` 后计算患者沿 OT 方向到 responder-like 的可达性增益=MRI；mICS=可达性×可干预性×证据 | OT 反事实 MRI；X-class 排序 |
| **B.v2** | 基于 Layer A 因果图的 **structural counterfactual**（SCM do-calculus），模块干预沿因果路径传播，预测下游模块与 response 变化 | 因果反事实 MRI（最强但依赖 A.v2 图质量） |

- **双轨耦合点**：B 依赖 PD-1 anchor。anchor 救援成功→B 用 anchor 流形（主叙事）；失败→B 退化为"HCC responder-like / 稳健模块方向"流形（HCC-first），方法不变，只换标签来源。

### Layer C — Multimodal Evidence Graph Adjudicator

| 级别 | 方法 | 产物 |
|---|---|---|
| **C.v0** | 异质证据**加权聚合 + 置换标定**：每模态（scRNA 关联 / TCR / spatial / perturb / bulk）给方向一致性分，按预设权重合成 L1/L2/L3，显式降级规则（任一关键冲突自动降级） | `candidate_module_evidence_matrix.csv` |
| **C.v1（推荐创新）** | **正交证据贝叶斯融合**：把各模态当条件独立证据，贝叶斯更新每个 candidate 为真机制的后验；置换/null 校准每模态似然；输出后验+可信区间作 L3 依据 | 概率化证据卡 |
| **C.v2** | 异质图神经网络（节点=模块/机制/靶点/细胞态，边=证据关系）做 candidate 排序，仅作**敏感性对照**，不作主结论 | GNN 排序（可选佐证） |

- **C 是"快速可发"的最优头条**：即便 A/B 只到 v0/v1，C 把现有 Phase4–7 全部证据重组为顶刊级证据等级体系，本身就是一个可发的"框架贡献"。

---

## 4. 三条可执行技术路线（含推荐与适用条件）

用户要"若干可尝试路线"。三条路线共享 PASCAR 栈，差别在**主叙事重心**与**风险偏好**：

### 路线 R1 —— 全栈 PD-1 中心（推荐主攻，押 anchor 救援成功）
- 配置：A.v1 + B.v1 + C.v1，PD-1 anchor 用 `GSE123813`+`GSE176021`+`GSE200996` 救援。
- 主叙事："PD-1 单药因果敏感性锚 → 共享/HCC 特异残障 → X 修复机制反事实提名 → 多模态证据裁定 + 湿验"。
- 适用：anchor 救援门通过（混杂降至 MEDIUM 以下、within-cohort 方向一致）。
- 卖点：完整回归 plan 原始野心，方法学+生物学双满。

### 路线 R2 —— HCC-first 因果转化（稳健退路，押 HCC 深描）
- 配置：A.v1（HCC 特异层为主）+ B.v0/v1（HCC responder-like 流形）+ C.v1，重仓 spatial 定量裁判。
- 主叙事："HCC 免疫耐受背景下的特异因果屏障（髓系/基质/血管）→ 空间裁判 → PD-1+X(anti-VEGF/TKI/lenvatinib) 转化提名"。
- 适用：anchor 救援失败，或 HCC 信号显著强于泛癌。
- 卖点：转化叙事强、对 anchor 不敏感、风险最低。

### 路线 R3 —— 证据图方法学头条（最快可发，押架构创新）
- 配置：C.v1 为头条方法 + A.v0/v1、B.v0/v1 为组件，把现有 Phase4–7 资产 + 已暂存增补层（TCR/spatial/pseudobulk）直接喂入。
- 主叙事："一个把因果结构、反事实可干预性、多模态正交证据统一为 L1/L2/L3 提名的机制发现框架，在多癌种 ICI + HCC 深描上演示"。
- 适用：时间最紧、要先锁一篇；A/B 创新作为后续升级。
- 卖点：复用资产最多、最快出结果，方法学贡献清晰。

> **推荐**：以 **R1 为目标、R3 为保底、R2 为 anchor 失败时的转向**。三者共栈，里程碑 C（机制合练）处依据 anchor 救援结果与 HCC/泛癌信号强弱一次性选定头条。这满足"三方向整合 + 三轨并行 + 双轨 PD-1"的全部要求。

---

## 5. 外部数据救援精确清单（要求 2）

> 原则：能给编号给编号；受控访问标注申请路径；无法定位的给具体数据要求。所有新数据按 INDEX §9 走 addendum，**不可直接替换 Phase4 主输入**。

### 5.1 干净 PD-1/PD-L1 单药 anchor scRNA（WP1 第一优先，决定 R1 成败）

| 编号 | 队列 | 治疗 | 关键属性 | 访问 | 用途 |
|---|---|---|---|---|---|
| **GSE123813** | Yost 2019 BCC/SCC | anti-PD1 单药 | pre/post, R/NR, **scRNA+scTCR** | GEO 开放 | anchor + TCR + 跨癌种 shared |
| **GSE176021** | Caushi 2021 NSCLC | 新辅助 nivolumab 单药 | MPR R/NR, **scRNA+scTCR**, 560K T | GEO 开放 | anchor + TCR（偏 T 细胞，髓系弱，需注明） |
| **GSE200996** | Luoma 2022 HNSCC | 新辅助 nivolumab（含 nivo+ipi 臂） | 配对 pre/post, **scRNA+scTCR**, 26 pt | GEO 开放 | anchor（**仅取单药臂**）+ TCR |
| Bi 2021 ccRCC | `phs002065.v1.p1` / SCP1288 | ICB(nivolumab) | pre/post ICB, R/NR | **dbGaP 受控**（需申请） | anchor 加强（可选上行） |
| Braun 2020 ccRCC | `EGAS00001004290/91/92` | nivolumab(CM-025) bulk | 大样本 R/NR + 生存 | **EGA 受控** | anchor 外部 bulk（可选） |

- **结论**：仅靠 3 个 GEO 开放队列（BCC+NSCLC+HNSCC 三癌种、均 anti-PD1 单药、均带 R/NR + TCR）即可满足 WP1"2+ 独立干净队列"门，**不依赖受控数据**。受控的 Bi/Braun 作为上行加强。
- **cleanliness 注记**：GSE176021/GSE200996 为新辅助（neoadjuvant）单药——属"较干净 ICI 队列"，符合 plan §3.1 "anchor 不必纯转移性单药"。GSE200996 须剥离 nivo+ipi 联合臂。

### 5.2 HCC ICI 外部锚（WP4 / R2 核心）

| 编号 | 治疗 | 类型 | 用途 |
|---|---|---|---|
| **GSE215011** | nivolumab 单药 HCC | bulk RNA-seq (n≈10) | HCC 单药外锚（小样本，方向验证） |
| **GSE235863** | anti-PD1+lenvatinib, HBV+ HCC | bulk, 配对 pre/post (n≈15) | **PD-1+X(TKI) HCC 锚**，直接服务 X 修复叙事 |
| **GSE140901** | nivo/pembro ± ipi HCC | NanoString 770 (42 pt) | HCC ICI 免疫基因外锚 |
| **GSE279750** | anti-PD-L1 联合 HCC | bulk (n≈10) | HCC 联合外锚 |
| IMbrave150（已有）| atezo+bev | bulk | 现有唯一锚，升级为模块级投影+CI |
| J Hepatol atezo+bev 422 样本 | atezo+bev/atezo/sorafenib | bulk | **疑非全公开**；若可得为最强 HCC 锚（需向作者/补充材料确认） |

### 5.3 泛癌 ICI bulk 外锚（WP4 / R1/R3 跨癌种验证）

| 编号 | 队列 | 治疗 |
|---|---|---|
| **GSE78220** | Hugo 黑色素瘤 | pembrolizumab |
| **GSE91061** | Riaz 黑色素瘤 | nivolumab |
| **PRJEB23709** | Gide 黑色素瘤 | anti-PD1 ± anti-CTLA4 |
| Van Allen 2015 | `phs000452.v2.p1` (dbGaP 受控) | anti-CTLA4（仅作 CTLA4 参考） |
| Liu 2019 黑色素瘤 | dbGaP 聚合（受控） | anti-PD1 |

### 5.4 空间 / TCR / 扰动

| 类别 | 资产 | 状态 / 要求 |
|---|---|---|
| HCC 空间（已暂存）| **GSE238264** | 增补层 manifest(7 样本, support_only)，待定量 |
| HCC scRNA 深描补强 | **GSE149614**（10 例 HCC scRNA） | 可纳入 HCC 特异层 |
| HCC Visium 补强 | **数据要求**：≥1 套 HCC 10x Visium，含核/侵袭前沿/基质分区，最好 ICI 相关 | Wang 2022 *Theranostics* HCC Visium 为候选（编号未公开，需向作者/补充材料索取）；或定向再检索 |
| 黑色素瘤 ICB scRNA+Visium+GRN | **EGAS00001006488**（Pozniak2024, scRNA `EGAD00001009291` / Visium `EGAD00001010921`）| **EGA 受控**；含 GRN 标杆，作 A 层方法对照（可选） |
| TCR（已暂存）| **GSE236581** TCR 克隆表（1.1M 记录, CRC, support_only）+ GSE123813/176021/200996 自带 scTCR | 构建 `tcr_clone_master_v6_1.csv` |
| 扰动先验（本地）| `GSE133344` / `GSE90063`(Dixit Perturb-seq) / `GSE193736` / `GSE306429` / LINCS-L1000 | **需确认本地是否已下载**（见第 11 节）；构建 perturb-prior 矩阵 |

---

## 6. 重构后的工作包与双轨门控

在原 WP0–WP6 基础上**新增方法学 WP（WM-A/B/C）**，并明确双轨里程碑。

| WP | 内容 | 优先级 | 退出门 |
|---|---|---|---|
| **WP0** 合同加固 | 维持 internal-review caveat、PD1_anchor support-only、Phase3.5 invalid | 立即 | 无越界 claim（已基本达成） |
| **WP1** anchor 救援 | 整合 GSE123813+GSE176021+GSE200996（剥离联合臂）→ Step1/2/3 anchor 分支 + within-cohort meta | **P0** | 混杂<HIGH 且 within-cohort 方向一致 → R1；否则 → R2 |
| **WP2** 表达救援回灌 | 把已建 GSE120575/GSE229772 pseudobulk 接入 Step2.8/2.9 → 再生 hotfix；移除幻影患者 `gse120575_2` | P1 | 无零表达主 cohort |
| **WM-A** SS-CGF | 实现 A.v0→A.v1（不变因果，扰动软约束），模块级 shared/HCC-specific DAG | **P0(方法)** | shared backbone 跨癌种稳定 + ≥1 HCC 特异因果屏障 |
| **WM-B** 反事实修复 | 实现 B.v0（正式 mIMS/MRI/mICS 文件）→ B.v1（OT 位移） | **P0(方法)** | M-signals 文件存在且每个 nomination 连回 anchor/HCC 流形 |
| **WM-C** 证据图裁判 | 实现 C.v0→C.v1，产出 `candidate_module_evidence_matrix.csv` + 降级规则 | **P0(方法)** | 每 candidate 有显式 L1/L2/L3 + 冲突降级 |
| **WP3** 定量空间 | GSE238264(+补强)做 spot 级模块投影 / exclusion / 邻近 / LRG 共定位 + 置换零控 | P1 | ≥1 核心机制过定量空间支持 |
| **WP4** 外锚激活 | GSE215011/235863/140901/279750(HCC) + GSE78220/91061/PRJEB23709(泛癌) 模块级投影+CI+失败案例 | P1 | 核心机制 ≥1 独立外锚复现方向 |
| **WP5** 扰动+TCR | 本地扰动→perturb-prior 矩阵；scTCR→`tcr_clone_master_v6_1.csv`；并入证据矩阵 | P2 | 证据矩阵含 perturb/TCR 字段 |
| **WP6** 投稿复现 | Supp.S1 + 参数附录 + 环境锁 + claim audit | P0(投稿前) | external_submission_ready=true |
| **WL** 湿实验 | 最小验证（见第 7 节） | 与 WP4 并行 | ≥1 核心靶点体外方向一致 |

### 双轨里程碑判决
- **Milestone A′（Step1/2/3 anchor 分支后）**：WP1 退出门 → 选 R1 或 R2 头条候选。
- **Milestone C′（WM-A/B/C v1 + WP3/WP4 部分后）**：依 anchor 稳定性 + HCC/泛癌信号强弱，**锁定头条路线**（R1/R2/R3），其余降为补充/敏感性。

---

## 7. 最小湿实验验证设计（湿实验已纳入范围）

目标：给 L3 一个真实闭环，匹配顶刊"至少一个机制有功能证据"的期待。按可达性从低到高：

1. **IHC / mIF 空间共定位验证**（最低成本，先做）：对 1–2 个核心机制（如髓系重编程 `CAND_P5MG_HCC_012`、`CAND_P5MG_PD1X_014`）在 HCC 组织切片验证关键细胞态空间邻近（如 TAM-肿瘤巢、CD8 排斥），直接支撑 Layer C 空间证据从"规则赋标"升为"组织学验证"。
2. **体外共培养 / 极化实验**：验证 X-class 方向（如髓系 M2→M1 重编程、IFN/APC 恢复）对 T 细胞功能的因果方向，校验 Layer B 反事实位移的方向预测。
3. **CRISPR / 小分子扰动**（算力+湿实验充足时）：对 Top-1 靶点做敲低/抑制，测 module score 与 responder-like 位移是否按 PASCAR 预测移动——这是把反事实"做实"的最强证据。

> 具体平台（细胞系 / 类器官 / 小鼠模型 / 抗体试剂可得性）需你确认后定稿（第 11 节）。

---

## 8. 与对齐报告 Gap 的逐项闭合

| Gap | 闭合方式 |
|---|---|
| S1/T1 PD-1 anchor 断裂 | WP1 三 GEO 开放队列救援 + within-cohort meta；双轨退路 R2 |
| S2 仅关联非因果 | **WM-A 不变因果 + WM-B 反事实**——这是方法学创新的正面回应 |
| S3 HCC 集中度 | WP2 回灌 GSE229772 + GSE149614 深描 + WP4 HCC 外锚 + WP3 空间 |
| S4 髓系脆弱 | 髓系注释复审（marker panel/参考映射）+ WL IHC/mIF 验证 |
| S5 / T?? L1/L2/L3 未正式化 | **WM-C 证据图把 L1/L2/L3 做成算法输出** |
| T2 表达不均 | WP2 回灌（基础设施已建） |
| T3 强 ML/scVI 基线缺 | A.v0 GGM + scVI latent simple classifier 补基线（仅作 baseline，非创新） |
| T4 空间规则化 | WP3 定量化 |
| T5 外锚窄 | WP4 多队列激活（编号见 §5.2/5.3） |
| T6 扰动策展非因果 | WP5 perturb-prior 矩阵 + 作 A 层方向软约束 |
| T7 TCR 主表缺 | WP5 从已建克隆原始表构建 master（基础设施已建） |
| T8 投稿未就绪 | WP6 |
| M-signals 文件缺 | WM-B v0 直接产出正式文件 |

---

## 9. 被排除 / 可选未采用方案及具体理由（要求 3）

> 不做笼统归因，逐条给具体技术/科学理由。

**继续排除（维持原 plan 禁止项）**：
- **RL 顺序用药主线**：当前无纵向多线治疗的密集 response 反馈数据（timepoint unknown 占 7,150 行），RL 需要的状态-动作-奖励轨迹根本不存在，强上=拟合噪声；且与"机制提名非临床推荐"定位冲突。→ 排除。
- **Flow/ODE 连续动力学主建模**：需高质量 pre/on/post 配对密集时序，本项目配对样本稀疏（has_paired_pre_post 极少）；ODE 在稀疏时序上不可识别。→ 后置，作为 pre/post 充足后的 secondary dynamic 分析（与 plan §7 一致）。
- **网页 / 临床推荐系统**：定位为机制提名，监管与伦理不允许无前瞻验证的推荐。→ 排除。
- **LINCS reversal 单独提名药物 / marker 均值单独提机制**：单证据源易假阳，PASCAR 要求 ≥2 正交证据过 Layer C。→ 排除为单独依据，仅作 perturb 软约束。

**可选未采用为主线（给理由，保留为变体/对照）**：
- **基因级因果图（vs 模块级）**：基因级 DAG 在万级维度下不可识别、算力爆炸、生物可解释性差，且易过拟合 batch。→ 采用**模块级**因果（A 层主推），基因级仅在单个已确认模块内部做局部精修（可选 v2）。
- **黑盒图神经网络作主排序（C.v2）**：顶刊审稿人对"黑盒决定机制等级"高度警惕，可重复性与可解释性弱于贝叶斯/加权融合。→ 不作主结论，仅作敏感性对照。
- **基础模型 embedding（scGPT/Geneformer）作主表征**：可作跨队列 batch 协调的**可选表征底座**（算力足时值得一试，能提升 anchor 跨癌种对齐），但**不作主创新**——理由：当前可解释模块体系已是叙事核心，foundation embedding 不可解释、且会引入"它学到的是癌种还是免疫态"的新混杂，需额外 leave-cancer 审计。→ 列为 A 层可选上行，须配混杂审计。
- **retrieval / memory bank（v6 旧件）**：仅作结果解释层展示，不进核心建模——理由：检索相似患者≠机制因果，且与证据图功能重叠。→ 降级为展示。
- **可微因果发现 A.v2 作主线**：方法最炫但识别性依赖强假设、对模块定义敏感、调参成本高、时间紧。→ 作野心选项，主推 A.v1 不变因果（更稳、更可解释、更易过审）。
- **纯泛癌大一统模型**：会稀释 HCC 转化叙事且加重癌种混杂。→ 用 shared/specific 分解替代"大一统"。

---

## 10. 时间线（时间紧，先锁可发）

| 阶段 | 周期(估) | 并行任务 | 产出 |
|---|---|---|---|
| **S1 救援底座** | ~2–3 周 | WP1(anchor 整合) ∥ WP2(回灌) ∥ WM-C.v0(证据矩阵) | anchor 分支判决 + L1/L2/L3 矩阵 → **R 路线初选** |
| **S2 方法学内核** | ~3–4 周 | WM-A.v1 ∥ WM-B.v0→v1 ∥ WP3(空间定量) | SS-CGF 图 + 正式 mIMS/MRI/mICS + 空间裁判 |
| **S3 证据闭环** | ~3 周 | WP4(外锚) ∥ WP5(TCR/perturb) ∥ WL(IHC/mIF) | 多模态证据卡 + 湿验方向 |
| **S4 头条锁定+投稿** | ~2–3 周 | Milestone C′ 锁路线 ∥ WP6 ∥ C.v1 | 6 主图 + 手稿 + external_submission_ready |

> R3（证据图头条）在 S1 末即具备"可发雏形"；R1/R2 在 S2–S3 成形。即"先有保底、再冲主攻"。

---

## 11. 待你确认（要求 4：不确定即停）

以下 4 点会影响可执行性，我需要你拍板再细化对应分支（其余我已自行决策）：

1. **受控数据是否申请**：dbGaP（Bi RCC `phs002065`、Van Allen `phs000452`、Liu2019）与 EGA（Braun `EGAS00001004290-92`、Pozniak `EGAS00001006488`）需数据访问委员会审批、周期不定。**默认假设：不申请，仅用 GEO 开放队列**（已足够过 WP1 门）。如要申请请告知，我把它们排进上行计划。
2. **本地扰动数据是否就位**：报告把 `GSE133344/GSE90063/GSE193736/GSE306429/L1000` 记为"本地资源"。请确认是否已下载到 `data/`？若无，WP5 perturb 矩阵需先补下载（我可给清单）。
3. **湿实验平台可得性**：HCC 细胞系 / 类器官 / 小鼠模型 / IHC-mIF 抗体试剂中，哪些现成可用？这决定第 7 节湿验从"IHC 验证"到"CRISPR 扰动"能走到哪一档。
4. **HCC Visium 补强**：是否需要我对"ICI 相关 HCC 10x Visium"再做一轮定向精检索以补到精确编号（GSE238264 之外）？还是先用已暂存的 GSE238264 起步？

确认这 4 点后，我可把选定路线（R1/R2/R3）展开为带脚本清单与输出契约的 Step 级施工方案。

---

## 附录：精确数据编号速查

- **anchor scRNA(开放)**：GSE123813、GSE176021、GSE200996
- **HCC ICI bulk**：GSE215011、GSE235863、GSE140901、GSE279750（+IMbrave150 已有）
- **泛癌 ICI bulk**：GSE78220、GSE91061、PRJEB23709（+dbGaP: phs000452 / Liu2019）
- **RCC anchor(受控)**：dbGaP phs002065 / SCP1288；EGA EGAS00001004290-92
- **HCC 空间/深描**：GSE238264(已暂存)、GSE149614；EGA EGAS00001006488(Pozniak, 受控对照)
- **TCR**：GSE236581(已暂存)+ 三 anchor 自带 scTCR
- **扰动(本地待确认)**：GSE133344、GSE90063、GSE193736、GSE306429、LINCS-L1000
