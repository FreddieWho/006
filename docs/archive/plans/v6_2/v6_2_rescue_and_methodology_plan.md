# v6.2 修复与方法学创新方案（Rescue & Methodology Innovation Plan）

**生成时间**：2026-06-17
**取代**：`reports/v6_1_rescue_and_methodology_plan.md`（v6.1 版仍保留作沿革）
**依据**：`reports/v6_1_science_engineering_alignment_report.md`（二次修订）+ `docs/plan/v6_2/{plan,roadmap}.md` + `results/v6_1/` 实际进度 + 外部数据精检 + 学界近半年方向
**性质**：把 v6.1 救援方案按 v6.2 升级（因果+反事实+空间 niche 共主轴+证据三角化+可选虚拟细胞 ISP），保留全部精确数据编号与双轨门控。

---

## 0. 与 v6.1 救援方案的差异（一页速读）

| 维度 | v6.1 救援方案 | v6.2 升级 | 升级理由（近半年学界） |
|---|---|---|---|
| 方法学骨架 | PASCAR（因果/反事实/证据图） | PASCAR 不变，但**反事实直接 adapt CINEMA-OT/scDRP，因果 adapt GRIT 思路** | 这些是 2023–2026 已发表、可复现的 OT 反事实/因果 GRN 方法，降实现风险 |
| 空间 | 裁判层 | **升为共主轴，niche 级定量（ISCHIA/GraphTME 风格）** | 近半年 ICI 抗性前沿明确转向空间 niche 与细胞互作 |
| 虚拟细胞/scFM | 仅"可选上行" | **明确为可选 ISP 引擎 + 三角验证项，配 41% 批判性防御** | scFM/虚拟细胞爆发（STATE/Tahoe-x1/VCC），但有"未必超基线"的强批判 |
| 证据 | L1/L2/L3 矩阵 | **新增"三角化"为命名方法学贡献** | 顶刊对 robustness/可重复性审查趋严 |
| HCC 屏障 | 泛泛"髓系/基质/血管" | **具体 marker 假设（S100A9⁺单核/TREM2⁺巨噬/C1QC niche/SIRPα/CAF/VEGF）** | 近半年 HCC 抗性单细胞工作给出具体分子 |
| PD-1 双轨 | 路线建议 | **制度化为判决协议 + 阈值门（G5）** | v6.1 anchor 断裂的直接教训 |

---

## 1. 战略决策（沿用上一轮用户确认）

三方向整合为同一框架且三轨可并行；目标生物医学顶刊；PD-1 双轨按救援结果决定；算力足/湿实验可纳入/时间紧需快速可发；排除项须给具体理由。v6.2 在此基础上据近半年方向再优化。

---

## 2. 统一框架 PASCAR（v6.2 强化版）

> **PASCAR**：PD-1-**A**nchored **S**hared-specific **CA**usal module graph with counterfactual **R**epair scoring, **niche-adjudicated and triangulated**。

```
输入：module_score_matrix(已产) + PD-1 anchor 标签 + 多模态层(TCR/spatial-niche/perturb/ISP/bulk)
  │
  ▼ Layer A 结构层  共享=因果不变性 / HCC特异=环境残差边（模块级；扰动软约束）
  │                 adapt：不变因果发现(ICP/IRM) + GRIT 思路；产 causal DAG
  ▼ Layer B 干预层  module-ITE：do(repair) 沿 responder-like 流形的可达位移 = MRI
  │                 adapt：CINEMA-OT / scDRP 条件最优传输反事实；产正式 mIMS/MRI/mICS
  ▼ Layer C 裁判层  多模态正交证据贝叶斯融合 + 三角化 → L1/L2/L3
  │                 triangulation：{因果方向, 反事实ITE, 扰动/ISP, 空间niche, 外部bulk} ≥2 一致
  ▼ 输出：PD-1+X mechanism nomination + 最小湿验
  └─ 横切：空间 niche（共主轴，既是发现也是证据）；scFM ISP（可选第二引擎，三角验证）
```

**对生物医学顶刊的竞争力**：新颖性=问题建模+架构组合（因果不变性×反事实ITE×空间niche×多源三角化，以 PD-1 单药生物学为锚，HCC 转化闭环），而非单算法。所有组件踩在近半年最热的三条方法线（虚拟细胞ISP / 因果-OT反事实 / 空间niche）上，且每条都用已发表方法 adapt，审稿人既觉得"前沿"又觉得"可信可复现"。

各层 v0/v1/v2 分级与接口契约同 `docs/plan/v6_2/plan.md §6`，此处不赘。

---

## 3. 三条技术路线（v6.2 更新）

| 路线 | 配置 | 主叙事 | 适用 | 卖点 |
|---|---|---|---|---|
| **R1 全栈 PD-1 中心**（推荐主攻） | A.v1 + B.v1 + C.v1 + niche 共主轴 + ISP 三角；anchor 用 GSE123813/176021/200996 救援 | PD-1 单药因果敏感性锚 → 共享/HCC 特异因果屏障 → 反事实 X 修复 → niche 裁判 + 多源三角 → 湿验 | anchor 救援门(G5)过 | 完整野心，方法+生物双满，全踩近半年前沿 |
| **R2 HCC-first 因果转化**（稳健退路） | A.v1(HCC 特异层)+B.v0/v1(HCC responder-like)+C.v1，重仓 niche | HCC 免疫耐受下的特异因果屏障(髓系S100A9/TREM2、CAF、VEGF) → niche 解构 → PD-1+X(anti-VEGF/lenvatinib) 转化提名 | G5 不过 / HCC 信号强 | 转化强、对 anchor 不敏感、风险最低；直接对接近半年 HCC 抗性热点 |
| **R3 证据图+niche 方法学头条**（最快保底） | C.v1 头条 + niche 定量 + A/B 为组件 + ISP 三角 | 一个把因果结构、反事实可干预性、空间 niche、多源证据三角化为 L1/L2/L3 的机制提名框架，多癌种+HCC 演示 | 时间最紧先锁一篇 | 复用资产最多、最快出；方法学贡献清晰 |

**推荐**：R1 为目标、R3 为保底、R2 为 anchor 失败转向。三者共栈，Milestone C 处依 G5 + HCC/泛癌信号 + 时间一次性锁头条。

---

## 4. 外部数据救援精确清单（要求 2，沿用并补 niche/ISP）

> 走 addendum，不替换主输入。

### 4.1 干净 anti-PD1 单药 anchor scRNA（WP1，全 GEO 开放，足够过 G5）
`GSE123813`(BCC/SCC,+scTCR) · `GSE176021`(NSCLC 新辅助 nivo,+scTCR) · `GSE200996`(HNSCC 新辅助 nivo 剥联合臂,+scTCR)。三癌种交叉即满足"2+ 独立干净队列"，不依赖受控数据。受控上行：dbGaP `phs002065`(Bi RCC)、EGA `EGAS00001004290-92`(Braun)。

### 4.2 HCC ICI 外锚（WP4 / R2）
`GSE215011`(nivo 单药) · `GSE235863`(PD1+lenva,配对,直接服务 X-TKI 叙事) · `GSE140901`(NanoString) · `GSE279750`(PD-L1 联合) · IMbrave150(已有,升级模块投影+CI)。J Hepatol atezo+bev 422 样本若可得为最强（疑非全公开，需向作者/补充确认）。

### 4.3 泛癌 ICI 外锚（R1/R3 跨癌种）
`GSE78220`(Hugo) · `GSE91061`(Riaz) · `PRJEB23709`(Gide) · dbGaP `phs000452`(Van Allen,仅CTLA4参考)/Liu2019。

### 4.4 空间 niche（升为共主轴）
已暂存 `GSE238264`(HCC,support_only) · `GSE149614`(HCC scRNA 深描) · **需补 ≥1 套 ICI 相关 HCC 10x Visium**（含核/侵袭前沿/基质分区；Wang 2022 Theranostics HCC Visium 为候选，编号未公开需索取）· 可选受控 `EGAS00001006488`(Pozniak,scRNA+Visium+GRN，作 A 层方法对照)。

### 4.5 TCR / 扰动 / ISP
TCR：`GSE236581`(已暂存,1.1M 克隆)+ 三 anchor 自带 scTCR → `tcr_clone_master_v6_2.csv`。
扰动(本地待确认)：`GSE133344`/`GSE90063`/`GSE193736`/`GSE306429`/LINCS-L1000 → perturb-prior 矩阵。
ISP(可选)：公开 scFM（如 STATE/scOTM 类）对核心模块做 in-silico 扰动，与 OT 反事实+扰动先验三角，配 leave-cancer 审计。

---

## 5. 工作包重构（v6.2，含新增 niche/ISP 方法学 WP）

| WP | 内容 | 优先级 | 退出门 |
|---|---|---|---|
| **WP0** 合同加固 | internal-review caveat、PD1_anchor support-only、Phase3.5 invalid、scFM 仅可选条款 | 立即 | 无越界 claim |
| **WP1** anchor 救援 | GSE123813+176021+200996(剥联合臂)→Step1/2/3 anchor 分支+within-cohort meta | **P0** | G5：混杂<HIGH 且 within-cohort 方向一致→R1；否则 R2 |
| **WP2** 表达回灌 | 已建 GSE120575/229772 pseudobulk 接 Step2.8/2.9→再生 hotfix；移除幻影患者 gse120575_2 | P1 | 无零表达主 cohort |
| **WM-A** SS-CGF 因果 | A.v0→A.v1 不变因果（扰动软约束），模块级 shared/HCC-specific DAG | **P0(方法)** | shared backbone LOCO 稳定 + ≥1 HCC 特异因果屏障 |
| **WM-B** 反事实修复 | B.v0(正式 mIMS/MRI/mICS)→B.v1(CINEMA-OT/scDRP 风格 module-ITE) | **P0(方法)** | M-signals 文件存在 + 每 nomination 连回 anchor/HCC 流形 |
| **WM-N** 空间 niche（新增，升主轴） | niche 构建(ISCHIA 风格)+细胞互作图+模块映射+置换零控；GSE238264(+补强) | **P0(方法)** | ≥1 核心屏障过 niche 级定量支持 |
| **WM-C** 证据图三角化 | C.v0→C.v1 贝叶斯融合 + 三角化规则 → candidate_module_evidence_matrix.csv | **P0(方法)** | 每 candidate L1/L2/L3 + 三角化命中数 + 冲突降级 |
| **WP4** 外锚激活 | §4.2/4.3 多队列模块级投影+CI+失败案例 | P1 | 核心机制 ≥1 独立外锚复现方向 |
| **WP5** 扰动+ISP+TCR | perturb-prior 矩阵；可选 scFM ISP 三角；tcr_clone_master_v6_2.csv | P2 | 证据矩阵含 perturb/ISP/TCR 字段；ISP 不稳则退可选 |
| **WP6** 投稿复现 | Supp.S1+参数附录+环境锁(含 scFM/OT/因果方法版本)+claim audit | P0(投稿前) | external_submission_ready=true |
| **WL** 湿验 | IHC/mIF(髓系 niche,如 S100A9/TREM2/C1QC 共定位)→共培养/极化→CRISPR/小分子 | 与 WP4 并行 | ≥1 核心靶点体外方向一致 |

**双轨里程碑**：Milestone A′(WP1 后)→选 R1/R2 候选；Milestone C(WM-A/B/N/C v1 + WP4 部分后)→锁头条 R1/R2/R3。

---

## 6. 最小湿实验（湿验已纳入，按近半年 HCC 生物学定靶）

按可达性排序，优先验证近半年高频 HCC 髓系/基质/血管屏障：
1. **IHC/mIF 空间共定位**（先做、最低成本）：验证髓系-肿瘤巢 niche（S100A9⁺CD14⁺单核、TREM2⁺巨噬）、CAF 屏障、CD8 排斥带——直接把 WM-N 的 niche 证据从计算升为组织学。
2. **体外共培养/极化**：验证 X-class 方向（M2→M1 重编程、SIRPα 阻断释放髓系抑制、IFN/APC 恢复）对 T 细胞功能的因果方向，校验 WM-B 反事实方向。
3. **CRISPR/小分子扰动**（算力+湿验充足）：Top-1 靶点敲低/抑制，测 module score 与 responder-like 位移是否按 PASCAR 预测——把反事实"做实"，也为 scFM ISP 提供真值对照。

平台可得性需确认（见第 8 节）。

---

## 7. 被排除 / 可选未采用方案及具体理由（要求 3，含 v6.2 新评估）

**继续排除**：
- **RL 顺序用药主线**：无纵向密集 response 反馈轨迹（timepoint unknown 7,150 行），RL 状态-动作-奖励数据不存在，强上=拟合噪声；与机制提名定位冲突。
- **Flow/ODE 连续动力学主建模**：配对密集时序稀缺（has_paired_pre_post 极少），ODE 不可识别。后置为 pre/post 充足后的 secondary dynamic。
- **网页/临床推荐系统**：定位机制提名，监管伦理不允许无前瞻验证推荐。
- **单证据源提名**（LINCS reversal 单独提药 / marker 均值单独提机制）：易假阳，PASCAR 三角化要求 ≥2 正交证据。

**可选未采用为主线（给具体理由，保留为变体/对照）**：
- **scFM 作主表征/主 ISP**：近半年虽爆发，但有强批判证据——大规模 CRISPRi 仅约 41% 扰动有可测全转录组效应，scFM 在多项基准未稳超简单基线，且会引入"学到癌种还是免疫态"的新混杂。→ 仅作可选第二 ISP 引擎 + 跨癌种表征底座，配 leave-cancer 审计；三角化以 OT 反事实+扰动先验为准（关闭条件 5）。这是对前沿的"采纳但不押注"的防御性设计。
- **基因级因果图（vs 模块级）**：万级维度不可识别、算力爆炸、易过拟合 batch。→ 模块级为主，基因级仅在已确认模块内部局部精修(v2)。
- **黑盒 GNN 作主排序（C.v2）**：顶刊审稿人警惕黑盒定机制等级，可解释性/可复现性弱于贝叶斯融合。→ 仅敏感性对照。
- **可微因果发现 A.v2 作主线**：最炫但识别性依赖强假设、对模块定义敏感、调参贵、时间紧。→ 野心选项；主推 A.v1 不变因果（更稳、更可解释、更易过审）。
- **空间 GNN 主模型**：与 plan 禁止项一致，且 niche 构建(ISCHIA 风格)+互作图已足够支撑 Q4，黑盒空间 GNN 解释性差。→ 排除为主模型，niche 用可解释统计。
- **retrieval/memory bank（v6 旧件）**：检索相似患者≠因果，与证据图重叠。→ 降级展示层。
- **纯泛癌大一统模型**：稀释 HCC 叙事、加重癌种混杂。→ 用 shared/specific 因果分解替代。

---

## 8. 时间线（时间紧，先锁可发）

| 阶段 | 周期(估) | 并行 | 产出 |
|---|---|---|---|
| **S1 救援底座** | ~2–3 周 | WP1(anchor) ∥ WP2(回灌) ∥ WM-C.v0(证据矩阵) | G5 判决 + L1/L2/L3 矩阵 → R 路线初选（R3 此时已具可发雏形） |
| **S2 方法学内核** | ~3–4 周 | WM-A.v1 ∥ WM-B.v0→v1 ∥ WM-N(niche) | 因果图 + 正式 mIMS/MRI/mICS + niche 定量裁判 |
| **S3 证据闭环** | ~3 周 | WP4(外锚) ∥ WP5(扰动/ISP/TCR) ∥ WL(IHC/mIF) | 三角化证据卡 + 湿验方向 |
| **S4 头条锁定+投稿** | ~2–3 周 | Milestone C 锁路线 ∥ WP6 ∥ C.v1 | 6 主图 + 手稿 + external_submission_ready |

---

## 9. 与对齐报告 Gap 的逐项闭合（v6.2）

| Gap | v6.2 闭合 |
|---|---|
| S1/T1 anchor 断裂 | WP1 三 GEO 开放队列 + G5 双轨门 |
| S2 仅关联非因果 | WM-A 不变因果 + WM-B 反事实 ITE（直接 adapt 已发表方法） |
| S3 HCC 集中度 | WP2 回灌 + GSE149614 深描 + WP4 HCC 外锚 + WM-N niche |
| S4 髓系脆弱 | niche-ready 注释 + S100A9/TREM2/C1QC 假设 + WL IHC/mIF |
| S5/L1L2L3 | WM-C 证据图三角化做成算法输出 |
| T2 表达不均 | WP2（基础设施已建） |
| T3 强基线 | streamlined 六基线 + scVI latent |
| T4 空间规则化 | WM-N niche 级定量（共主轴） |
| T5 外锚窄 | WP4 多队列（编号见 §4） |
| T6 扰动策展 | WP5 perturb 矩阵 + 可选 scFM ISP 三角 |
| T7 TCR 主表 | WP5 从已建克隆原始表构建 master |
| T8 投稿 | WP6 |
| M-signals 文件缺 | WM-B v0 直接产文件 |

---

## 10. 待你确认（要求 4，沿用上一轮 4 点 + v6.2 新增 2 点）

1. **受控数据是否申请**：dbGaP(`phs002065`/`phs000452`/Liu2019) + EGA(`EGAS00001004290-92`/`EGAS00001006488`)。默认：不申请，仅 GEO 开放（已足过 G5）。
2. **本地扰动数据是否就位**：`GSE133344/90063/193736/306429/L1000` 是否已下载到 `data/`？无则先补。
3. **湿验平台可得性**：HCC 细胞系/类器官/小鼠/IHC-mIF 抗体哪些现成？决定湿验走到 IHC 还是 CRISPR。
4. **HCC Visium 补强**：是否要我定向精检索补精确编号，还是先用 GSE238264 起步？
5. **（新增）scFM ISP 是否启用**：算力足但有 41% 批判证据。默认：作可选三角项，不进关键路径；若你要把"虚拟细胞 in-silico 扰动"做成卖点之一，我把 WP5 的 ISP 升为 P1 并加基准对照（vs OT 反事实 vs 扰动先验 vs 简单基线）。
6. **（新增）头条路线偏好**：默认按 G5+信号在 Milestone C 自动裁决（R1>R3>R2 优先级）。若你已倾向某条（如直接押 R2 HCC-first 以求稳发），告知我提前收敛施工。

确认后，我把选定路线展开为带脚本清单与输出契约的 Step 级施工方案（对接 `docs/plan/v6_2/roadmap.md`）。

---

## 附录：精确数据编号速查（v6.2）
- anchor scRNA(开放)：GSE123813、GSE176021、GSE200996
- HCC 深描：GSE149614；HCC 空间：GSE238264(已暂存)
- HCC ICI bulk：GSE215011、GSE235863、GSE140901、GSE279750（+IMbrave150）
- 泛癌 ICI bulk：GSE78220、GSE91061、PRJEB23709（+dbGaP phs000452/Liu2019）
- 受控上行：dbGaP phs002065、EGA EGAS00001004290-92、EGAS00001006488
- TCR：GSE236581(已暂存)+ 三 anchor 自带 scTCR
- 扰动(本地待确认)：GSE133344、GSE90063、GSE193736、GSE306429、LINCS-L1000

## 附录：方法学 adapt 来源（近半年，可直接复用降风险）
- 反事实 ITE：CINEMA-OT(Nat Methods 2023)、scDRP(2026 条件 OT 反事实/ITE)
- 因果 GRN：GRIT(Bioinformatics 2025, OT+ODE 因果)、不变因果(ICP/IRM)
- 空间 niche：ISCHIA(niche 构建)、GraphTME(2025 空间互作图预测 ICI)
- 虚拟细胞 ISP：STATE(Arc)、scOTM、Tahoe-x1、Virtual Cell Challenge(Cell 2025)（可选）
- HCC 生物学靶：S100A9⁺CD14⁺单核、TREM2⁺巨噬、巨噬源血管生成、SIRPα/PD-L1⁺髓系、C1QC niche、CAF
