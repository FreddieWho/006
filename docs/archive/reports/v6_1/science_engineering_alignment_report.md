# v6.1 科学目标与工程目标对齐性报告（二次修订版）

**生成时间**：2026-06-17（二次修订：在 2026-06-17 初版 + 2026-06-15 Gap Report 基础上，补充 2026-06-11 新数据增补包的产物核验）
**分析范围**：`docs/plan/v6_1/plan.md`、`docs/plan/v6_1/roadmap.md`、`results/v6_1/INDEX.md`、`results/v6_1/V6_1_PLAN_ROADMAP_GAP_AND_RESCUE_REPORT_20260615.md`、`results/v6_1/V6_1_PHASE_SYSTEMATIC_ASSESSMENT_20260611.md`、`results/v6_1/new_data_integration_addendum_20260611/ADDENDUM_COMPLETION_REPORT_20260611.md` 及其 `05_phase4_merge_decision/*.yaml`
**方法**：CodeGraph 索引 + 关键文档直读 + Gap Report 系统比对 + 增补包产物逐项核验（merge_decision YAML 与 completion report 双源交叉确认）

> **二次修订要点**：项目主流水线自 2026-06-15 起未推进，故 06-17 初版的主链结论与对齐评分整体维持。本次修订的唯一新增信息来自对 2026-06-11 增补包（`new_data_integration_addendum_20260611/`）产物的逐项核验——初版与 Gap Report 将增补工作笼统记为"进行中"且把 TCR 层记为"整体缺失"，**实际上 P0 批次已完成裁定并落地了 TCR 克隆原始表、两个表达救援队列的 pseudobulk、以及第二个空间队列 manifest**。但这些产物按治理规则（INDEX §9、§6.4）**全部为 `support_only` / `primary_addendum_candidate`，尚未并入 Step2 hotfix / Step3 / Phase4 主流水线**，因此不改变主链对齐评分，只把若干 Gap 从"缺失"精确化为"基础设施已建、主流水线未并入"。详见新增第 5.4 / 6.4 节与第 13 节。

---

## 1. 执行摘要

| 维度 | 结论 |
|---|---|
| **总体对齐度** | 中等（约 60–70%）。项目治理结构优秀、证据包内部就绪，但科学主链在 PD-1 anchor 处断裂，多个关键方法学目标未实现。 |
| **科学目标对齐** | 方向一致，核心链断裂。"PD-1 anchor → residual barrier → X repair"因果链在第一环断裂（PD1_anchor 仅有 25 名有响应标签患者来自 1 个队列，混杂不可消除）。已降级为"多队列关联证据包"而非"PD-1 中心提名框架"。 |
| **工程目标对齐** | 流程骨架强对齐，交付深度不足。Step1–Phase10 门控流水线已完整执行；但 PD1 anchor 模块输出缺失、强 ML 回滚、空间裁判规则化而非定量化、外部锚点覆盖不足、TCR 层缺失、mIMS/MRI/mICS 正式输出文件不存在。 |
| **关键不对齐** | ① PD1 anchor 主监督被禁止；② Causal GRN 未实现；③ Phase3.5 强 ML 回滚；④ 空间裁判为规则赋标而非定量模型；⑤ 外部锚点仅 IMbrave150；⑥ TCR 层整体缺失；⑦ mIMS/MRI/mICS 文件未生成。 |
| **最终状态** | Phase10 **CONDITIONAL_PASS**：内部评审就绪，外部投稿未就绪，3 个 Methods 阻塞项 + 更深层科学门槛问题。 |
| **最适当前定位** | "经过门控审计的多队列单细胞 ICI 关联证据框架，含 HCC 优先转化提名候选"，而非"PD-1 中心 TME 机制提名框架"。 |

---

## 2. v6.1 规划文档中的科学目标

### 2.1 核心科学问题

`docs/plan/v6_1/plan.md` 明确 v6.1 要回答三个问题：

1. **Q1**：跨癌种可复现的 PD-1 敏感性基础模块是什么？
2. **Q2**：PD-1 不够时，主要残余障碍是什么？
3. **Q3**：哪些 X 机制类最可能修复这些障碍？

总目标是：**在多癌种 ICI 数据中学习可复现的共享/特异 TME 免疫模块，以 PD-1 单药/ICI 单药反应作为基础免疫敏感性锚点，进一步解释 PD-1+X 中 X 所修复的残余免疫障碍，并输出可被空间组织证据、扰动先验和最小实验支持的机制类 X 与候选靶点。**

输出名称为 **PD-1+X hypothesis nomination framework**（机制假设提名系统），明确排除"临床推荐系统"。

### 2.2 关键科学原则

| 原则 | 规划要求 |
|---|---|
| **PD-1-centered continuum** | PD-1/ICI monotherapy 是 anchor（学习基础敏感性），PD-1+X 是 extension（解释 X 修复的残余障碍），高混杂联合治疗只做敏感性分析。三层构成连续整体，不可割裂。 |
| **强基线前置（硬门槛）** | fraction、signature、pathway、elastic net/XGBoost、static module、scVI latent 六类基线必须先于任何复杂模型结果。复杂模型须在预测性能/泛化稳定/机制解释/外部一致性/湿实验候选清晰度上至少一项超过强基线。 |
| **三档证据等级** | L1 Robust association → L2 Mechanistic support → L3 Translational nomination。L3 必须建立在 L1+L2 上；任何关键冲突自动降级。 |
| **空间裁判核心地位** | 空间/组织数据是 biological adjudicator，无空间支持的 TME 机制 claim 不得升为 L2/L3 主结论。 |
| **机制提名而非临床推荐** | 输出叫 PD-1+X hypothesis nomination，不叫临床推荐系统。 |
| **MVP 审慎继承** | 只继承数据底座；手工 MOA、heuristic counterfactual、demo ranking 不得作为主结果。 |

### 2.3 规划中的主结果

- 至少 1 个 shared PD-1 sensitivity module（L1 认证）
- 至少 1 个 HCC-specific resistance/barrier module（L1 认证）
- 至少 1 个经空间/组织裁判支持的机制链（L2 认证）
- 至少 1 个可实验落地的候选靶点（L3 认证）
- 6 张主图 + 完整复现包
- 正式 mIMS/MRI/mICS 表文件

---

## 3. v6.1 规划文档中的工程目标

### 3.1 分步执行路线（Roadmap）

`docs/plan/v6_1/roadmap.md` 将项目拆分为 11 步 + 5 个合练里程碑：

| 步骤 | 工程目标 | 关键交付物 |
|---|---|---|
| Step 0 | 范围锁定 | `scope_lock_v6_1.md`、`analysis_contract.md`、`terminology_dictionary.md` |
| Step 1 | 数据资产盘点 | `cohort_registry_v6_1.csv`、`sample/patient_metadata_master.csv`、`treatment_context_flags.csv` |
| Step 2 | 免疫状态测量 | cell fraction、pseudobulk、signature/pathway/TF scores、myeloid QC |
| Step 3 | 强基线包 | 六类基线（fraction/signature/pathway/elastic net·XGBoost/static module/scVI latent） |
| Step 4 | 共享/特异模块发现 | `shared_immune_modules.csv`、`hcc_specific_modules.csv`、stability report |
| Step 5 | PD-1 monotherapy anchor | `PD1_anchor_modules.csv`、sensitivity/resistance axes、patient scores |
| Step 6 | PD-1+X repair logic | `PD1X_repair_logic_table.csv`、residual barriers、mechanism class mapping |
| Step 7 | 空间/组织裁判 | spatial module scores、neighborhood scores、LRG spatial support（定量模型） |
| Step 8 | 扰动支持 + mIMS/MRI/mICS | `mIMS_scores.csv`、`MRI_scores.csv`、`mICS_scores.csv`、candidate targets |
| Step 9 | Bulk/external anchor | bulk projection scores、external validation results（非 trend demo） |
| Step 10 | Evidence cards + nomination | `PD1X_nomination_cards.csv`、`wetlab_priority_list.csv`、L1/L2/L3 显式标注 |
| Step 11 | 论文主图 + 复现包 | 6 张主图、reproducibility manifest |

### 3.2 硬性执行门

- **G0**：未完成 registry/context lock，不得进入建模。
- **G1**：未完成 strong baseline report，不得提交复杂模型主结果。
- **G2**：未形成 PD-1 anchor，不得进入 PD-1+X repair nomination。
- **G3**：无 spatial adjudication，不得发布 L2/L3 机制主结论。
- **G4**：所有主结论必须可追溯到 output_manifest 与中间表。

---

## 4. 当前项目实际进度

### 4.1 总体状态

| 指标 | 状态 |
|---|---|
| Verdict | **Conditional Go**（非完整 Go） |
| 内部评审包 | 就绪 |
| 外部投稿就绪 | **否**（3 个 Methods 阻塞项 + 深层科学门槛） |
| Phase3.5 强 ML | **已回滚，不可作为证据（有效行=0）** |
| PD1 anchor 主监督建模 | **被禁止**（cohort-response 混杂 p=3.6×10⁻⁵，高风险） |
| PD1 anchor 实际规模 | 25 名已知标签预处理患者，仅来自 1 个可执行队列（TASK01） |
| 可用主宇宙 | `HCC_specific`、`PD1X_extension`、`pan_cancer_shared` |
| mIMS/MRI/mICS 文件 | **未生成**（results/v6_1 下不存在） |

### 4.2 Step1–Step3 实际状态

| 项目 | 数值/状态 |
|---|---|
| cohort 总数 | 68（33 main、11 sensitivity、6 external、1 spatial_only、5 perturb_prior_only、12 excluded） |
| 样本行数 | 9,839 |
| 患者数 | 1,538（main branch: 642 患者） |
| 主分析样本 | 约 1,798 行（33 个 main cohorts） |
| 宽响应标签已知 | 1,307 行 / 12 cohorts |
| 严格响应标签已知 | 706 行 / 7 cohorts |
| response_binary=unknown | 1,461 |
| response_binary=not_applicable | 7,352 |
| timepoint=unknown | 7,150（主因 GSE140228_smartseq2，sensitivity-only） |
| treatment_context | PD1X_extension:27, external_anchor_only:21, high_confounding_support:14, PD1_ICI_anchor:6 |
| Step2 细胞总数 | 14,565,583（main 纳入 8,514,050） |
| Step2 特征矩阵 | 2,292 sample × 1,485 feature；1,129 patient × 8,925 feature |
| 零表达 cohort | GSE120575、GSE229772（无法提供 pseudobulk/signature/pathway/TF 特征） |
| signature/pathway 缺失 cohort | 6 个（GSE120575、GSE225063、GSE229772、GSE272993、GSE273718、krishna_2021_rcc） |
| Step3 有效基线行 | 108 行（Phase3.4 可解释基线，A:48、B:26、C:24、D:10） |
| Step3 强 ML 有效下游行 | **0**（Phase3.5 回滚） |

### 4.3 Phase4–Phase7 实际状态

| Phase | 关键产出 | 数量 |
|---|---|---|
| Phase4 | Tier A/B 主模块（shared:60, HCC:54, PD1X:61）；Tier C/support 模块 | 175 主 + 32 支持 = 207 总 |
| Phase5 | 机制组（Priority 1:15、Priority 2:24、Priority 3:12） | 51 |
| Phase6 | 输入机制组 → 靶点预筛行（P1:120、P2:15、探索:60）；**无最终靶点推荐** | 39 输入 / 195 预筛 |
| Phase7 | 主候选（core:4、support:21、support-only:12、pending:2）；外部支持 18/冲突 9 | 39 行总计 |
| Phase7 空间支持 | spatially_supported:10、tissue_supported:2、weakly_supported:9、pending:6 | 27 主候选中 |

### 4.4 Phase8–Phase10 实际状态

| 项目 | 状态 |
|---|---|
| Phase8 | 6 张主图冻结 + 3 个最小验证计划锁定 |
| Phase9 | 手稿初稿、图例、补充包、审稿风险响应表 |
| Phase10 | 投稿候选包（SVG+spec），内部评审就绪；外部投稿未就绪 |
| adversarial review fatal blockers | 0 |
| claim-audit forbidden claims | 0 |
| Methods 外部投稿阻塞项 | 3（Supp. Table S1 缺失、Phase4/Step3 参数补充缺失、软件环境锁缺失） |

### 4.5 新增并行工作（2026-06-11 增补包，已核验）

| 工作流 | 内容 | 状态 |
|---|---|---|
| `data_completeness_strengthening_20260611/` | P0/P1/P2 数据候选注册、本地未整合资产清单、Step1 增补 source records | 评估完成 |
| `new_data_integration_addendum_20260611/` | GSE236581、GSE238264、GSE120575_rescue、GSE229772_rescue 的元数据冻结 + h5ad 转换 + TCR 克隆表 + pseudobulk + QC + Step3 addendum audit + merge decision | **P0 批次完成裁定；P1 批次未启动；产物未并入主流水线** |

**P0 批次四个数据集的实际裁定与落地产物**（经 `merge_decision.*.yaml` 与 completion report 双源核验）：

| 数据集 | 裁定 | 样本/患者/细胞 | 响应标签 | 已落地关键产物 |
|---|---|---|---|---|
| GSE236581（CRC anti-PD-1 + TCR） | `support_only` | 169 / 22 / 975K | unknown（待 PMID 38981439 补充材料） | `h5ad_addendum` + `tcr_clone_addendum.GSE236581.tsv`（**1.1M TCR 记录 / 605K 唯一克隆 / 182 样本**）+ pseudobulk（36K 基因×169）+ QC |
| GSE120575_rescue（黑色素瘤 PD1X） | **`primary_addendum_candidate`** | 43 / 32 / 16K | Responder/Non-responder（已恢复） | TPM pseudobulk（55K 基因×43）+ Phase4 audit；待移除主线幻影患者 `gse120575_2` |
| GSE229772_rescue（HCC 相关） | `support_only` | 31 / 11 / 58K | unknown（CASCADE 象限/生存终点，非 RECIST） | norm-count pseudobulk（15K 基因×31）+ QC |
| GSE238264（HCC 空间） | `support_only` | 7 / 7 / 空间 | inferred R/NR | `spatial_sample_manifest.GSE238264.tsv`（7 样本） |

**关键治理事实**：以上全部产物按 INDEX §6.3「不可将未整合的新数据集作为验证证据」与 §9「新数据整合请走 addendum 流程，不可直接替换 Phase4 主输入」**未并入** Step2 hotfix 主矩阵、Step3 宇宙、Phase4 模块发现或 Phase7 验证。因此它们**降低了多个 Gap 的救援成本（基础设施已就位），但当前不改变任何主链 claim 与对齐评分**。GSE236581 三个被阻断门为：response labels（BLOCKED）、confounding audit（BLOCKED，因无响应）；GSE229772 同。GSE236581 为 CRC，主要服务 pan-cancer shared 与 TCR 支撑，非 HCC 深描。

---

## 5. 科学目标 vs 现状对齐性分析

### 5.1 对齐良好的方面

| 规划科学目标 | 当前实现 | 对齐度 |
|---|---|---|
| PD-1-centered continuum 治疗结构 | 四类治疗上下文分层到位；Anchor/Extension/Confounding/External_anchor_only 语义清晰 | ✅ 高 |
| 强基线前置精神 | Phase3.4 可解释基线保留；复杂模型未强行主叙事 | ✅ 高 |
| 机制提名而非临床推荐 | claim audit 0 违规；Phase10 明确禁止 clinical/drug recommendation 用语 | ✅ 高 |
| MVP 审慎继承 | 手工 MOA、demo ranking 未作为主结果；archetype/axis 降级为 proxy | ✅ 高 |
| HCC-first + pan-cancer shared 双轨 | HCC_specific / pan_cancer_shared / PD1X_extension 三个宇宙均已进入主分析 | ✅ 中高 |
| 模块发现（inferred module 层） | 175 Tier A/B 模块、51 机制组已生成，跨队列稳定性有报告 | ✅ 中高 |

### 5.2 科学 Gap（S1–S5）

#### Gap S1：PD-1-centered continuum 科学链第一环断裂

规划核心逻辑：PD-1 anchor 学习基础敏感性 → PD-1+X 解释残余障碍。当前 PD1_anchor 仅有 25 名有标签预处理患者（1 队列 TASK01），主导队列样本比 0.788，cohort-response chi-square p=3.6×10⁻⁵，混杂等级 HIGH。

**当前可防守措辞**：
> We identify audited HCC/ICI and pan-cancer immune mechanism candidates, with PD1_anchor used only as support/sensitivity biological context.

**不可防守**：
> We learned a robust PD-1 monotherapy sensitivity axis and nominate X mechanisms that repair residual PD-1 failure.

#### Gap S2：机制组是关联候选，非因果机制

Phase4/5/6 均明确禁止因果证明。模块稳定可解释，但未建立因果方向（无 causal GRN、无正式伪 bulk 混合模型、无 within-cohort meta-analysis 公式化输出）。

#### Gap S3：HCC 特异性受 cohort 集中度制约

HCC_specific 为 primary_allowed，但部分预览中 HCC_specific 特征宇宙低至 112 样本/44 患者/1 队列。GSE229772 在主流水线仍为表达特征缺失，直接削弱 HCC 覆盖。

> **二次修订**：GSE229772_rescue 已在增补包内产出 norm-count pseudobulk（15K 基因×31 样本，support_only），但响应标签仍 BLOCKED（CASCADE 象限/生存终点，非 RECIST），且未并入主流水线——故主流水线 HCC 集中度风险**未实质缓解**。增补包另落地 GSE238264 HCC 空间 manifest（7 样本，support_only），为 HCC 空间裁判提供了第二个候选数据源但尚未量化。

#### Gap S4：髓系 claim 高价值但脆弱

核心机制 CAND_P5MG_HCC_012、CAND_P5MG_PD1X_014 涉及髓系重编程，但 Step2 报告显示髓系注释不确定性与冲突较高，未经 IHC/mIF 独立验证。

#### Gap S5：L1/L2/L3 证据等级未正式对标

Phase10 使用 PASS/CONDITIONAL_PASS、Priority A/B/C、core/support/pending 等分级，与规划的 L1/L2/L3 证据层级不一一对应。缺少显式的 `candidate_module_evidence_matrix.csv`（含 L1/L2/L3 字段与降级规则）。

### 5.3 科学目标对齐性矩阵

| 科学目标 | 规划强度 | 实际强度 | 对齐评级 | 说明 |
|---|---|---|---|---|
| PD-1-centered continuum | 核心 | 治理层满足 | 🟡 中对齐 | 分层到位，但科学链在 anchor 处断裂 |
| PD-1 anchor 主监督 | 核心 | 被禁止 | 🔴 弱对齐 | 最大计划外降级 |
| Shared+specific causal GRN | 核心（v5.3遗留） | 未实现 | 🔴 弱对齐 | Phase4 为关联/聚类模块，非因果图学习 |
| 强基线前置 | 硬门槛 | 部分满足 | 🟡 中对齐 | Phase3.4 有效；Phase3.5 回滚；scVI latent 基线缺失 |
| 空间裁判 | 关键裁判 | 规则化支持 | 🟡 中对齐 | 非定量 spot-level 模型，仅分类赋标 |
| mIMS/MRI/mICS | 主读数 | 文件未生成 | 🔴 弱对齐 | 无正式输出文件，Phase6 规则分数为雏形 |
| TCR 层 | 支撑证据 | 原始表已建/未并入 | 🔴 弱对齐 | 主表 `tcr_clone_master_v6_1.csv` 不存在；但增补层有 GSE236581 TCR 克隆原始表（1.1M 记录，support_only） |
| 外部锚点 | 关键验证 | 仅 IMbrave150 | 🟡 中对齐 | TCGA/GTEx/ICGC/GSE179994/GSE200996 未激活 |
| 机制提名（非临床推荐） | 核心 | 强满足 | 🟢 强对齐 | claim boundary 控制良好 |
| L1/L2/L3 显式对标 | 重要 | 未正式化 | 🟡 中对齐 | 内部等级系统存在但未与规划三档对应 |

### 5.4 增补层产物 vs 主流水线并入状态对账（二次修订新增）

本节是二次修订的核心增量。它把"救援所需基础设施"与"主链是否真的用上了"分开记账，避免把"文件存在"误读为"主张成立"。

| 救援目标 (WP) | 规划要求 | 增补层是否已建产物 | 是否并入主流水线 | 对主链 claim 的影响 |
|---|---|---|---|---|
| WP1 PD1 anchor 升级 | 2+ 独立干净 PD-1 预处理队列 | 否（仅 GSE120575_rescue 含 R/NR，但为黑色素瘤 PD1X，非干净 monotherapy anchor） | 否 | **无**：anchor 仍 support-only，主链第一环仍断 |
| WP2 表达特征救援 | GSE120575/GSE229772 等非零表达 | **部分是**：两队列 pseudobulk 已产出 | 否（未回灌 Step2/Step3/hotfix） | 无：主流水线仍零表达 |
| WP3 定量空间裁判 | spot-level 模块投影 + 置换控制 | 否（仅 GSE238264 manifest，7 样本，support_only） | 否 | 无：空间仍为规则赋标 |
| WP4 外部锚点激活 | TCGA/GTEx/ICGC/GSE179994/GSE200996 | 否 | 否 | 无：外部仍仅 IMbrave150 |
| WP5 扰动+TCR | perturb-prior 矩阵 + TCR master + mIMS/MRI/mICS | **部分是**：GSE236581 TCR 原始克隆表已建（1.1M 记录） | 否（无 TCR master、无 perturb 矩阵、无 M-signals 文件） | 无：TCR/扰动/M-signals 仍不在主证据 |
| WP6 投稿复现性 | Supp.S1 + 参数附录 + 环境锁 | 否 | 否 | 无：3 个 Methods 阻塞项仍在 |

**结论**：增补层使 **WP2、WP5（TCR 分支）、并对 WP3 提供了第二空间数据源**——救援路径已被部分"预铺设"，工程上显著降险；但**没有任何增补产物跨过治理门并入主链**，故 06-17 初版的科学/工程对齐评分（第 12.2 节）整体维持不变。换言之：地基旁边已经堆好了部分钢筋，但还没有一根浇进主体结构。

---

## 6. 工程目标 vs 现状对齐性分析

### 6.1 对齐良好的方面

| 规划工程目标 | 当前实现 | 对齐度 |
|---|---|---|
| 分阶段门控流水线 | Step1–Phase10 完整执行，每个阶段有 FINAL_DECISION.yaml | ✅ 高 |
| 数据治理与 registry | 68 cohorts、sample/patient metadata masters、treatment_context_flags、frozen splits 均已生成 | ✅ 高 |
| 特征工程标准化 | Step2 产生 14.5M 细胞、统一特征矩阵、Repair B2 hotfix 已冻结 | ✅ 高 |
| 输出契约与 manifest 文化 | output_manifest_template.yaml、RUN_REGISTRY.csv、各 phase package index 已建立 | ✅ 高 |
| 审计与 scope lock 习惯 | 每个 phase 有报告、决策文件、adversarial review、claim boundary audit | ✅ 高 |
| 测试与可复现性（内部） | scripts/v6_1/、reproducibility manifest 已建立 | ✅ 中高 |

### 6.2 工程 Gap（T1–T8）

#### Gap T1：PD1_anchor 是 Roadmap 骨干但已被封锁

Roadmap Step5 未完成为主分析。Step6 PD-1+X repair logic 无法锚定到可防御的 PD-1 轴。下游所有 PD1X claim 必须保持"repair-like mechanism class"而非真正的"PD-1 residual barrier logic"。

#### Gap T2：特征完整性不均匀（主流水线仍零表达；增补层已部分救援）

主流水线中 GSE120575、GSE229772 在 Step2-to-Step3 阶段零表达特征；六个 cohort（+GSE225063、GSE272993、GSE273718、krishna_2021_rcc）缺少 signature/pathway 可用性。GSE120575 有宝贵响应标签但主流水线无法提供表达衍生特征。

> **二次修订**：增补包已对两个最高价值队列产出 pseudobulk——GSE120575_rescue（TPM，55K×43，**已升级为 `primary_addendum_candidate`**，含 Responder/NR 标签）与 GSE229772_rescue（norm-count，15K×31）。即 WP2 的表达救援基础设施在增补层**已部分完成**，但尚未回灌 Step2.8/2.9、未再生 hotfix 矩阵、未进入 Step3 宇宙。**主流水线零表达事实不变**；待办包括移除 GSE120575 主线幻影患者 `gse120575_2`。

#### Gap T3：强 ML/scVI 基线不完整

Phase3.5 强 ML 全面回滚（有效下游行=0）。规划要求的 elastic net/XGBoost/scVI latent simple classifier 基线未全部到位。无有效的复杂模型 vs 基线比较。

#### Gap T4：空间裁判为规则赋标，非定量模型

`run_phase7_validation_package.py` 通过 `spatial_label()` 按机制类规则分配 spatially_supported / tissue_supported / weakly_supported 标签，非 spot-level 邻近统计/置换控制的定量分析。Fig4 类空间 claim 必须维持"支持性空间/组织证据"而非"定量空间验证"。

#### Gap T5：外部验证范围窄

仅 IMbrave150 作为门控外部锚点（方向一致性，非全模块投影 + 置信区间 + 失败案例报告）。TCGA/GTEx/ICGC/GSE179994/GSE200996 仍是声明的缺口。Phase7 中 18 外部支持 / 9 冲突候选，冲突案例必须保持 support/pending 而非被忽略。

#### Gap T6：扰动证据为策展前筛，非因果支持

Phase6 使用规则权重计算靶点分数（priority + curated perturb support + druggability + direction clarity + rank penalty）。Xatlas 缺失。无正式 perturb-prior 矩阵。扰动支持计数：25 moderate curated prior、12 pending direction、2 limited curated prior。

#### Gap T7：TCR 主表/特征整合缺失（原始基础设施已建于增补层）

主流水线的 `tcr_clone_master_v6_1.csv` **仍不存在**，TCR 衍生特征（克隆扩增/多样性/克隆型共享/耗竭·细胞毒连接）未进入任何 Phase4–7 证据，故 T 细胞 dysfunction/细胞毒/修复 claim 仍缺少克隆支持。

> **二次修订（重要纠正）**：初版/Gap Report 称"TCR 层整体缺失"过于绝对。增补包已落地原始 TCR 基础设施 `tcr_clone_addendum.GSE236581.tsv`（**1.1M 条 TCR 记录、605K 唯一克隆、182 样本**）+ TCR join audit，merge 裁定 `support_only`，允许用途含 `TCR_clone_support`/`PD1_directional_validation`。准确表述应为：**TCR 原始克隆表已建（support_only，未并入），但 v6.1 TCR 特征主表与 TCR-机制联动证据尚未构建**。WP5 的 TCR 前置因此已显著降险——剩余工作是 leakage/split 审计、构建 master、并入证据矩阵；注意 GSE236581 为 CRC，服务 pan-cancer/TCR 支撑而非 HCC。

#### Gap T8：手稿包内部就绪，投稿未就绪

3 个已知 Methods 阻塞项：Supp. Table S1（队列/样本/患者精确计数）、Phase4/Step3 参数补充、软件环境锁/包版本。

### 6.3 工程目标对齐性矩阵

| 工程目标 | 规划要求 | 实际状态 | 对齐评级 |
|---|---|---|---|
| Step0 范围锁定 | 5 个治理文档 | scope_lock、analysis_contract、terminology_dictionary 已生成 | 🟢 强 |
| Step1 数据资产盘点 | 3 主表 + 4 标记表 | 68 cohorts、metadata masters、flags、leakage log 已生成 | 🟢 强 |
| Step2 免疫状态测量 | 5 类特征 + myeloid QC | 14.5M 细胞、统一矩阵、Repair B2 hotfix；2 个 cohort 表达缺失 | 🟡 中高 |
| Step3 强基线包 | 6 类基线 | Phase3.4 可解释有效；Phase3.5 强 ML 回滚；scVI latent 基线缺失 | 🟡 中 |
| Step4 模块发现 | shared/HCC-specific 模块 | 175+32 模块已生成，稳定性报告存在 | 🟢 强 |
| Step5 PD-1 anchor | anchor 模块/轴/分数 | **被禁止**；`PD1_anchor_modules.csv` 未生成 | 🔴 弱 |
| Step6 PD-1+X repair | repair logic + mechanism mapping | Phase5/6 机制组 + 候选已生成，但非严格 PD1 anchor → X 链 | 🟡 中 |
| Step7 空间裁判 | 5 类空间定量分数 + report | GSE238264 规则赋标，非定量 spot-level 模型 | 🟡 中 |
| Step8 扰动 + mICS/MRI/mIMS | 正式 M-signals 文件 | **`mIMS_scores.csv`、`MRI_scores.csv`、`mICS_scores.csv` 不存在** | 🔴 弱 |
| Step9 外部锚点 | bulk projection + external validation | 仅 IMbrave150；TCGA 等未激活 | 🟡 中 |
| Step10 nomination + L1/L2/L3 | evidence cards + wet-lab plan + 显式等级 | 4 核心/21 支持/2 pending + 3 验证计划；L1/L2/L3 未正式对标 | 🟡 中高 |
| Step11 手稿/复现 | 6 主图 + 投稿就绪 manifest | SVG + manuscript release candidate 内部就绪；外部投稿未就绪 | 🟡 中高 |
| 合同路径 S0/S1/S2 对齐 | no-bypass 规则可审计 | 路径漂移已记录，symlink 未执行 | 🟡 中 |

---

## 7. 当前安全声明 vs 不安全声明

基于 `V6_1_PLAN_ROADMAP_GAP_AND_RESCUE_REPORT_20260615.md` 明确划定：

### 7.1 当前安全声明（Safe now）

- v6.1 建立了经过门控多队列 scRNA/ICI 特征与 metadata 骨干。
- Step2 hotfix 修复了 TF activity，产出了条件注意事项下可用的下游特征资产。
- Step3 识别了主要允许的非 PD1 宇宙，封锁了不安全的 PD1_anchor 主监督建模。
- Phase4 从主证据池发现了经审计、与响应相关的免疫模块。
- Phase5 将模块压缩为机制组。
- Phase6 产出了靶点预筛和验证设计，不是最终靶点推荐。
- Phase7 产出了小规模核心集的门控关联性外部/空间/组织支持。
- Phase8–10 产出了内部评审手稿和图包。

### 7.2 当前不安全声明（Unsafe now）

- 稳健的 PD-1 单药监督预测器。
- 基于 PD-1 anchor 的 PD-1+X 修复提名作为主要 claim。
- 因果 TME 机制证明。
- 最终药物/靶点推荐。
- 临床治疗指导。
- Spot-level 空间临床模型。
- HCC 机制泛化到所有癌种。
- 复杂 ML 验证。

---

## 8. MVP/v5.3 承接分析

### 8.1 规划中的承接

`mvp/MVP_TO_V53_TRANSITION_ROADMAP_2026-04-17.md` 建议 9 个 Stage，核心是：保留数据底座，升级因果 GRN/M-signals/ICSmini/e_X-mini/验证框架，后置 Flow/ODE。

### 8.2 v6.1 实际承接

| MVP/v5.3 规划 | v6.1 实际 | 评价 |
|---|---|---|
| 保留 phase1–4 数据底座 | ✅ 已继承。registry、sidecar、metadata、pseudobulk、fraction features 复用或重建 | 强对齐 |
| 升级 phase5–7 算法层 | ⚠️ 部分升级。模块发现、机制裁定、扰动映射、外部验证已形成 Phase4–7，但未进入 causal GRN/SS-CGF | 中等对齐 |
| 新增 shared+specific causal GRN | ❌ 未实现。Phase4 是 role-aware module discovery，非因果图学习 | 弱对齐 |
| 新增 M-signals v1（FG/mICS/MRI） | ❌ 未完整。mechanism/priority scores 可视为雏形，但正式 mIMS/MRI/mICS 文件未生成 | 弱对齐 |
| 新增 ICSmini / e_X-mini / Tier ranking | ⚠️ 靶点预筛有，组合排序弱，正式文件未生成 | 中等对齐 |
| 新增论文级验证框架 | ⚠️ 有框架性输出，但外部锚点覆盖不足、LOCO 未系统报告 | 中等对齐 |
| 后置 Flow/ODE | ✅ 未启动，符合规划 | 强对齐 |

### 8.3 MVP 差距是否被 v6.1 填补

| MVP 差距 | v6.1 填补状态 | 说明 |
|---|---|---|
| 核心算法 SS-CGF / shared+specific causal GRN | **否** | v6.1 仍未进入因果网络学习 |
| 机制层以手工程序为主 | **部分** | Phase4/5 数据驱动，但非因果 |
| 验证体系停留在 demo 级 | **部分** | 新增了 adversarial review、claim audit、负对照，但 LOCO/外部锚点仍弱 |
| 范围收缩到 HCC demo | **部分** | v6.1 三宇宙并行，但 HCC 覆盖受 GSE229772 缺失限制 |
| M-signals 未正式化 | **否** | 正式 mIMS/MRI/mICS 文件仍未生成 |
| TCR 层缺失 | **部分** | 主表仍缺失；增补层 GSE236581 TCR 原始克隆表已建（support_only，未并入） |

---

## 9. 关键风险汇总

### 9.1 高风险

| 风险 | 影响 | 说明 |
|---|---|---|
| **PD1 anchor 主监督被禁止（Gap T1/S1）** | 整个 PD-1-centered nomination 框架失去第一锚点；论文核心主张不可防守 | 规划最大计划外降级 |
| **mIMS/MRI/mICS 文件不存在（Gap T6）** | M-signals 是 v6.1 关键输出；当前无正式文件可引用 | 手稿 Methods 无法引用实际文件 |
| **空间裁判为规则赋标（Gap T4）** | Fig4 的"空间裁判"不能声称定量空间验证 | 核心图的 claim 需降级 |

### 9.2 中风险

| 风险 | 影响 | 建议 |
|---|---|---|
| **Causal GRN 缺失（Gap S2）** | 论文只能走"稳健模块+机制解释"路径 | 明确命名为 inferred module，不升级 causal claim |
| **TCR 主表/整合缺失（Gap T7）** | T 细胞克隆证据链未进主证据 | GSE236581 TCR 原始克隆表已建（1.1M 记录，support_only）；剩余：leakage/split 审计 → 构建 `tcr_clone_master_v6_1.csv` → 并入证据矩阵 |
| **外部锚点仅 IMbrave150（Gap T5）** | L3 translational nomination 证据链薄弱 | 激活 TCGA/GSE179994/GSE200996；或降级为"探索性支持" |
| **特征完整性不均（Gap T2）** | HCC 宇宙 cohort 集中度风险 | GSE120575/GSE229772 rescue 已在 addendum 中推进 |
| **L1/L2/L3 未正式对标（Gap S5）** | 投稿后审稿人将质疑证据等级 | 生成 `candidate_module_evidence_matrix.csv` 与显式降级规则 |

### 9.3 低风险/已控制

- Claim boundary 控制良好（0 forbidden claims）
- 禁止临床推荐用语已落实
- MVP 手工 MOA 未进入主结论
- adversarial review 无 fatal blocker
- 合同路径漂移已记录（path_alignment_report.md 存在，待执行 symlink）

---

## 10. 救援工作包（Rescue Work Packages）

基于 `V6_1_PLAN_ROADMAP_GAP_AND_RESCUE_REPORT_20260615.md` 定义；**进度列为二次修订新增**，反映增补层已预铺设的基础设施：

| WP | 内容 | 优先级 | 退出门 | 进度（二次修订） |
|---|---|---|---|---|
| **WP0：合同加固** | 为当前版本加顶层 "internal-review / conditional" 警告；标记 PD1_anchor 所有下游输出为 support-only；标记 Phase3.5 输出为 invalid | **立即** | 手稿各章节无 PD1_anchor 主监督、无因果证明、无药物推荐 claim | 基本落实（claim audit 0 违规） |
| **WP1：PD1 anchor 升级** | 识别 2+ 个独立干净 PD-1/PD-L1 预处理队列（responder + non-responder 均有）；重跑 Step1/Step2/Step3 anchor 分支；添加 within-cohort meta-analysis | **P0** | PD1_anchor 混杂从 HIGH 降至 MEDIUM 以下，within-cohort 方向一致 | **未启动**（GSE120575_rescue 为黑色素瘤 PD1X，非干净 monotherapy anchor，不满足） |
| **WP2：表达特征救援** | GSE120575、GSE229772、GSE225063 等的 raw-count → gene-symbol → matrix orientation 完整审计；重跑 Step2.8/2.9；再生 hotfix 矩阵 | **P1** | 无零表达主 cohort；signature/pathway/TF 可用性改善无泄漏 | **基础设施部分完成**：GSE120575/GSE229772 pseudobulk 已产出（增补层），待回灌主线 |
| **WP3：定量空间裁判** | 将模块基因/signature 映射到 spatial spots；计算 T cell exclusion/infiltration、myeloid-tumor 邻近、APC niche 邻近、LRG 共定位；添加置换控制 | **P1** | ≥1 核心机制通过预定义方向 + 零值控制的定量空间支持 | 数据源扩展：增补层 GSE238264 第二空间 manifest（support_only）；定量模型未做 |
| **WP4：外部锚点激活** | 激活 TCGA/GTEx/ICGC/GSE179994/GSE200996；使用模块级投影 + 置信区间 + 失败案例 | **P1** | 核心机制在 ≥1 个独立外部锚点复现方向 | **未启动**（仍仅 IMbrave150） |
| **WP5：扰动+TCR 支持** | 从本地资源（GSE133344/193736/306429/90063/l1000）构建 perturb-prior 矩阵；从 GSE236581 构建 `tcr_clone_master_v6_1.csv`；生成正式 mIMS/MRI/mICS 文件 | **P2** | candidate evidence matrix 含 perturb/TCR 字段和降级规则；mIMS/MRI/mICS 文件存在 | **TCR 分支部分完成**：GSE236581 TCR 原始克隆表已建（1.1M 记录）；perturb 矩阵、TCR master、M-signals 文件均未做 |
| **WP6：投稿复现性** | 补 Supp. Table S1（精确计数）；补 Phase4/Step3 参数附录；补软件环境锁/包版本；re-run claim audit | **P0（外部投稿前）** | Phase10 external_submission_ready 变为 true | **未完成**（3 个 Methods 阻塞项仍在） |

---

## 11. 建议

### 11.1 立即行动（外部投稿前）

1. **执行 WP0 合同加固**：防止在救援期间过度 claim。
2. **解决 WP6 三个 Methods 阻塞项**：Supp. Table S1、参数补充、环境锁定。
3. **明确 PD1 anchor 降级声明**：在手稿 Abstract/Discussion 中写明 PD1 anchor 为 weak/auxiliary biological context，非主监督信号。
4. **修复合同路径漂移**：建立 `results/v6_1/step0/` symlink 或更新 analysis_contract.md，使 no-bypass 规则可审计。
5. **生成 `candidate_module_evidence_matrix.csv`**：为每个候选显式对应 L1/L2/L3 字段与降级规则。

### 11.2 近 30 天优先行动（Gap Report 推荐）

1. 启动 **WP1** PD1 anchor 队列搜索/救援（最高战略优先级）。
2. 启动 **WP2** GSE120575 和 GSE229772 表达特征救援。
3. 为四个核心机制实施真实空间评分（**WP3**）：
   - `CAND_P5MG_HCC_012`
   - `CAND_P5MG_PD1X_013`
   - `CAND_P5MG_PD1X_014`
   - `CAND_P5MG_SHARED_021`
4. 激活 ≥1 个 IMbrave150 以外的外部 HCC/bulk 锚点（**WP4**）。
5. **将增补层已建产物回灌主线**（二次修订重点，低成本高回报）：把 GSE120575/GSE229772 pseudobulk 接入 Step2.8/2.9 → 再生 hotfix（WP2 收尾）；对 GSE236581 TCR 原始克隆表做 leakage/split 审计后构建 `tcr_clone_master_v6_1.csv`（WP5 收尾）。这些不是从零开始，而是把已堆好的钢筋浇进结构。

### 11.3 短期增强（下一轮迭代）

1. **激活 TCGA/GTEx/ICGC 外部锚点（WP4）**：提升 L3 可信度。
2. **从已建 TCR 原始表构建 `tcr_clone_master_v6_1.csv`（WP5）**：GSE236581 克隆表（1.1M 记录）已就绪，剩余为审计 + 聚合 + 并入。
3. **生成正式 mIMS/MRI/mICS 文件（WP5）**：这些是规划的 Step8 核心交付物，目前全项目零文件。
4. **整合已建空间 manifest**：GSE238264 addendum manifest（support_only）已落地，下一步是做定量空间评分而非仅再找数据集（WP3）。

### 11.4 中期方法学升级（若目标是方法学顶刊）

1. **落地 shared+specific static GRN 最小版**：从 MVP 到 v5.3/v6.1 的方法学分水岭。
2. **正式化 M-signals v1**：定义 FG proxy / mICS / MRI 稳定 schema 与 patient-level 输出文件。
3. **实现 ICSmini / e_X-mini / Tier ranking v1**：补齐 v5.3 组合排序框架。
4. **建立统一 evaluator**：LOCO、leave-center-out、ablation、calibration、baseline comparison 一键可复现。

### 11.5 论文策略建议

**若 WP1 PD1 anchor 救援失败**，安全投稿路线为：
> "A gated multi-cohort single-cell framework identifies robust HCC/ICI and pan-cancer immune mechanism candidates, with associative external and support spatial evidence, and defines validation-ready hypotheses under strict claim boundaries."

**若 WP1 PD1 anchor 救援成功**，可回归原始论文主张：
> "A PD-1 monotherapy anchor defines base ICI sensitivity and resistance modules; PD-1+X cohorts and external/spatial/perturbation evidence nominate residual TME repair mechanisms."

当前最合理的近期投稿定位：
> **"稳健免疫模块 + HCC 优先转化提名 + 支持性空间证据"**

而非：
> 方法学创新（causal GRN / SS-CGF / 组合排序 / PD-1 中心主叙事）

---

## 12. 结论

### 12.1 总体判断

**v6.1 项目的科学目标与工程目标在框架方向上高度一致，在治理工程上执行优秀，但在科学主链深度和关键方法学交付上存在显著缺口。**

- **框架一致**：PD-1-centered continuum、强基线前置、空间裁判、机制提名而非临床推荐等核心科学原则均被遵守；工程上建立了完整的分阶段门控流水线、数据治理、特征工程、模块发现、手稿包。
- **主链断裂**：PD1 anchor 主监督建模被禁止，mIMS/MRI/mICS 正式文件未生成，causal GRN 未实现，空间裁判为规则化而非定量，TCR 层缺失，外部锚点覆盖有限。这些共同导致项目从"PD-1 中心因果机制提名框架"收敛到"多队列关联证据包 + HCC 转化叙事"。

### 12.2 对齐性评分（修订）

| 维度 | 评分（满分 10） | 说明 |
|---|---|---|
| 科学目标方向对齐 | 7.5/10 | 核心原则一致；但 PD-1 中心连续链第一环断裂 |
| 科学目标深度对齐 | 4.0/10 | PD1 anchor 降级、GRN/M-signals 缺失、TCR 缺失、空间定量缺失 |
| 工程目标流程对齐 | 8.5/10 | 分阶段门控、治理、输出契约极优秀 |
| 工程目标交付对齐 | 5.5/10 | 合同路径漂移、强基线不完整、mICS/MRI/mIMS 文件不存在、外部锚点弱 |
| MVP/v5.3 承接对齐 | 5.0/10 | 数据底座继承良好，算法核心（causal GRN/M-signals）未升级 |
| **总体对齐度** | **6.0–6.5 / 10** | 内部评审就绪，外部投稿有条件，PD-1 中心主张需要救援才能成立 |

*注：相较之前 6.5–7.0 分估计略有下调，反映 Gap Report 对空间裁判（规则化）和 M-signals（文件不存在）的更精确评估。*

### 12.3 最终结论

当前项目已经形成一个**内部评审就绪、外部投稿有条件（CONDITIONAL_PASS）**的 v6.1 候选包，从 MVP 的"演示骨架"升级为"有门控审计证据的 ICI 关联分析包"。

**但项目尚未达到 v6.1 原规划的核心主张**：PD-1-centered nomination framework 在第一科学锚点处失效，mIMS/MRI/mICS 正式 M-signals 体系未落地，空间裁判未达到定量 adjudicator 要求，TCR 证据链缺失。

正确的救援路径不是增加建模复杂度，而是：
1. 补充干净的 PD1 anchor 数据（WP1）；
2. 修复表达特征缺失（WP2）；
3. 实现定量空间/外部验证（WP3/WP4）；
4. 生成正式 mIMS/MRI/mICS 输出文件（WP5）；
5. 降级当前不可防守的 claim（WP0）。

若这些门通过，v6.1 可回归完整 PD-1 中心路线图。否则，应以"审计多队列 ICI 免疫模块候选框架 + HCC 转化提名"路线发表，并诚实处理所有降级 caveat。

**二次修订补记**：好消息是这条救援路径并非从零开始——2026-06-11 增补包已在 `support_only` 状态下预铺设了 WP2（GSE120575/GSE229772 pseudobulk）、WP5-TCR（GSE236581 1.1M 克隆原始表）与 WP3 第二空间数据源（GSE238264 manifest）的基础设施。剩余工作主要是**跨治理门的回灌与整合**（构建 TCR master、回灌表达特征、做定量空间评分、生成 M-signals 文件），而非重新采集与处理数据。坏消息是最关键的 WP1（干净 PD-1 monotherapy anchor 队列）与 WP4（IMbrave150 以外外部锚点）仍未启动，主链第一锚点与 L3 外部证据链依旧是项目的真正瓶颈。在这两者未解决前，对齐评分与可防守 claim 维持本报告第 7 节与第 12.2 节判断不变。

---

## 附录 A：参考文件清单

| 文件 | 作用 |
|---|---|
| `/home/huyudi/006/docs/plan/v6_1/plan.md` | v6.1 科研方案（最新权威规划） |
| `/home/huyudi/006/docs/plan/v6_1/roadmap.md` | v6.1 分步施工总览 |
| `/home/huyudi/006/results/v6_1/V6_1_PLAN_ROADMAP_GAP_AND_RESCUE_REPORT_20260615.md` | **核心 Gap 分析（2026-06-15，含 T1-T8 + S1-S5 + WP0-WP6）** |
| `/home/huyudi/006/results/v6_1/V6_1_PHASE_SYSTEMATIC_ASSESSMENT_20260611.md` | Phase 系统评估（2026-06-11，详细技术基础） |
| `/home/huyudi/006/results/v6_1/INDEX.md` | 当前项目总索引（2026-06-15）|
| `/home/huyudi/006/docs/v6_1/scope_lock_v6_1.md` | v6.1 范围锁定 |
| `/home/huyudi/006/docs/v6_1/analysis_contract.md` | 分析输入输出合同 |
| `/home/huyudi/006/docs/v6_1/PD1_anchor_decision_note.md` | PD1 anchor 决策 |
| `/home/huyudi/006/docs/v6_1/path_alignment_report.md` | 合同路径漂移 |
| `/home/huyudi/006/.codegraph/swarm_v6_1_scientific_goals.md` | 科学目标综合整理 |
| `/home/huyudi/006/results/v6_1/phase10_submission_release_candidate_20260615/PHASE10_FINAL_REPORT.md` | Phase10 最终报告 |
| `/home/huyudi/006/results/v6_1/data_completeness_strengthening_20260611/DATA_COMPLETENESS_STRENGTHENING_20260611.md` | 数据完整性强化评估 |
| `/home/huyudi/006/results/v6_1/new_data_integration_addendum_20260611/ADDENDUM_COMPLETION_REPORT_20260611.md` | 新数据增补完成报告 |
| `/home/huyudi/006/mvp/MVP_GAP_ANALYSIS_2026-04-14.md` | MVP 与 v5.3 差距 |

---

## 附录 B：术语对照

| 术语 | 含义 |
|---|---|
| SS-CGF | Shared–Specific Causal Graph Field（未在 v6.1 实现） |
| mIMS | module immune mechanism score（主读数 1，正式文件未生成） |
| MRI | mechanism restoration index（主读数 2，正式文件未生成） |
| mICS | module intervention candidate score（主读数 3，正式文件未生成） |
| ESR | edge / evidence support rate（辅读数） |
| LRG-spatial | ligand-receptor gain with spatial support（辅读数） |
| FG | Flow Gain / 流增益（已暂时下线） |
| ICSmini | Intervention Candidate Score mini（v5.3 遗留，v6.1 未正式实现） |
| e_X-mini | Expected combination effect mini（v5.3 遗留，v6.1 未正式实现） |
| LOCO | Leave-One-Cancer-Out |
| WP0–WP6 | 救援工作包（见第 10 节） |
| Gap T1–T8 | 技术差距（见第 6.2 节） |
| Gap S1–S5 | 科学差距（见第 5.2 节） |
