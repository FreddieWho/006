# 外部补充数据与阻塞点审计（2026-08-23）

## 结论先行

本轮完成了本地候选重核验和公开数据的 metadata-first 筛查。结果是：

1. **没有发现一个现在就能同时解除阻塞点 1 和 2 的新队列。** 现行 G1 门槛没有降低，也没有把不同终点、不同治疗轴或 bulk 数据直接混在一起。
2. **TASK01 仍是唯一当前可直接进入 WP6B 数字可行性审计的主锚**：T0/pre 为 25 位患者，13 responder、12 non-responder；其 72 个 sample×timepoint 键已经与 `~/004` 恢复出的源元数据逐一对齐。
3. **最值得先修的是本地 GSE301741 的数据接口**：原始队列有 16 位患者，但冻结宽矩阵只保留 2 位可直接使用；后续投影虽得到 11 位 support 患者，只有 4 个桥接配对且不同 frozen module 的桥接不稳定，因此尚不能升级。
4. **外部数据能补的是“验证/桥接”而不是立即替代主锚**：GSE120575、GSE123813 可作为单细胞跨癌种 support；GSE78220、GSE91061 是 bulk-to-module bridge 候选，必须先做预注册的跨模态桥接。

## 两个阻塞点的准确含义

### 阻塞点 1：缺少第二个可比较的独立 response environment

这里的“第二个环境”不是再找一个有 responder/non-responder 标签的数据集，而是要同时满足：

- 独立患者和独立队列；
- 人体肿瘤、治疗前取样，并能证明是在第一剂 PD1/PD-L1 前；
- 单药或可完全拆分的治疗臂；
- 预先冻结的同一类 response 定义和观察窗口；
- 可以在不看 response 的情况下计算 8 个 frozen FM；
- 两个 response 类别都存在，不能由极端分离或标签缺失“制造”效果。

当前形式化门槛是 `n>=15` 且最小类别 `>=5`；实际要做稳健外部复制，建议目标提高到总样本约 30、少数类约 10。这个建议是设计上的稳健性目标，不是偷偷修改当前门槛。

### 阻塞点 2：没有重复的同一 endpoint × treatment × baseline 轴

阻塞点 2 比“再找一个癌种”更严格。必须在至少两个独立环境中重复同一估计单元，然后才可分别估计环境效应并讨论 shared direction。癌种变化可以是异质性来源，但不能用癌种、药物、终点、时间点同时变化来假装复制。

因此，单独增加一个 HNSCC pathologic-response 队列，可能帮助阻塞点 1，却不能自动解除阻塞点 2；单独增加一个 bulk melanoma 队列，也不能替代 scRNA 的直接主分析。

## 本地数据核验

### TASK01：源链已恢复，仍是“条件性主锚”

当前绑定表给出：

| 项目 | 当前事实 |
|---|---|
| T0/pre 患者 | 25 |
| responder / non-responder | 13 / 12 |
| endpoint | mRECIST |
| treatment | PD1 |
| T1/T2 | 仅 support / temporal sensitivity |
| 当前冻结边界 | `temporal_sensitivity_only` |

脚本 `scripts/preprocessing/script_sc/task01.R` 指向的 `~/004.1` 路径已不存在，但在 `~/004/_from_004.1/zz.META/` 找到了源元数据镜像。`task01_raw.csv` 有 73 行，其中 72 行有效；与当前 006 表的 72 个 sample×timepoint 键在 `WYZ→WYZ1` 别名处理后完全一致，T0/T1/T2 数量均为 25/24/23。response 表中当前源患者均能匹配，另有 JGY/LW/ZJF 三条不在当前样本源中的记录，未被纳入。样本目录还保留 `Pre`、`Tace`、`ICI` 标签，部分样本名含日期片段；但文件仍没有每位患者的首剂 PD1 日期、正式活检日期或 response assessment 日期。

所以需要向导师明确的一句话是：

> **TASK01 的 T0 选择现在已经有源元数据支撑，25 例主锚可以继续保留；但现有文件仍不能逐患者证明 T0 早于第一剂 PD1，因此只能写成 source-backed conditional baseline，不能写成 provenance-clean baseline。**

这不是要推翻 TASK01，也不是要重新选择患者。要把 `temporal_sensitivity_only` 升级为正式 baseline，需要补齐或核实：

1. T0 活检早于第一剂 PD1；
2. response 是治疗后的 mRECIST 结果；
3. T1/T2 没有被混入 baseline 特征；
4. 证明后冻结一版新的 baseline provenance manifest，再重新做 metadata-only gate。

源链核对明细见 [TASK01 source reconciliation](/home/huyudi/006/results/v6_2/scientific_engineering_realignment_v1/task01_source_reconciliation_v1.yaml)。在日期证据补齐前，TASK01 可用于当前的条件性数字审计，但不宜把“完全无时间混杂”写成已经验证的结论。

### GSE206325：不是新增主锚

本地有效绑定显示它目前只有 post-treatment 样本，endpoint 为 `unknown_or_mixed`，`supervised_use_allowed=no`。这与旧候选注册表中“clean anchor”的文字标签冲突；以可执行绑定表为准，不能把它加入 baseline 主分析。

### GSE301741：首要上游修复对象

公开 [GSE301741 GEO 记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741) 是 16 位 HNSCC 患者的 pembrolizumab 新辅助前/后单细胞/TCR 队列。当前项目的损失轨迹是：原始 16 位 → phase4a 交接 5 位 → phase7 共同可用 4 位；冻结 WP6B 宽矩阵中的直接 baseline 只有 2 位。后续 support 投影得到 11 位（6 responder/5 non-responder），但只有 4 个 mid bridge pairs，FM01/04/07 的桥接相关不稳定，故仍为 `support_only_unresolved`。

这说明阻塞主要是**数据接口/表示桥接损失**，不是原始队列不存在。修复优先级高于继续盲目搜更多 accession：先找出 phase4a→phase7 丢失的 22 个样本，并逐个验证 frozen FM 的基因映射、样本 ID、coverage 和 module score 分布。

### GSE286827：治疗臂拆分后仍不够

Durvalumab 单药 pre arm 当前为 13 位（4 responder/9 non-responder），低于 `n>=15`；durvalumab+tremelimumab 为组合治疗，只能 support。即使修正注册表的臂级标签，也不能单独解除两个阻塞点。

## 公开数据筛查结果

### GSE120575：有标签，但不能作为当前 clean primary

官方 GEO 记录该研究包含 48 个肿瘤样本、19 个治疗前样本，并提供患者 ID、时间点、治疗和 response 元数据。解析后，治疗前总计 19 位患者（9 responder/10 non-responder）；其中 anti-PD1 单药只有 12 位（4/8），其余混有 CTLA4 单药或 CTLA4+PD1 联合。它适合做外部 support 或跨模态 bridge 候选，但按当前门槛不能成为 clean primary。

官方来源：[GSE120575 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575)。

### GSE123813：跨癌种价值高，但样本太小且标签需外接

官方 GEO 是 BCC 的 anti-PD1 前/后单细胞/TCR 数据；公开 cell metadata 没有 response 列。原始论文报告 11 位患者，其中 6 responder、5 non-responder。它能作为 BCC 的独立验证候选，但低于 `n>=15`，且必须先把论文标签与稳定患者 ID 对齐。

官方来源：[GSE123813 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813)；response 分组需以对应原始论文的患者表为准。

### GSE78220：最适合做 bulk bridge，不是 scRNA 主锚

官方 GEO 提供 28 个 melanoma pre-treatment bulk RNA-seq 样本的 FPKM workbook；公开分析将其按 15 responder/13 non-responder 使用（例如 [TIP 案例说明](https://aacrjournals.org/cancerres/article/78/23/6575/543748/TIP-A-Web-Server-for-Resolving-Tumor)）。样本数足够做候选验证，但它不是单细胞，不能直接替代 TASK01；必须先冻结 gene-set 级 FM 映射、归一化和 response-blind bridge 规则，再决定能否作为独立方向复制。

官方来源：[GSE78220 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220)。

### GSE91061：目前最强的 bulk bridge 候选

官方 GEO 设计为 109 个 bulk RNA-seq 样本、65 位患者，其中 51 个 pre-treatment、58 个 on-treatment；官方 family SOFT 的样本字段给出 pre-treatment 中 10 个 PR/CR、39 个 SD/PD、2 个 UNK。原始论文将该队列定义为 nivolumab 治疗的 melanoma 队列。它因此是目前最强的 bulk bridge 候选，但仍不能把 51 个 pre 样本直接算成 scRNA clean primary：必须先冻结 bulk-to-FM 的 gene-set 映射、归一化和 response-blind QC。

官方来源：[GSE91061 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061)；治疗信息见原始 [nivolumab 论文记录](https://pubmed.ncbi.nlm.nih.gov/29033130/)。

## 这批数据对两个阻塞点实际贡献什么

| 候选 | 对阻塞点 1 | 对阻塞点 2 | 当前处理 |
|---|---|---|---|
| TASK01 | 已提供一个条件性主锚 | 不能提供独立重复 | 保留，先补时间线 |
| GSE301741 | 若修复交接和 FM bridge，可提供第二癌种 | pathologic endpoint 与 mRECIST 不同，不能直接重复 | 上游修复优先 |
| GSE286827 单药臂 | 样本不足 | PD-L1/pathologic 轴不同 | validation/support |
| GSE120575 | anti-PD1 mono 子集不足，混合治疗 | 终点和跨模态桥接未冻结 | external bridge |
| GSE123813 | 样本不足，标签需外接 | 仅 support，不能同轴复制 | label join + support |
| GSE78220 | bulk 数量足，但不是直接 scRNA primary | 只有完成跨模态 bridge 才能谈方向复制 | bridge candidate |
| GSE91061 | 51 例 pre，10/39/2 标签，适合 bridge | 仍需跨模态桥接，不能直接成为 scRNA 同轴复制 | freeze bridge + response-blind QC |

## 建议的解锁顺序

### 第一优先级：修复 GSE301741 上游接口

这是最可能用现有数据增加一个可比较环境的点。验收条件不是“投影出了更多患者”，而是：

- 22 个丢失样本的每一步原因可追溯；
- 至少 15 位 baseline 患者重新进入同一 frozen representation；
- 8 个 FM 都完成 gene mapping、coverage、score 分布和 response-blind QC；
- 4 个桥接配对扩充并且模块方向不再出现明显不一致；
- role、endpoint、timepoint 和 support/primary 状态在任何 response effect 前冻结。

### 第二优先级：恢复 TASK01 临床时间线

向数据提供者索取三类最小字段：`patient_id`、`biopsy_collection_date`、`first_PD1_dose_date`，最好再有 `response_assessment_date`。只要这三类字段能逐患者对齐，就可以把 T0 的“操作性 baseline”升级为 provenance-clean baseline，或明确剔除少数不合格患者。

### 第三优先级：建立外部 bulk bridge

先用 GSE78220，再评估 GSE91061（后者现在已具备样本级 response 字段）。步骤必须是：

1. 只用冻结的 8 个 FM 基因集，不按 response 选基因；
2. 预先固定 bulk normalization、缺失基因处理和 score aggregation；
3. 在不读取 response 的情况下完成样本/队列 QC；
4. 再冻结治疗臂和 endpoint；
5. 最后才读取 response 做独立方向性验证。

### 第四优先级：只把小队列用于支持，不把它们“拼大”

GSE123813 和 GSE286827 单药臂可以检验方向是否大致一致，但不能因为样本少就与 TASK01 直接合并，也不能用跨癌种、跨终点的合并结果宣称解除阻塞点 2。

## 本轮授权与边界

本轮已按用户明确授权创建外部候选筛查解锁：

`results/v6_2/scientific_engineering_realignment_v1/external_anchor_search_unlock_v1.yaml`

该解锁允许读取公开 response/endpoint 元数据、下载公开处理文件、做候选可行性计数和来源审计；**没有**解锁 response-driven module selection、effect model、shared direction、WP7 或 Phase8。审计明细在：

`results/v6_2/scientific_engineering_realignment_v1/external_anchor_search_audit_v1.tsv`

G0/G1 冻结文件未修改。当前状态仍是：`WP6B_STOPPED_FEASIBILITY`、`WP7 NOT_STARTED_BLOCKED`、`Phase8 NOT_RUN`。

## 需要导师明确的 TASK01 决策

建议导师只需确认以下最小决策：

> **在原始给药/采样日期尚未恢复前，是否同意把 TASK01 T0 继续作为“条件性主锚”，并把“时间线 provenance 补齐”列为下一步上游任务，而不是现在强行把它改写为完全 clean baseline？**

我的建议是“同意”。这样既保留现有 25 例的研究价值，又不把 T0 标签过度解释成已证明的首剂前采样；同时把主要精力放到 GSE301741 接口修复和外部 bridge，而不是继续无边界地搜 accession。
