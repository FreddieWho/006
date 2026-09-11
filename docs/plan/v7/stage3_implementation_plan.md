# v7 Stage 3 实施方案：空间基础设施与最小跨平台回放

状态：**COMPLETE_WITH_LIMITATIONS**（D3 通过，2026-09-04）

版本日期：2026-09-03

适用阶段：v7 Stage 3

上游基线：Stage 1 数据注册表；Stage 2 共同词汇、计量合同与 D2 结果

执行结果：已完成 8 个逻辑 pilot 单元下的 11 个真实技术捕获（7 个患者身份包，Visium 与 Xenium 两类平台）。counts、坐标、患者分组、Stage2 投影、section-local 图、99 次 section 内置换和患者级 bootstrap 均已落盘；1 个共同特征在两类平台均完整可测并产生非退化空间统计。D3 为 `D3_PASS_WITH_LIMITATIONS`，Stage 4 可消费。Xenium 上游 h5ad 的 `uns/response` 仅在响应盲副本中移除，原文件未修改；GT 只完成 locator/hash 封存审计，未解析标签值。

## 1. 本阶段要解决什么

Stage 3 不直接寻找“免疫屏障”、TLS 或 PD-1 疗效机制，而是先回答：

> Visium、Xenium 和既有 h5ad 等不同空间数据，能否在保留各自技术差异的情况下，被整理成同语义、可复核、可汇总到患者层的空间测量？

本阶段完成后，Stage 4 才能在可靠的计数、坐标、患者身份和空间统计上做 response-blind 结构发现。Stage 3 本身只形成方法与资源层结论，不形成疗效预测、共享空间机制或潜在场数量结论。

## 2. 与前两阶段的接口

### 2.1 从 Stage 1 继承

- 数据集、患者、组织块、切片和技术捕获之间的身份关系；
- 数据角色、重复血缘和可用性状态；
- 空间数据的真实文件位置与 metadata provenance；
- validation 数据的只读角色，不因本阶段接入而改变。

### 2.2 从 Stage 2 继承

- `module_dictionary.tsv` 和 `module_membership.tsv`；
- cell-state ontology；
- gene identifier 解析和 score 计算语义；
- “缺失基因保留为 NA，不补 0”的观测掩码规则；
- coverage、technical status、uncertainty 和来源记录字段。

Stage 2 的 `D2_class` 只代表 scRNA 计量证据。空间结果必须新增 `spatial_projection_status`，不能把 scRNA 的可测性直接冒充为空间平台的可测性。

### 2.3 不从 `/013_spatial` 继承的内容

可以复用其输入政策、患者分组、计数读取和 provenance 思路，但不得继承：

- 未完成的真实数据 `K` 选择；
- K=0/K=3 wiring diagnostic 的科学结论；
- 被 CPU 运行时间阻塞的 K=2 bridge 结果；
- 任何旧 checkpoint、潜在场、TLS 或屏障结构结论。

`/013_spatial` 在本阶段保持只读；Stage 3 的可执行代码和结果在 `/006` 内自包含。

## 3. 最小回放面板

首轮采用 8 个高价值空间单位，覆盖两类平台、三种输入形态和三个癌种。

| 数据集 | 癌种 | 平台/输入 | 计划单位 | 本阶段用途 | 明确不做 |
|---|---|---|---:|---|---|
| HTAN CRC | CRC | spatial h5ad | 4 个技术捕获 | h5ad adapter、同患者多捕获嵌套、已知结构锚点的隔离测试 | 不用 GT 选特征或调统计方法 |
| GSE238264 | HCC | Visium 10x H5 + positions + scalefactors | 2 个患者单位 | 原生 Visium 文件链、坐标/尺度/图像元数据回放 | 不读取 R/NR、疗效或治疗后机制 |
| GSE291246 | BCC | Xenium cell-level h5ad | 2 个不同患者单位 | cell identity、panel gene 缺失、segmentation 引用与 Xenium 回放 | 不读取 pre/post 或疗效，不做 HCC 证据 |

HTAN 暂定使用 `7003_AS_1`、`7003_AS_2`、`7003_AS_3` 和 `6723_KL_1`。其中前两个用于检验同一患者/组织块下多个技术捕获不能被误计为独立患者；若 metadata 不能证明它们是独立生物切片，统一称为“技术捕获”。

GSE238264 和 GSE291246 的具体单位在 S3.0 由代码确定：先按纯技术资格过滤，再对不含 response、treatment、timepoint 的患者键做稳定哈希排序。选择过程不能检查文件名中的 R/NR 或 pre/post 字样；实际入选单位与选择依据写入 manifest。

以下数据继续保持锁定，不参加 Stage 3 方法选择：

- GSE175540/Meylan ccRCC：external validation；
- USZ TLS：external validation；
- ST_CRC_CMS：internal/serial-section validation。

已知结构锚点只要求“存在、可定位、被密封”。Stage 3 仅在 discovery 输出冻结后，以独立 evaluator 检查 GT join 和隔离机制，不检验生物学富集。

## 4. 数据防火墙

### 4.1 双 manifest

建立两个互不混用的清单：

1. `model_safe_input_manifest.tsv`
   - 只包含 counts、spot/cell/bin ID、坐标、必要尺度、非 GT tissue mask 和技术元数据；
   - 不包含 response、endpoint、治疗结局、结构 GT、GT distance 或 GT-derived mask。
2. `sealed_validation_manifest.tsv`
   - 只保存 GT locator、文件哈希、允许的验证任务和角色；
   - 主回放 pipeline 不得导入或读取此文件。

### 4.2 禁止进入发现输入的内容

- response、endpoint、疗效分组和由其派生的字段；
- TLS、tumour boundary 等结构 GT 及其距离、mask、人工区域；
- 当 GT 由 H&E 定义时，H&E-derived predictive features；
- 患者、中心、批次、治疗或文件名编码本身作为预测特征；
- 已归一化、插补或推断表达冒充 raw counts。

患者/中心/平台信息仍可用于分组、QC、分层汇总和技术审计，但不能成为生物学 score 的输入。

### 4.3 可验证的隔离规则

- 主 pipeline 只接受 model-safe manifest；
- GT evaluator 放在独立模块，主 pipeline 无 import path；
- discovery 输出及其哈希冻结后才能启动 GT join 测试；
- 修改 sealed GT 文件不应改变 discovery 输出哈希；
- 日志、异常和表格中不得泄露被禁止字段的取值。

## 5. 软件结构

计划新增：

```text
config/v7/stage3.yaml
config/v7/stage3_pilot.tsv
scripts/v7/spatial/run_stage3.py

src/v7/spatial_io/
  contracts.py
  registry.py
  locators.py
  gene_mapping.py
  gt_guard.py
  qc.py
  materialize.py
  adapters/
    base.py
    h5ad_counts.py
    visium_10x.py
    xenium.py

src/v7/spatial_stats/
  scores.py
  graphs.py
  geometry.py
  autocorrelation.py
  adjacency.py
  composition.py
  nulls.py
  patient_summary.py
  pipeline.py

src/v7/spatial_validation/
  sealed_gt_evaluator.py

tests/v7/
  test_spatial_contract.py
  test_spatial_adapters.py
  test_spatial_identity.py
  test_spatial_gt_isolation.py
  test_patient_block_split.py
  test_spatial_projection.py
  test_spatial_stats.py
  test_stage3_integration.py
```

实现优先使用当前已有的 `numpy`、`pandas`、`scipy`、`h5py`、`anndata`、`networkx`、`shapely` 和 `pyarrow`。不把 `scanpy`、`squidpy`、`libpysal`、`esda` 或 `scikit-learn` 作为必需依赖；当前环境中这些依赖缺失或不可稳定导入，而本阶段统计可由稀疏矩阵和 `scipy.spatial.cKDTree` 实现。

## 6. 统一数据合同

每个空间单位至少输出以下四类表。

### 6.1 Unit 表

```text
dataset_id, opaque_patient_id, opaque_block_id, opaque_section_id,
capture_id, modality, platform, cancer, source_locator_hash,
counts_semantics, coordinate_system, scale_status,
image_status, segmentation_status, tissue_mask_status,
identity_status, adapter_status, qc_flags
```

### 6.2 Observation 表

```text
observation_id, capture_id, native_barcode,
native_x, native_y, native_coordinate_unit,
analysis_x, analysis_y, analysis_coordinate_unit,
in_tissue, segmentation_id, observation_qc
```

原始 barcode 和后缀必须保留。不同切片不得连边；同名 barcode 也不能跨 capture 合并。

### 6.3 Feature projection 表

```text
observation_id, feature_id, native_value, standardized_value,
expected_genes, observed_genes, gene_coverage,
technical_status, source_scRNA_D2_class, D2_scope,
spatial_projection_status, uncertainty, qc_flags
```

### 6.4 Hierarchy 表

```text
opaque_patient_id, opaque_block_id, opaque_section_id, capture_id,
relationship_provenance, relationship_confidence, leakage_group
```

无法解析的身份保持 unknown，不根据表达或空间邻近反推患者关系。

## 7. 分步实施

### S3.0 合同、资源和运行前检查

1. 校验 Stage 1、Stage 2 manifest 和关键结果哈希；
2. 检查 `/013_spatial` 的 sample registry、structure registry、R04 input manifest 和 `input_policy.tsv`；
3. 生成 8-unit pilot 选择表以及选择暴露审计；
4. 建立双 manifest，并对禁止字段做列名、值域和路径扫描；
5. 检查输入文件可读性、空间占用和输出预算；
6. 记录 Python、包版本、Git commit、dirty state、配置哈希和源文件哈希。

运行前硬性停止条件：源文件不可读、身份链断裂、raw count 语义无法判定、model-safe manifest 含有禁止字段、Stage 2 合同哈希不一致或可用磁盘低于现有 R04 的 1.2 TB 安全底线。

当前磁盘约 1.5 TB 可用，只允许引用原始文件，不复制大型原始矩阵或图像。

### S3.1 Adapter 与合成测试

先以小型合成 fixture 固定以下行为，再接真实数据：

- 稀疏矩阵方向、gene/barcode 对齐和整数 UMI 检查；
- 10x H5 与 positions 的一一 join；
- native array、pixel 和可验证 micron 坐标分别保存；
- h5ad layer/X 选择有显式依据，不能自动猜测；
- Xenium panel 未检测基因记为 unavailable/NA；
- 重复 observation、barcode suffix 丢失和跨 section 合并直接失败；
- 禁止 GT/response 输入时 fail closed。

各 adapter 返回同一合同，但保留平台原生字段。平台差异通过状态列和 sidecar 表表达，不用强制填成虚假的共同单位。

### S3.2 真实文件最小读取与 QC

按数据集逐个完成：

- count conservation：读取前后每个 observation 和 feature 的总量关系可解释；
- coordinate audit：坐标有限、方向明确、同一 capture 内有变化；
- scale audit：仅在有 scalefactor 或可信标尺时声明物理距离；
- tissue mask audit：记录来源，不能用结构 GT 替代；
- image audit：主流程只登记文件、尺寸、哈希和 scalefactor，不读像素作为发现特征；
- segmentation audit：仅验证 observation-to-segment join；缺失时标为 N/A，不伪造 segmentation；
- identity audit：patient/block/section/capture 层级闭合，重复单位进入同一 leakage group。

Visium 使用整数 UMI 和 Space Ranger 文件链。Xenium 优先使用 cell-feature matrix；若以后从 transcript 构建矩阵，必须作为新的明确 adapter 路径，不能在当前 adapter 内静默切换。

### S3.3 Stage 2 词汇投影

1. 复用 Stage 2 gene mapping 与模块成员；
2. 按 observation 分块计算，不把完整矩阵复制到结果目录；
3. 对每个 feature 记录平台面板实际可见的成员和覆盖；
4. 缺失成员保持 NA，不做跨平台插补；
5. standardized score 只在同一平台/数据集内使用明确的标准化方法；
6. 输出 section-sharded Parquet 和 unit-level coverage summary；
7. 对 spatial projection 单独给出 `measurable / limited_with_uncertainty / not_estimable`。

本步骤只证明“同一个功能对象如何被各平台观测”，不要求所有平台覆盖相同基因数，也不据此挑选对疗效最有利的模块。

### S3.4 空间图和基础统计

所有图在 section/capture 内构建，禁止跨切片连边。

#### 邻域图

- 使用 `cKDTree` 构建稀疏 kNN 和 radius graph；
- 不生成全距离矩阵；
- 原生距离除以该 section 的中位最近邻距离，形成可跨平台描述的无量纲距离；
- 只有 scale audit 通过时才额外报告微米距离。

#### 连续 score 统计

- 邻域均值与局部差异；
- Moran 型空间自相关；
- 距离分箱 variogram；
- graph edge 上两个 score 的关联；
- 保持 section 几何和观测掩码的空间置换 null。

#### composition、adjacency 和 boundary

- 有可信 cell composition 或 cell label 时计算 local composition 和 adjacency；
- 没有 composition 输入时明确输出 `NOT_ESTIMABLE_NO_COMPOSITION_INPUT`；
- 没有非 GT boundary 时不计算 boundary distance；
- 不用结构 GT 补齐缺失的 composition 或 boundary。

Stage 3 不设生物学显著性门槛，也不按结果方向筛特征。输出效应量、null 分布、可重复性与数据限制，留给 Stage 4 做探索性比较。

### S3.5 患者层汇总与不确定性

- observation → capture → section → block → patient 逐层汇总；
- patient 是外层独立单位，多个 section 不能增加患者样本量；
- 同一患者内的 block/section/capture 作为嵌套重采样层；
- 数据集层汇总时患者等权，避免 observation 多的样本支配结果；
- split 以 leakage group/patient 为单位，用稳定哈希生成，不按 response 分层；
- 输出每层单位数、缺失率、有效患者数和 bootstrap 区间。

此处 bootstrap 用于表达患者与切片采样不确定性，不把 8-unit pilot 当成总体人群验证。

### S3.6 GT 防火墙验证

1. 完成并冻结 S3.1–S3.5 的输出及哈希；
2. 独立启动 `sealed_gt_evaluator.py`；
3. 仅用 HTAN training GT 检查 locator、identity join、允许任务和泄漏防护；
4. 不计算候选结构的生物学富集，不比较模型优劣；
5. 用 sentinel 测试证明 GT 内容变化不会改变主 pipeline 输出；
6. GSE175540、USZ TLS 和 ST_CRC_CMS 的 GT 保持未读。

### S3.7 全面回放、D3 判定与冻结

- 对 8 个单位执行完整 replay；
- 每个 unit 有独立状态文件，可断点续跑；
- 对相同输入和配置重跑一次轻量确定性检查；
- 汇总 adapter、identity、counts、coordinate、scale、coverage、GT isolation 和统计结果；
- 生成 D3 gate 与最终报告；
- 只有 D3 通过后，才把 Stage 3 标为完成并允许 Stage 4 消费结果。

## 8. D3 技术准入标准

D3 只判断数据和计算是否可用于下一阶段，不根据是否得到“漂亮”的生物学结果来放行。

### 8.1 整体最低条件

- 真实完成 6–10 个空间单位，本方案目标为 8 个；
- 覆盖至少 2 类平台、2 个癌种；
- 至少一条 Visium 原生文件链和一条 Xenium 全链路成功；
- 至少一个已知结构锚点被正确密封并通过隔离测试；
- 至少一个 Stage 2 feature 在各主要平台均非退化地可计算；
- patient/block/section/capture 关系和 leakage group 可审计；
- 主 pipeline 未读取 response、endpoint 或结构 GT；
- 输入、代码、配置和输出都有可复核哈希。

### 8.2 单位级状态

- `eligible`：身份、counts、坐标和 GT isolation 均通过；
- `limited_support`：主链可用，但 scale、image、segmentation 或部分 feature 不可测；
- `reject`：核心语义无法确认或存在无法消除的泄漏。

image、segmentation、composition 和物理尺度允许按平台标为 N/A；counts、observation identity、基础坐标和 GT isolation 不允许用替代数据补过。

### 8.3 必须停止并保持 Stage 3 未完成的情形

- 实际成功单位少于 6，或只剩单一平台/单一癌种；
- Xenium 或非 Visium 链路全部失败；
- counts 语义只能靠猜测；
- matrix、barcode 与坐标无法可靠 join；
- response/GT/GT-derived 信息进入主 pipeline；
- locked external/internal validation GT 被提前读取；
- 用 normalized、imputed 或 metadata-only 数据替代核心 counts；
- 没有任何共同且非退化的 feature 可用于跨平台统计。

发生这些情况时，报告具体失败单位、失败原因和对结论的影响；不得缩小到单平台后仍宣称 Stage 3 完成。

## 9. 预期结果文件

```text
results/v7/spatial_foundation/
  model_safe_input_manifest.tsv
  sealed_validation_manifest.tsv
  pilot_selection_audit.tsv
  adapter_qc.tsv
  count_conservation.tsv
  coordinate_scale_audit.tsv
  identity_leakage_audit.tsv
  feature_coverage.tsv
  spatial_projection_summary.tsv
  graph_statistics.parquet
  spatial_null_statistics.parquet
  patient_level_summary.parquet
  split_audit.tsv
  gt_isolation_audit.tsv
  platform_replay_manifest.yaml
  STAGE3_RUN_MANIFEST.yaml
  D3_GATE.json
  SPATIAL_FOUNDATION_REPORT.md
```

本次实际落盘的核心汇总文件与上述清单一致，并额外保留了
`pilot_manifest.tsv`、`technical_eligibility_audit.tsv`、
`source_asset_audit.tsv`、`spatial_variogram.parquet`、
`composition_edge_status.tsv` 和 `patient_bootstrap_ci.tsv`。大矩阵的
observation-level score 仅作为未纳入 Git 的 `scores/*.parquet` 分片；其
文件哈希和生成方式由运行清单记录。

大体积 observation-level score 和 graph edge 分片不提交 Git；manifest 保存路径、大小、哈希和生成命令。代码、配置、汇总表、报告和小型测试 fixture 可以提交。

## 10. 验证策略

这是跨模块行为变更，实施时采用分层验证：

1. **合同测试**：字段、类型、禁止列和状态枚举；
2. **adapter 单元测试**：矩阵方向、barcode join、坐标单位、panel 缺失；
3. **泄漏测试**：GT/response fail closed、跨 section 连边禁止、patient split；
4. **统计测试**：在已知合成空间梯度、随机场和断裂图上检查方向与退化行为；
5. **真实小样本 smoke**：每个平台一个单位；
6. **8-unit integration replay**：完整输出和 manifest；
7. **确定性复核**：相同输入、版本和配置产生一致的汇总哈希。

核心真实数据回放不能由合成测试或 smoke test 替代。若只完成前五项，状态仍为 `STAGE3_IN_PROGRESS`。

## 11. 资源与执行顺序

```text
S3.0 运行前审计
  ↓
S3.1 adapters + synthetic tests
  ↓
S3.2 三个平台各 1 个真实单位 smoke
  ↓
S3.3 Stage 2 feature projection
  ↓
S3.4 spatial stats + nulls
  ↓
S3.5 patient nesting / bootstrap
  ↓
S3.6 sealed GT isolation test
  ↓
S3.7 8-unit replay + D3 freeze
```

本阶段预计 CPU 可完成，不依赖 GPU。优先稀疏、分片和按 unit 处理，避免重复保存 raw data。若 Xenium 单位过大，允许调整 chunk 大小和并行度，但不允许抽样 observation 后替代正式回放。

## 12. 阶段完成后能够支持的结论

若 D3 通过，可以说：

- 已建立覆盖 Visium、Xenium 和 h5ad 的可复现空间输入与统计底座；
- 不同平台的原始计数、坐标、尺度和缺失信息被保留，没有被强行抹平；
- Stage 2 功能词汇可按平台观测范围投影，并携带 coverage 与 uncertainty；
- 空间统计可正确汇总到患者层，且 GT、疗效标签和验证集未进入发现输入；
- 合格数据可进入 Stage 4 的 response-blind 空间结构发现。

仍不能说：

- 已发现跨癌种共同免疫屏障；
- TLS、NET、CAF 或其他结构与 PD-1 failure 有关；
- 某个潜在场数量、GNN 或 architecture surrogate 更优；
- 8-unit pilot 已构成独立临床验证。

## 13. 当前是否需要用户决策

当前没有必须由用户决定的阻塞点。默认采用上述 8-unit 面板、`/013_spatial` 只读、三个锁定验证集不触碰的方案即可开始实施。

只有出现以下情况才需要重新决策：

- 8-unit 面板因实际文件或身份问题不能满足 D3 最低覆盖；
- 必须提前动用 locked validation 数据才能补足平台链路；
- 需要下载新的大型原始数据或显著扩大磁盘/计算预算；
- 发现现有 counts 实际为归一化矩阵，必须更换输入来源。
