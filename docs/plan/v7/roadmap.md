# v7 执行路线图

版本：v7.0
日期：2026-08-31
状态：Stage 0 完成；Stage 1 待启动；尚未运行 v7 科学计算
科学依据：`docs/plan/v7/plan.md`

---

## 0. 路线图定位

v7 的工程目标不是继续修补旧 Phase8，也不是先造一个复杂模型。它要建立一条新的、可审计的主线：

```text
多模态数据与身份接口
  → 共同分子/细胞状态词汇
  → response-blind 空间场与屏障结构发现
  → PD-1 临床锚定和 context 分解
  → HCC 深入分析
  → PD1+X 观测重排
  → perturbation / bulk / 外部组织证据
  → Evidence Card 与 X-class 外推
```

本路线图是探索性研究路线。节点、模型和阈值可以根据证据调整，但数据身份、患者级独立性、模态语义、结果可追溯和结论边界不能被跳过。

---

## 1. v6 继承与废弃规则

### 1.1 直接复用

- `data_interface_repair_v1` 的 patient/sample/timepoint/treatment/response 接口；
- 原始 counts、表达层、QC、gene mapping 和隔离记录；
- v6 response-blind 模块发现与全对象投影基础；
- Phase7 measurement、missingness、route disagreement、coarse merge 和 identifiability 审计；
- TASK01、TASK02、LAMBRECHT_HCC 等患者级/纵向来源，经 v7 重新核验后使用；
- GSE301741 身份交接与 coverage 根因报告；
- 真实 bulk、perturbation、TCR 和 spatial 数据资产；
- negative controls、patient-grouped split、provenance 和 Evidence Card 思路；
- `~/013_spatial` 的 sample/structure registries、患者/组织块切分、ground-truth 隔离、R-04 loader 和 checkpoint 审计。

### 1.2 复用但重新解释

- 8 个 FM：经验分子程序基线，不是 failure barrier；
- FM07：旧条件下的独立归因候选，不是唯一有信息的 module；
- coarse barrier groups：联合表示和不确定性对象，不是删除特征的理由；
- v6 Phase5 混杂结果：风险证据，不是 v7 模型性能；
- R-04：方法与审计资产，不继承任何已成立 latent-field 结论；
- HCC：主要器官 context 和 PD1+antiangiogenic 深入案例，不代表全部癌种。

### 1.3 版本化废弃

- SRB 必须作为核心模型；
- SNGM 只在末端充当 sidecar；
- FM01/FM04/FM07 三模块合同和 FM07-only response 主线；
- identifiability 作为 representation hard filter；
- v6 旧 Phase 编号和含义；
- 预注册、response unlock、固定 gate 等治理语义；
- 历史 Phase8 effect 作为现行结果；
- v6.1 hash-derived perturbation 连续分数；
- 没有多 X-class 时的 zero-shot 强结论；
- 将所有空间数据或所有模态拼成一个统一 embedding。

这些资产不物理删除，进入只读历史归档或资产处置清单。

---

## 2. v7 目录与输出约定

```text
docs/plan/v7/                 # 活跃科学方案与路线图
docs/v7/                      # 当前状态、资产处置、接口说明
config/v7/                    # 数据源、角色与运行配置
scripts/v7/                   # v7 可复现入口，按任务分层
src/v7/                       # 可复用实现模块
tests/v7/                     # v7 接口、泄漏和统计不变量测试
results/v7/
  STATUS.md                   # 当前科学状态
  inheritance/                # v5/v6/013 继承清单
  registry/                   # 多模态数据与患者/切片账本
  ontology/                   # cell state、module、barrier vocabulary
  spatial_foundation/         # 坐标、图像、分割、平台 adapter
  spatial_discovery/          # response-blind fields/circuits
  clinical_anchor/            # 环境内 PD-1 关联
  context/                    # shared/context/HCC 分解
  repair/                     # PD1+X paired/trajectory
  perturbation/               # 真实 perturbation mapping
  external_validation/        # bulk/scRNA/spatial 外部验证
  evidence_cards/             # claim-level 证据与冲突
  figures/                    # 可复现图表
```

大型矩阵、h5ad、图像、raw counts 和模型 checkpoint 保留在本地数据/结果目录，通过 manifest、路径和 SHA-256 引用，不进入普通 Git。

---

## 3. 全局执行原则

1. 患者、样本、组织块、切片、spot/cell、timepoint、treatment 和 response 逐层唯一映射。
2. 同一患者的所有切片和时间点始终在同一 split。
3. 每张空间切片独立构图，禁止跨患者直接连边。
4. counts、normalized、imputed 和 image-derived feature 明确分层。
5. 先跑透明基线，再决定 latent-field、OT、GNN 或深模型是否必要。
6. discovery/calibration 可探索迭代；但对同一 claim，参与选择的数据不得再升级为独立验证，claim-bearing validation 必须保持未触碰并记录角色来源。
7. 无 response 的空间数据可做发现与结构验证，不能冒充疗效标签。
8. 非等价治疗、endpoint 和时间点先分环境估计，跨环境只做分层汇总。
9. 结果同时报告方向、量级、患者级不确定性、异质性和失败案例。
10. 核心计算不能静默跳过、替代或降级。

---

# Stage 0｜v7 Reset、仓库重组与资产继承

状态：**COMPLETE**

## 科学目的

确保 v7 从真实证据起步，而不是继承旧版本的模型名称、阶段语义或未完成结论。

## 任务

1. 建立 v7 plan、roadmap 和当前状态入口；
2. 对 v6 与 `~/013_spatial` 资产做复用/重解释/废弃分类；
3. 历史文档与临时诊断脚本进入可逆归档；
4. 大型数据和结果不搬动，以 manifest 记录逻辑角色；
5. 更新 README、INDEX、STRUCTURE 和 agent 工作规则；
6. 建立轻量 Git 提交边界，排除原始数据、大型二进制、scratch 和缓存；
7. 核对全局 `Infra/bioinf-data-index` 已覆盖 006 与 013；若新增外部数据再重建索引；
8. 提交并推送 v7 reset。

## 输出

- `docs/plan/v7/plan.md`
- `docs/plan/v7/roadmap.md`
- `docs/v7/V6_ASSET_DISPOSITION.md`
- `docs/v7/SPATIAL_ASSET_INTAKE.md`
- `docs/v7/REPOSITORY_MIGRATION.md`
- `config/v7/data_source_registry.tsv`
- `results/v7/inheritance/asset_manifest.tsv`
- 更新后的根目录入口

## 完成标准

- GitHub 主分支不含大于 100 MB 的文件；
- 所有历史资产仍可定位；
- v7 尚未运行的科学计算明确标记 `NOT_RUN`；
- 根目录所有活跃链接指向 v7；
- 旧结果没有被包装成 v7 结论。

---

# Stage 1｜多模态数据注册与空间资格审计

状态：**NEXT**

## 科学问题

哪些现有数据能支持空间发现、临床锚定、repair 或独立验证？有效样本量到底是多少？

## 输入

- `/006` 的 68-cohort/48-object identity surface 与当前有效 metadata；
- `~/013_spatial/infra/sample-registry/`、structure registry 和数据 briefing；
- 全局 bioinformatics index；
- raw/processed 文件及数据自带 metadata。

## 任务

### 1.1 建立跨仓库 logical-unit registry

字段至少包括：

```text
source_project, dataset_id, patient_id, block_id, section_id, region_id,
sample_id, modality, platform, resolution, coordinate_system, scale,
counts_layer, image_available, segmentation_available, structure_gt_available,
cancer, tissue, treatment, timepoint, response, endpoint,
identity_confidence, metadata_provenance, duplicate_lineage,
n_patients, n_blocks, n_sections, n_auditable_units,
claim_id, primary_role, role_source, permitted_role
```

### 1.2 空间平台资格分层

- Visium/ST/Stereo-seq spot/bin；
- Visium HD；
- Xenium cell/transcript；
- imaging/proteomics；
- h5ad、Space Ranger、matrix+coordinate 等存储格式。

### 1.3 数据角色与 claim 隔离

- HTAN CRC：开发/训练；
- GSE175540 ccRCC、USZ TLS：沿用只读 role ledger，作为 locked external validation；
- ST_CRC_CMS：沿用只读 role ledger，作为 internal/serial-section validation；
- HEST/10x/STOmics：atlas/reference/压力测试，临床 metadata 不足时不进 response；
- GSE238264：HCC 治疗后 spatial support；
- GSE291246：BCC Xenium 支持集，不得计入 HCC；Mendeley：11 个 HCC object，metadata 核验后决定；
- GSE211956：response-linked candidate，需先闭合患者/表达/坐标/结构语义；
- TASK01/TASK02 等非空间队列：临床和纵向 anchor。

### 1.4 重复与泄漏审计

识别镜像、重复 accession、同患者多 section、相邻切片、同一组织块的多格式副本。

## 输出

- `results/v7/registry/multimodal_logical_units.tsv`
- `results/v7/registry/spatial_physical_units.tsv`
- `results/v7/registry/patient_block_section_crosswalk.tsv`
- `results/v7/registry/treatment_response_completeness.tsv`
- `results/v7/registry/duplicate_lineage.tsv`
- `results/v7/registry/data_role_assignment.tsv`
- `results/v7/registry/REGISTRY_REPORT.md`

## 决策点 D1

每个数据集针对每个 claim 获得一个 primary role：`discovery / calibration / internal_validation / external_validation / support / negative_control / reject`。角色可随新 claim 重新登记，但参与过 ontology、候选、阈值或模型选择的数据，对同一 claim 不得再升级为独立验证；secondary role 只能是不承担结论升级的 support。身份不清的数据保持 unknown，不根据文件名或表达模式补标签。

若临床关联空间队列仍不足，空间主线可继续做 response-blind barrier architecture，但“空间结构与 PD-1 failure 相关”必须等待 response-linked spatial 数据或 Stage 4B 通过的 architecture surrogate。二者均无时，临床 response 只能锚定 module/cell state，空间结果保持独立的 Branch A；不得把两条平行证据拼成空间疗效结论。

---

# Stage 2｜共同生物学词汇与跨模态测量

状态：PENDING

## 科学问题

能否在 scRNA、spatial、bulk 和 perturbation 中稳定测量同一组功能对象，而不强迫原始分辨率一致？

## 任务

### 2.1 经验程序基线

- 重放 v6 8 FM 的 gene coverage、carrier、ambient 和 cohort dependence；
- 保留 scoreable modules，不继承 failure/barrier 名称；
- 对随机基因集、matched-expression null 和 leave-cohort-out 稳定性做审计。

### 2.2 机制轴补齐

至少覆盖：

- VEGF/angiogenesis/endothelial；
- CAF/ECM/exclusion；
- hypoxia；
- suppressive myeloid；
- APC/IFN/antigen presentation；
- cytotoxic T/NK；
- progenitor/terminal exhaustion；
- TLS/B；
- neutrophil/TAN/NET；
- WNT/β-catenin exclusion。

机制轴可来自 curated signatures、response-blind empirical programs 或二者的显式组合。不得根据最终 response 效果挑成员后再把同一数据当验证。

### 2.3 Cell-state ontology

建立跨数据的 coarse/mid/fine 映射，保留原始 label、mapping confidence、unknown、mixed 与 low-quality。

### 2.4 跨模态投影

输出 cell/spot/bin/region 和 patient 级分数，同时记录 gene coverage、platform、resolution 和 uncertainty。

同时定义 `architecture_surrogate` 接口：输入只能是非空间队列也可测量的 module、cell state 和临床前可用协变量；目标是后续已发现的 patient-level topology，不使用 response。这里仅定义 schema 与可评估性，不在看见空间结构或 response 前选择最终代理。

## 输出

- `results/v7/ontology/module_dictionary.tsv`
- `results/v7/ontology/module_membership.tsv`
- `results/v7/ontology/cell_state_ontology.tsv`
- `results/v7/ontology/cross_modal_mapping.tsv`
- `results/v7/ontology/measurement_reliability.tsv`
- `results/v7/ontology/architecture_surrogate_contract.yaml`
- `results/v7/ontology/ONTOLOGY_AND_MEASUREMENT_REPORT.md`

## 决策点 D2

输出三类集合：

- 可测量；
- 可用于联合表示但带 uncertainty；
- 仅注释/拒绝。

独立归因资格不作为联合表示的预先删除条件。

---

# Stage 3｜空间基础设施与最小回放

状态：PENDING

## 科学问题

不同空间平台能否在不抹平技术差异的前提下，生成同语义的患者级空间统计？

## 任务

### 3.1 Adapter

至少实现：

- matrix/counts 读取；
- spot/cell/bin identity；
- coordinate system 与单位；
- tissue mask；
- image/segmentation 可选输入；
- gene mapping；
- region/structure GT 隔离；
- patient/block/section 嵌套。

### 3.2 Ground-truth policy

复用 `/013_spatial` 的 `input_policy.tsv`：GT、GT distance、GT-derived mask 和结构定义性变量只能用于允许的验证任务，不能循环进入对应发现输入。

### 3.3 最小平台回放

优先选择 6–10 个高价值 logical units，覆盖至少两类平台、两个癌种和一个已知结构锚点。先证明计数、坐标、尺度、图像和 patient split 可重放，再扩到大库存。

### 3.4 基线空间统计

- kNN/radius graph；
- local composition；
- distance/boundary；
- adjacency/enrichment；
- Moran/variogram；
- spatial permutation；
- 以 patient 为外层的分层 bootstrap；block/section 仅作患者内嵌套层，不增加独立样本数。

## 输出

- `src/v7/spatial_io/`
- `src/v7/spatial_stats/`
- `tests/v7/test_spatial_identity.py`
- `tests/v7/test_spatial_gt_isolation.py`
- `tests/v7/test_patient_block_split.py`
- `results/v7/spatial_foundation/platform_replay_manifest.yaml`
- `results/v7/spatial_foundation/SPATIAL_FOUNDATION_REPORT.md`

## 决策点 D3

只有通过 identity、coordinate、scale、counts semantics 和 GT isolation 的数据进入 Stage 4。平台适配失败的数据保持 support/reject，不用插值或静默替代核心计算。

---

# Stage 4｜Response-blind 空间场与屏障结构发现

状态：PENDING

## 科学问题

是否存在跨患者可重复、超越 abundance 和已知 marker 的空间结构？

## 任务

### 4.1 可解释结构 primitives

建立：

- tumor–immune infiltration/penetration；
- tumor–stroma boundary organization；
- myeloid–tumor suppressive adjacency；
- CAF/ECM/vascular exclusion；
- APC–T-cell support；
- TLS/B-cell neighborhood；
- progenitor/terminal CD8 spatial balance；
- ligand–receptor proximity。

### 4.2 开放潜在场

在不使用 response、结构 GT 和定义性 marker 的发现配置下，比较：

- `K=0`/无潜在场；
- 可解释低秩空间场；
- R-04/MNSF 候选；
- 简单 spatial factor/Gaussian process baseline。

评估优化平台、held-out gene/patient likelihood、空间置换、seed/fold 子空间稳定性和 `K_eff`。当前 R-04 仅完成 restart 稳定性诊断，尚未正式选择真实数据 K；不得继承为 v7 结果。

### 4.3 结构锚定

HTAN、Meylan、USZ、ST_CRC_CMS 的 TLS 和边界 GT 用于独立解释/定位 benchmark，不决定候选场的发现范围。

### 4.4 Topology 增量

比较：

1. QC/platform-only；
2. composition-only；
3. module abundance-only；
4. region-only；
5. spatial smoothing；
6. topology-aware statistical model；
7. latent-field；
8. graph model candidate。

### 4B Response-blind 共享结构对齐与代理验证

在读取任何 response 前：

1. 用 discovery/training 数据将候选结构对齐为 shared-function、context-specific 或 unresolved；
2. 固定一个版本化候选表，记录 carrier、尺度、方向、允许的模态投影和不确定性；
3. 在同时具有空间与转录测量的患者中，拟合 architecture surrogate，并按 patient 做外层留出；
4. 用未参与 ontology、候选、阈值或模型选择的 claim-specific validation 数据检验 surrogate；
5. 分别报告 direct topology、surrogate-predicted topology 和 module/cell-state abundance，禁止把三者混称为空间结构。

## 输出

- `results/v7/spatial_discovery/niche_primitives.tsv`
- `results/v7/spatial_discovery/latent_field_registry.tsv`
- `results/v7/spatial_discovery/barrier_architecture_candidates.tsv`
- `results/v7/spatial_discovery/topology_increment.tsv`
- `results/v7/spatial_discovery/response_blind_shared_architectures.tsv`
- `results/v7/spatial_discovery/architecture_surrogate_validation.tsv`
- `results/v7/spatial_discovery/leave_dataset_out.tsv`
- `results/v7/spatial_discovery/SPATIAL_DISCOVERY_REPORT.md`

## 决策点 D4

- topology 有稳定增量：进入 clinical/context；
- 只有 abundance 稳定：保留 abundance story，停止空间机制升级；
- latent field 不稳定：不阻止可解释结构主线；
- GNN 无增量：删除 GNN；
- 无结构达到跨数据复现：报告 context-specific 或 unresolved，不制造共享结构；
- shared/context 候选在读取 response 前版本化；后续 clinical signal 只能注释或否定，不能反向重定义后再复用同一队列验证；
- surrogate 未在整患者、未触碰数据上证明可转移：非空间队列不得承担 topology claim。

---

# Stage 5｜PD-1 临床锚定

状态：PENDING

## 科学问题

哪些分子和空间结构在可比较临床环境内与 PD-1 response/failure 稳定相关？

## 任务

### 5.1 Patient-level clinical table

统一 treatment arm、endpoint、response、timepoint、sample provenance 和 patient grouping。TASK01 的 biopsy/first-dose/assessment 日期缺口保留 sensitivity 标签；GSE301741 先完成 marker/harmonized bridge 验证。

### 5.2 环境内模型

每个 treatment×endpoint×cancer 分开分析，并明确三条互不替代的测量路径：

- cell fraction；
- molecular programs；
- direct spatial architecture（仅 response-linked spatial 患者）；
- validated architecture surrogate（仅 Stage 4B 通过者，且明确标为 surrogate）；
- clinical/QC；
- abundance vs topology increment。

优先 logistic/ordinal/Cox/mixed-effect 等低容量模型，使用整患者和整队列留出。

### 5.3 跨环境汇总

使用 effect-size meta-analysis 或 hierarchical model，输出 shared、context-modulated、specific、conflicting 或 unresolved。

### 5.4 Negative controls

- cohort-only、treatment-only、endpoint-only；
- response permutation；
- patient/timepoint shuffle；
- random module；
- sampling/QC-only；
- leave-cohort/cancer-out。

## 输出

- `results/v7/clinical_anchor/patient_level_table.*`
- `results/v7/clinical_anchor/environment_effects.tsv`
- `results/v7/clinical_anchor/confounding_baselines.tsv`
- `results/v7/clinical_anchor/responder_compatible_ecology.tsv`
- `results/v7/clinical_anchor/CLINICAL_ANCHOR_REPORT.md`

## 决策点 D5

若 cohort-only 仍显著强于生物模型，禁止 shared PD-1 claim；可以保留环境内描述或回到 metadata/measurement 修复。不存在 shared direction 也是合法结果。

若既无 response-linked spatial 数据、也无通过 Stage 4B 的 surrogate，则 spatial-architecture response association 标记 `NOT_IDENTIFIABLE`；只能形成 module/cell-state 临床锚和独立空间图谱，不能声称该空间结构与 PD-1 failure 相关。

---

# Stage 6｜共享结构的临床注释、器官 context 与 HCC 深入分析

状态：PENDING

## 科学问题

Stage 4B 已响应盲定义的共享功能结构如何与临床结局相关？HCC 如何通过肝脏耐受、髓系、血管和基质改变其实现？

## 任务

1. 对 Stage 4B 已经 response-blind 版本化的共同功能结构做临床与器官注释，不用 response 重新聚类或选结构；
2. carrier、尺度、区域和方向的 context decomposition；
3. leave-one-cancer/cohort/platform；
4. HCC scRNA、spatial 与 bulk 的独立复现；
5. etiology 只作覆盖充分时的 HCC sensitivity；
6. 检查 pure abundance、ambient、取材区域和组织纯度解释。

## 输出

- `results/v7/context/shared_architecture_clinical_annotations.tsv`
- `results/v7/context/context_implementations.tsv`
- `results/v7/context/hcc_residual_architectures.tsv`
- `results/v7/context/context_heterogeneity.tsv`
- `results/v7/context/SHARED_AND_HCC_CONTEXT_REPORT.md`

## 决策点 D6

每个对象标为：shared、shared_with_context_modulation、HCC/context-specific、conflicting 或 unresolved。HCC 结果不自动外推到泛癌。

---

# Stage 7｜PD1+X 纵向重排与 repair 判断

状态：PENDING

## 科学问题

PD1+X 是否与 residual barrier 的削弱和 responder-compatible ecology 的建立相容？

## 输入优先级

- TASK01：PD-1 单药 reference trajectory；
- TASK02：PD1+lenvatinib paired/longitudinal；
- LAMBRECHT_HCC：atezolizumab+bevacizumab 外部 antiangiogenic-class 支持；
- GSE238264：PD1+cabozantinib post-treatment spatial response support；
- 其他经 Stage 1 核验的 mono/combo paired assets。

## 任务

### 7.1 治疗保留的预处理

不能使用已回归掉 timepoint/treatment 的 residualized matrix。只去除技术 nuisance，保留 arm、time、arm×time 和 patient pairing；composition-inclusive 与 composition-adjusted 并列。

### 7.2 患者内位移

计算 T1−T0、T2−T0 或 pre→post 的 module 与 cell-state displacement。空间部分分开处理：同患者纵向空间数据计算 `observed_spatial_rewiring`；Stage 4B 验证通过的代理计算 `surrogate_inferred_architecture_rewiring`。GSE238264 等 post-only 空间队列只能用于横断面结构兼容性。

### 7.2A 空间 repair 可识别性硬边界

- 同患者 pre/on/post spatial 存在：允许直接检验空间结构重排；
- 没有纵向 spatial，但 architecture surrogate 在独立患者中通过：只允许代理推断，并保留 surrogate 标签；
- 两者均无：`architecture_rewiring=NOT_IDENTIFIABLE`，Branch C 阻塞；模块/细胞状态纵向 repair 可继续，但不能写成空间结构重排。

### 7.3 Repair compatibility

比较：

- responder centroid/direction；
- regularized linear displacement；
- mixed-effect trajectory；
- entropic OT；
- 必要时的小型非线性 transport。

### 7.4 对照

- pairing permutation；
- timepoint shuffle；
- composition-only；
- cohort/QC-only；
- PD-1 mono vs combo 的非因果敏感性；
- 机制不匹配/临床失败组合。

## 输出

- `results/v7/repair/paired_displacements.tsv`
- `results/v7/repair/observed_spatial_rewiring.tsv`
- `results/v7/repair/surrogate_inferred_architecture_rewiring.tsv`
- `results/v7/repair/repair_compatibility.tsv`
- `results/v7/repair/negative_combination_controls.tsv`
- `results/v7/repair/REPAIR_REPORT.md`

## 决策点 D7

- 有 paired/外部一致证据：升级 repair-compatible hypothesis；
- 只有组合队列横断面关联：只称 combination-associated state；
- mono/combo 不可交换：不得声称 X 增量效应；
- 位移方向冲突或由 QC/composition 解释：候选降级或否定。
- 无纵向 spatial 且 surrogate 未通过：空间 repair 保持 `NOT_IDENTIFIABLE`，不以 post-only 空间支持替代。

---

# Stage 8｜真实扰动、靶点可达性与外部患者证据

状态：PENDING

## 科学问题

候选 X 是否能在合理 carrier 中沿预期方向移动 barrier architecture 的分子组成？患者级外部数据是否支持？

## 任务

### 8.1 真实 perturbation mapping

使用 GSE133344、GSE193736、GSE306429、GSE90063、XAtlas、L1000 等真实测量，重算 perturbation→module/cell-state direction。v6.1 hash-derived score 禁止使用。

### 8.2 Target/carrier/spatial reachability

检查 target/pathway 是否在 tumor、myeloid、T cell、CAF、endothelial 等合理 carrier 和空间区域中表达；缺覆盖不自动解释为机制为负。

### 8.3 External bulk/scRNA

将固定 module/architecture surrogate 映射到 IMbrave150、GSE274975、IMvigor210 等合格队列，先核验 treatment、endpoint 和 mapping 后再分析。

### 8.4 Negative controls

- random perturbation；
- target-label permutation；
- unrelated cell context；
- failed combination；
- gene coverage matched null。

## 输出

- `results/v7/perturbation/perturbation_direction.tsv`
- `results/v7/perturbation/target_carrier_reachability.tsv`
- `results/v7/external_validation/patient_effects.tsv`
- `results/v7/external_validation/negative_controls.tsv`
- 对应报告

## 决策点 D8

spatial、clinical、perturbation 和 external 方向冲突时，不做加权平均消除冲突；Evidence Card 明确列出并降低 claim。

---

# Stage 9｜X-class mechanism space 与 class-held-out 外推

状态：PENDING；需 Stage 7–8 提供至少三类 X-class

## 科学问题

barrier–repair mapping 能否推广到未见 X-class，而不是记住药名和训练组合？

## 任务

1. 版本化 X-class ontology：target、pathway、direction、carrier、spatial reach、known constraints；
2. 纳入成功、失败和机制不匹配组合；
3. leave-one-X-class / leave-one-combination；
4. 比较 pathway overlap、nearest prototype、additive baseline；
5. semantic leakage audit；
6. calibrated abstention。

## 输出

- `results/v7/repair/x_class_ontology.tsv`
- `results/v7/repair/class_held_out_results.tsv`
- `results/v7/repair/abstention_calibration.tsv`
- `results/v7/repair/X_CLASS_GENERALIZATION_REPORT.md`

## 决策点 D9

不足三类 X-class 或没有可信负例时，节点保持 `NOT_IDENTIFIABLE`；不把单一 antiangiogenic 案例包装成 zero-shot。

---

# Stage 10｜Evidence Card、独立复核与论文闭合

状态：PENDING

## 科学问题

哪些结论已经达到可发表最低充分证据？哪些仍是线索、冲突或未回答？

## 任务

1. 每个 claim 建立 Evidence Card；
2. 整患者/整队列复现；
3. data leakage、pseudo-replication、causal overclaim hostile review；
4. 失败和负结果单独成表；
5. 形成数据、代码、环境和图表复现包；
6. 更新 Data Availability、数据索引和资产清单；
7. 根据实际证据选择论文主线。

## 可形成的论文分支

### Branch A：泛癌空间屏障结构

当空间 topology、潜在场和 context transfer 成立，但临床 repair 较弱时使用。

### Branch B：PD-1 failure 与 HCC 空间实现

当 clinical anchor、shared/context 和 HCC replication 成立时使用。

### Branch C：PD1+X repair-compatible rewiring

当 paired trajectory、空间结构和 perturbation/external evidence 闭合时使用。空间结构必须来自同患者纵向空间观测，或来自通过独立整患者验证并明确标注的 architecture surrogate；否则 Branch C 阻塞。

### Branch D：X-class 外推框架

只有多 X-class class-held-out 证据成立时使用。

分支不是降级路线，而是由数据决定的结论边界。可合并的前提是各自证据均独立成立。

## 输出

- `results/v7/evidence_cards/*.yaml`
- `results/v7/evidence_cards/claim_matrix.tsv`
- `results/v7/evidence_cards/conflict_register.tsv`
- `results/v7/figures/`
- `results/v7/FINAL_SCIENTIFIC_REPORT.md`
- 可复现代码与轻量 manifest

---

## 4. 最小验证矩阵

| 分析 | 必须比较的基线 | 独立单位 | 最低验证 |
|---|---|---|---|
| module/cell-state measurement | random/matched genes、cohort/QC | patient/cohort | leave-cohort-out |
| spatial structure | composition、abundance、region、smoothing | patient；block/section 为 nested repeat | leave-patient/dataset-out + spatial permutation |
| latent field | K=0、简单 factor/GP | patient/fold/seed | optimization + held-out + subspace stability |
| response association | cohort/treatment/endpoint/QC-only | patient | nested/grouped CV + external cohort |
| HCC context | composition/platform/region | cohort/patient | independent HCC dataset |
| paired rewiring | unpaired、time/pair permutation | patient-pair | independent combination cohort；空间结论另需 longitudinal spatial 或 validated surrogate |
| perturbation | random perturb/target labels | perturb-condition | replicate/context transfer |
| X-class | pathway overlap/prototype | X-class/cohort | class-held-out + abstention |

---

## 5. 高风险点与提前处置

### 风险一：空间库存大但临床标签少

处置：让无标签数据承担发现、结构复现和平台泛化；疗效仍由少量可比患者队列锚定。

### 风险二：患者身份和重复来源不完整

处置：Stage 1 建立 duplicate lineage；unknown 不进入独立验证计数。

### 风险三：跨平台分辨率被误认为生物学

处置：平台专属 adapter、患者级统计、platform-held-out 和 scale sensitivity；不强制原始 embedding 一致。

### 风险四：空间结构被 cell fraction 解释

处置：abundance-only 是正式强基线；无增量则否定空间机制主张。

### 风险五：模型过度复杂

处置：透明统计→latent field→graph 的逐层必要性比较；复杂模型失败不阻止生物学主线。

### 风险六：治疗与队列绑定

处置：环境内估计、患者内轨迹和整队列外部验证；不用 mono/combo 绝对均值推导因果 X effect。

### 风险七：response-guided discovery 泄漏

处置：response-blind foundation、嵌套选择、未触碰 cohort、完整 permutation。

### 风险八：013 的 R-04 长任务与 v7 冲突

处置：不修改、不停止、不复用未完成科学结论；v7 仅通过只读 registry/manifest 接口继承其稳定资产。

---

## 6. 当前立即执行顺序

本次 v7 reset 完成后，下一轮只应启动以下三项：

1. **Stage 1 cross-repository registry**：先知道真实患者、block、section、治疗和 response 面；
2. **Stage 2 module/cell-state vocabulary audit**：补足 vascular/CAF/myeloid/exhaustion 等空间 repair 轴；
3. **Stage 3 minimal spatial replay**：用少数高价值数据证明多平台接口、GT 隔离和患者级统计成立。

在这三项完成前，不启动全量空间训练、PD1+X 复杂模型、GNN 或 zero-shot X-class。

---

## 7. 当前状态声明

截至 2026-08-31：

- v7 科学方案与路线图已形成；
- v6/013 资产已完成初步只读继承审计；
- 仓库正在重组并准备 GitHub 提交；
- v7 registry、ontology、spatial foundation、clinical anchor、repair 和 X-class 计算均为 **NOT_RUN**；
- 当前没有 v7 新生物学结论。
