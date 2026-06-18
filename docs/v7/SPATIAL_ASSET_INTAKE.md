# v7 空间资产接入说明

日期：2026-08-31
来源：`/home/huyudi/013_spatial` 只读审计与 `/006` 现有空间资产

## 判断

`~/013_spatial` 约有 1.9 TB 空间资产，足以把空间分析提升为主线，但数据量不等于独立证据量。多数大规模 HEST、10x、STOmics 资产缺治疗、时间点或 response；它们适合 response-blind atlas、跨平台压力测试和候选发现，不适合直接训练疗效模型。

## 已审计的高价值层

| 资产 | 当前可审计规模 | v7 优先角色 | 主要边界 |
|---|---:|---|---|
| HTAN CRC | R-04 使用 30 patient / 31 block / 47 section | 开发、空间基础、内部复现 | 不是系统治疗 response cohort |
| GSE175540 ccRCC | 24 sample，18 个可审计 TLS source | TLS/跨癌外部验证 | 按患者切分，不能按 TLS component 计独立样本 |
| ST_CRC_CMS | 7 patient、14 section、12 个可审计 source | 边界和跨切片稳定性 | 相邻切片不是独立患者；2 张切片无 boundary audit row |
| TLS_VISIUM_USZ | 肾癌 3 + 肺癌 5 | 跨癌/跨平台 TLS 验证 | 小患者数 |
| Heiser/structure registry | TLS 与 tumor-stroma boundary GT | 结构解释和 benchmark | GT 与派生距离不得进入对应发现输入 |
| HEST-1k 子集 | 492 sample，跨 Visium/ST/Xenium/HD | atlas/reference/平台泛化 | treatment metadata 几乎缺失 |
| 10x 近年癌症空间集 | 73 dataset | adapter 与平台压力测试 | 需要去重、患者和许可审计 |
| STOmicsDB 资产 | 78 dataset | Stereo-seq/跨平台候选 | 临床字段不完整 |
| GSE211956 | 可见 response token 的候选 | 临床关联候选 | 表达、坐标、patient crosswalk、结构语义尚未闭合 |
| `/006` GSE238264 | HCC 7 patient，post PD1+cabozantinib，4R/3NR | HCC 治疗后空间支持 | 不能作 baseline prediction 或拆分 PD1 与 X 效应 |
| `/006` GSE291246 | BCC Xenium；35 section，其中 17 个 h5ad PASS | 跨癌种空间结构支持 | 18 张缺 transcripts；无 response 标签，不能作 HCC 证据 |
| `/006` Mendeley | 11 个 HCC spatial RDS object | HCC 组织结构候选 | patient/section、response、图像与 scale contract 均待闭合 |

## 可复用接口

- `~/013_spatial/infra/sample-registry/physical_units.tsv`
- `~/013_spatial/infra/structure-registry/structure_instances.tsv`
- `~/013_spatial/infra/structure-registry/input_policy.tsv`
- R-04 的 patient-grouped split、counts+coordinate loader、checkpoint 和审计框架
- R01 注册/去重与 R02 ground-truth 注册代码

R-04 最新诊断状态为 `RESTART_STABILITY_DIAGNOSTIC_COMPLETE_NOT_FORMAL_K_SELECTION`。这说明多 restart 稳定性诊断已完成，但真实数据的正式 K 选择和科学场命名仍未完成；v7 可复用接口与诊断框架，不能继承已成立的 latent field 或有效 K。

## v7 首个空间合同

每个 logical unit 至少登记：

```text
patient_id, block_id, section_id, region_id, modality, platform,
resolution, coordinate_system, coordinate_unit, scale_factor,
counts_layer, normalized_layer, image_available, segmentation_available,
ground_truth_available, treatment, timepoint, response, endpoint,
metadata_source, identity_confidence, duplicate_lineage,
n_patients, n_blocks, n_sections, n_auditable_units,
claim_id, primary_role, role_source, permitted_role
```

## Claim 角色边界

GSE175540 与 TLS_VISIUM_USZ 保持未触碰的 external-validation 角色，ST_CRC_CMS 保持 internal/serial-section validation 角色。它们不得参与同一 claim 的 ontology、候选、阈值或模型选择；一旦使用，只能降为 development/support，不能在看过结果后重新命名为独立验证。该边界只防止选择后验证，不限制探索性假设和模型迭代。

## 当前 repair 数据缺口

当前已核对资产中没有闭合的同患者 pre/on/post spatial PD1+X 链。TASK02 与 LAMBRECHT_HCC 提供纵向 scRNA，GSE238264 只提供 post-treatment spatial。因此，目前可直接做的是分子/细胞状态纵向位移与横断面空间兼容性；直接空间结构重排为 `NOT_IDENTIFIABLE`。只有 Stage 1 找到真实纵向空间资产，或 Stage 4B 建立并独立验证 architecture surrogate 后，才能升级空间 repair 分支。

## 主线职责

1. 定义 barrier 的组织形式，而不是只画 module heatmap；
2. 判断 topology 是否提供 abundance 之外的信息；
3. 定位分子程序的 carrier、边界、邻接与空间尺度；
4. 识别共享功能结构与器官特异实现；
5. 支撑 HCC 中 vascular/myeloid/CAF/exclusion context；
6. paired 数据存在时检验 barrier rewiring；
7. 约束 X 的 target reachability 与 repair direction；
8. 在结果冲突时提供反证，而不是只作正向验证。

## 禁止事项

- 不把 TB、文件、spot、cell、section 或 GT component 当患者样本量；
- 不跨患者建边；不同切片也分别构图，不因属于同一患者而直接连边；
- 不从文件名静默补 patient/treatment/response；
- 不把横断面空间邻近解释为治疗因果；
- 不全量启动 1.9 TB 训练后再补接口；
- 不修改或干扰 `~/013_spatial` 当前运行任务。
