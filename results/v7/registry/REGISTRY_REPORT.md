# Stage 1 registry report

版本：`v7-stage1-registry-1.0`  
范围：跨 `/006` 与只读 `/home/huyudi/013_spatial` 的元数据、身份、角色、重复和治疗/疗效字段审计。  
运行边界：未打开表达矩阵、图像、h5ad 内容或大型二进制；未运行 v7 科学模型。

## 生成与输入

可复现命令：

```bash
python scripts/v7/registry/build_stage1_registry.py
```

输入文件哈希（本次仅读取元数据）：

- `assets`: `/home/huyudi/013_spatial/infra/sample-registry/source_assets.tsv` (29342 bytes, SHA256 `a423137d8967568aa830e91b46f39d0ac280f502fa8bfdd66d3e7fd9e0f8288a`)
- `cohort`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/metadata/cohort_registry.effective_v1.csv` (18405 bytes, SHA256 `671f2c5948c03b58cefffb6301eb89c65eb96e2328e1eb82747deccb54cb78ad`)
- `config`: `/home/huyudi/006/config/v7/data_source_registry.tsv` (5527 bytes, SHA256 `fce4577d8da22955b4cdcfd46af9e929b6909ded12061608248a50b3d4f51b0b`)
- `context_sidecar`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/analysis_context_sidecar.csv` (349470 bytes, SHA256 `84c68860d565d4b033252d303eb333330306e671ec48b45615be274f9c8b8174`)
- `cross_sections`: `/home/huyudi/013_spatial/infra/structure-registry/cross_section_links.tsv` (472 bytes, SHA256 `c2dcae79d4402e0f07ecefac6f90b552f3655623d480a4e1be0ea26528124ddb`)
- `duplicates`: `/home/huyudi/013_spatial/infra/sample-registry/duplicate_groups.tsv` (40222 bytes, SHA256 `2ffd0523bcf7b9a1a1340513ea0d5a9eff2916ea9596580bafb8605fc061438f`)
- `feature`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/metadata/dataset_role_and_feature_eligibility.effective_v1.csv` (16630 bytes, SHA256 `0569f6eda95ee1497909904a0792a00153d80b2145b2c868e0b3deee511a1579`)
- `gse238_manifest`: `/home/huyudi/006/results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.GSE238264.tsv` (1817 bytes, SHA256 `29aadc370246867c45aebddb25de3960a87dfb2a263d4d237277a5521b077296`)
- `gse238_qc`: `/home/huyudi/006/results/v6_2/phase2_5_data_onboarding/02_conversion_qc/spatial_sample_manifest.GSE238264.tsv` (1817 bytes, SHA256 `29aadc370246867c45aebddb25de3960a87dfb2a263d4d237277a5521b077296`)
- `gse238_spot_qc`: `/home/huyudi/006/results/v6_2/phase2_5_data_onboarding/02_conversion_qc/per_sample_spot_stats.GSE238264.tsv` (1779 bytes, SHA256 `1a17f1c70eab7fac76d1e9247c5e019a5d87bc4ca3e1f4f22be15775464a2f04`)
- `gse291_manifest`: `/home/huyudi/006/results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.GSE291246.tsv` (12039 bytes, SHA256 `a6e5d44d74550a775fdefe26404983bf8afa3ec5f860be49883c9634a3f7ee70`)
- `gse291_qc`: `/home/huyudi/006/results/v6_2/phase2_5_data_onboarding/02_conversion_qc/per_sample_xenium_stats.GSE291246.latest.tsv` (6263 bytes, SHA256 `ffd940332302731d014201af04d7defebbf4e3be216dc608a51ae5b2f90e52dc`)
- `gt_audit`: `/home/huyudi/013_spatial/infra/structure-registry/gt_source_audit.tsv` (83359 bytes, SHA256 `54c6dc42db9af60f7f24e5007b9fc2eef8c2de63c08ee8ebbfea04b872b2720d`)
- `identity_evidence`: `/home/huyudi/013_spatial/infra/sample-registry/identity_evidence.tsv` (2028680 bytes, SHA256 `e5c46a05f2bf4edfb403a34b7aa2aa86a71e36572a109653eb661c5f93587a03`)
- `input_policy`: `/home/huyudi/013_spatial/infra/structure-registry/input_policy.tsv` (1444 bytes, SHA256 `520c0aa01430f78dfb1569fbcd6e9797f8e9ee5dc79a18280f0e4753298efa79`)
- `object_audit`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/object_metadata_audit.csv` (4881 bytes, SHA256 `c1af40b16fc5c252496ea57f96bc090f100bfc6fd34f8af9c1a1303588aaa038`)
- `ontology`: `/home/huyudi/013_spatial/infra/structure-registry/structure_ontology.tsv` (1095 bytes, SHA256 `117a769374411ca6cc6b35a56dec3b88b13f615534733f38e97a8ecb81f02c6e`)
- `outer_splits`: `/home/huyudi/013_spatial/infra/structure-registry/outer_splits.tsv` (49405 bytes, SHA256 `ad3807002d39b4ba3f48b7fa39ec41fe17b192621ea5a7a81774b4a5fbd3e26e`)
- `patient`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/metadata/patient_metadata_master.effective_v1.csv` (362492 bytes, SHA256 `a55f5d43e45191c8d36d55bbd4ff2a79d232ecc2834c5e3d91a4759923eb6e2c`)
- `physical`: `/home/huyudi/013_spatial/infra/sample-registry/physical_units.tsv` (218079 bytes, SHA256 `367d94e50e43c32722d1aa2e343963c6cbeb5b7f6f6f8d189d3b7de22092d219`)
- `r04_input_manifest`: `/home/huyudi/013_spatial/infra/r04/input_manifest.json` (52217 bytes, SHA256 `2ab11adc3b15deb68a14a20666cf7eb42df939db05479ff7731603b217bb1e4e`)
- `r04_restart_panel`: `/home/huyudi/013_spatial/infra/r04/restart_stability_panel_20260831.json` (56059 bytes, SHA256 `2a711cc254ba87db9f646cf468bfd1a6e3477d38519777fdf7fa2059ef40c276`)
- `raw_task`: `/home/huyudi/004/TASK/meta_task_raw.csv` (12412 bytes, SHA256 `f2068fa971e50f764314fd5b260f43457d4aacd6faec7c2035520c5417db6f1c`)
- `raw_tcr`: `/home/huyudi/004/Component/meta_tcr.csv` (1734 bytes, SHA256 `274bfe84701c2300d5b99e2d2437dad72bfdbb0c378b678170a8297a765ba037`)
- `replay`: `/home/huyudi/013_spatial/infra/structure-registry/h5ad_replay_index.tsv` (19512 bytes, SHA256 `b2855c76fe426b95cd2df0488048f2022aa3c88f8fb542d2faeb94eefae88242`)
- `response`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/metadata/response_label_environment.effective_v1.csv` (18495 bytes, SHA256 `f5bb25dc31b067d238ba3f6cb7d1fd43c83e0653e926e9dbfcafd4d4a8571869`)
- `role_freeze`: `/home/huyudi/013_spatial/infra/sample-registry/role_freeze.tsv` (2800 bytes, SHA256 `f02c9f7a862bda345b1bb913a3d495b118ea2f0dcd41df4d38ebebe8429ff96e`)
- `sample`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/metadata/sample_metadata_master.effective_v1.csv` (3810950 bytes, SHA256 `4b91eca9834640bfa422805e5ea1905991288a97a2aa55e3250812a970d7158a`)
- `sample_sidecar`: `/home/huyudi/006/results/v6_2/data_interface_repair_v1/sample_identity_sidecar.csv` (2443349 bytes, SHA256 `ba06a636378fec1bf80f59cc35ca187a657e057a4e711cb2390ec91ce21e56cd`)
- `source_inventory`: `/home/huyudi/013_spatial/infra/sample-registry/source_metadata_inventory.tsv` (8513 bytes, SHA256 `12cfd2ed8cc0035b483479550fda724705702a2bae98dc3b8ba95a51b8e3b0b4`)
- `structures`: `/home/huyudi/013_spatial/infra/structure-registry/structure_instances.tsv` (907809 bytes, SHA256 `166e042fd926caf2c4ef81efdd55fc1f5909d8868c91aa8ecd5f0f2be536e8aa`)

全局 `Infra/bioinf-data-index` 已覆盖 `/006` 与 `/013_spatial`；本阶段未下载新增外部数据，因此不改写全局索引。

## 账本规模

| 账本 | 数量 | 解释 |
|---|---:|---|
| `/006` effective cohort rows | 76 | 队列级元数据 |
| `/006` effective patient rows | 1670 | 患者级候选，不等于可监督患者 |
| `/006` effective sample rows | 10166 | 样本/时间点级记录 |
| `/013` physical rows | 1002 | 含排除记录的完整审计面 |
| `/013` resolved included candidates | 167 | 可进入后续空间候选池的物理单位 |
| `/013` confirmatory structure instances | 889 | 仅已登记结构 GT，不代表疗效证据 |
| Stage 1 logical-unit rows | 11191 | 样本、物理单位和 dataset-level 账本 |
| Stage 1 spatial physical rows | 1044 | `/013` + `/006` 两个空间 manifest |
| Stage 1 patient/block/section rows | 1044 | patient 是外层独立单位 |

按 source 的 logical-unit 行数：

| source | rows |
|---|---:|
| `006_GSE238264` | 8 |
| `006_GSE286827` | 58 |
| `006_GSE291246` | 36 |
| `006_GSE301741` | 28 |
| `006_LAMBRECHT_HCC` | 113 |
| `006_Mendeley_HCC` | 13 |
| `006_TASK01` | 73 |
| `006_TASK02` | 67 |
| `006_bulk` | 478 |
| `006_perturbation` | 1 |
| `006_sc_identity` | 9302 |
| `013_10x_cancer` | 3 |
| `013_ATLAS_TABLE_S2` | 341 |
| `013_GSE175540` | 25 |
| `013_GSE211956` | 9 |
| `013_GSE226997` | 5 |
| `013_GSE274103` | 6 |
| `013_GSE274557` | 56 |
| `013_HEST` | 493 |
| `013_HTAN_CRC` | 51 |
| `013_STOmicsDB` | 1 |
| `013_ST_CRC_CMS` | 15 |
| `013_TLS_USZ` | 9 |

## Stage 1 关键发现与资格结论

1. `/006` 的有效主表优先于旧摘要：TASK01 为 27 患者/72 样本行，TASK02 为 23/66；Lambrecht HCC 为 44 患者/112 样本行（pre 65、post 47）。`config/v7/data_source_registry.tsv` 中的 Lambrecht “25 pre/22 paired”是旧口径，已作为冲突保留，不能直接沿用。
2. GSE238264 是 7 名 HCC 患者的治疗后空间队列（nivo+cabo；manifest response 计数：non_responder=3, responder=4），可作治疗后空间关联和 support；无基线，不能分离 PD-1 与 cabozantinib，也不能建立纵向空间重排因果链。
3. GSE291246 是 BCC 而非 HCC：35 个 Xenium section，17 个 h5ad/QC 可重放记录；有效主表的时间点口径为 `post_PD1=27, treatment_naive=8`，response 未提供，故只作跨癌种空间 support。QC 表中患者列可能为空，本账本以 metadata-freeze manifest 的患者映射为准。
4. `/013_spatial` 的身份控制面纳入 1,002 条 physical rows；其中 167 条 `RESOLVED_INCLUDED_CANDIDATE` 才是后续空间候选池。block/section/spot/cell 是患者内嵌套观测，不能增加独立 n；同患者、同块、serial section 通过 leakage group 约束拆分。
5. GSE211956 的样本标题 response token 仍处于语义待核验状态；Mendeley 的有效主表目前是 6 名患者/12 个样本行，且本地 h5ad 被标为 scRNA（旧摘要中的 11 个 object 口径不作为 v7 计数）。两者均未升级为 response-bearing 或独立验证。
6. `/013_spatial/infra/r04` 只作为输入定位与重启稳定性诊断的 provenance；本阶段不继承 latent-field 的 K 或结构结论，且最新面板仍不是正式 K 选择。
7. 物理单位的来源反向覆盖已闭合：HEST、ATLAS_TABLE_S2、TENX 和 3 个额外 GEO 空间队列分别进入受控 source；所有 spatial `source_id` 都能在 config 和 role ledger 中找到，未再使用无角色的 `013_other_spatial`。

## 角色账本（D1）

角色行数（8 个 claim × source）：

| primary role | rows |
|---|---:|
| `BCC_cross_cancer_spatial_support` | 8 |
| `HCC_PD1_plus_X_longitudinal_rewiring` | 6 |
| `HCC_PD1_reference_and_longitudinal_anchor` | 8 |
| `HCC_external_combination_trajectory_validation` | 6 |
| `HCC_post_treatment_spatial_response_support` | 6 |
| `HCC_spatial_structure_support` | 8 |
| `HNSCC_anchor_candidate` | 8 |
| `adapter_and_platform_stress_test` | 8 |
| `atlas_reference_and_identity_audit` | 8 |
| `atlas_reference_and_platform_stress_test` | 8 |
| `candidate_pending_response_semantics` | 8 |
| `cross_platform_candidate_discovery` | 8 |
| `external_direction_support` | 8 |
| `external_patient_level_direction_and_clinical_scale` | 8 |
| `external_validation` | 10 |
| `external_validation_locked` | 6 |
| `internal_validation` | 8 |
| `intervention_to_program_direction` | 8 |
| `pan_cancer_response_blind_program_and_carrier_discovery` | 8 |
| `spatial_context_support` | 16 |
| `spatial_discovery_support` | 8 |
| `support_context_or_repair` | 6 |
| `training` | 8 |

角色来源分开记录：013 旧 `role_freeze.tsv` 只作为 provenance 和历史泄漏边界；同一 v7 claim 尚未运行选择/阈值/模型，因此 `selection_exposure` 统一写为 `not_run_v7_stage1_metadata_audit`。locked external 数据不进入候选或调参。

## 治疗/疗效完整性边界

`treatment_response_completeness.tsv` 以 dataset×patient 为行，保留 raw/normalized 值、端点、时间点和冲突标志。只有 treatment、timepoint、endpoint、response 都存在且无冲突的行才标记 `complete_for_response=yes`；这只是字段完整，不自动获得临床主张资格。

- TASK01：PD-1 anchor 环境；仅在自身治疗/终点环境内使用，不能与 PD1+X 直接 pooled。
- TASK02：PD1+lenvatinib 环境；用于 paired/context repair，不作 PD-1 单药交换对照。
- Lambrecht HCC：治疗组合混杂，作为纵向 support，先解决有效主表与旧摘要口径冲突。
- GSE238264：post-only PD1+cabozantinib spatial support。
- GSE291246：post_PD1 BCC，response missing。

## 阻隔与下一阶段

- 直接纵向空间 PD1+X rewiring：`NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`；当前没有同患者空间 pre/on/post 链，post-only 不能替代。
- response-linked spatial：GSE211956 保持 `response_semantics_pending`，待治疗药物、时间点、endpoint 和表达/坐标/结构链接闭合。
- Mendeley HCC：patient/section/response/image-scale provenance 未闭合。
- Stage 2 可启动：先用 response-blind 程序和 cell-state 词汇建立跨模态测量接口；不得把 Stage 1 账本当成生物学结果。

## 验收状态

- 主键和确定性排序：由生成器统一处理；重复运行应 byte-identical。
- 空 patient/block/section：保持空值并显式标记，不把 unknown 合并成共同 ID。
- `/013_spatial`：本阶段只读，未写入。
- 科学结论：`NOT_RUN`；本报告只证明数据资格、有效样本量和边界。
