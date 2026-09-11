# Stage 3 决策记录

## 必需任务清单

| requirement_id | original_requirement | mandatory | execution_status | result | evidence_locator | deviation_id | downstream_impact |
|---|---|---|---|---|---|---|---|
| S3-001 | Adapter（counts/identity/坐标/tissue mask/gene mapping/GT 隔离/嵌套层级） | 是 | COMPLETE | 11/11 捕获通过 | `adapter_qc.tsv`、`count_conservation.tsv`、`D3_GATE.json` | — | Stage 4 输入合同 |
| S3-002 | Ground-truth policy（GT/GT distance/GT-mask/定义性变量隔离） | 是 | COMPLETE | GT 隔离通过；锁定的验证 GT 未解析 | `gt_isolation_audit.tsv`、`sealed_validation_manifest.tsv` | — | 验证 GT 不得循环进发现 |
| S3-003 | 最小平台回放（6–10 高价值 units；≥2 平台/2 癌种/1 结构锚点） | 是 | COMPLETE | 8 pilot/11 捕获/7 患者；Visium＋Xenium；HTAN/GSE238/GSE291 | `pilot_manifest.tsv`、`platform_replay_manifest.yaml` | — | 回放语义闭合 |
| S3-004 | 基线空间统计（kNN/组成/距离/富集/Moran/variogram/置换/患者外层 bootstrap） | 是 | COMPLETE | section-local 图＋99 次置换＋1,000 次患者 bootstrap | `graph_statistics.parquet`、`spatial_null_statistics.parquet`、`patient_bootstrap_ci.tsv` | — | Stage 4 统计底座 |
| S3-005 | Stage 4 扩展资产连同一合同；新增/改变路径才回放；kNN≠真实接触 | 是 | COMPLETE | Stage 4 60 units 复用同一 adapter/身份/坐标合同（见 Stage 4 选择表 provenance）；kNN 语义显式限定 | `EVIDENCE_AND_CONFLICTS.md:扩展资产核对` | S3-D01 | 扩展资产不得伪造接触解释 |

结果规则：技术回放 COMPLETE；不判定屏障/疗效/潜在场。
范围变更：无。
reviewer/date：待第三方实际审查时填写（见文末模板）。

## 偏移表

| deviation_id | original | revised | reason | benefit_evidence | claim_loss | exposure_change | decision | decision_by | date |
|---|---|---|---|---|---|---|---|---|---|
| S3-D01 | 跨平台统一构图/同分布 embedding | 按平台独立构图＋统计层对齐；panel 缺失不补零 | 分辨率/坐标/误差不同；缺失补零伪造生物学 | 预期：语义保真；证据：`feature_coverage.tsv` 缺失保留＋1 共同可测特征非退化统计 | 无；跨平台 claim 限统计层 | 无 | ACCEPTED | repair-audit | 2026-09-11 |

## 第三方审查记录模板

| 项目 | 内容 |
|---|---|
| scope_version / run_id | plan_v7.0 / STAGE3_RUN_MANIFEST.yaml |
| reviewer / date / independence | 待填写（须未承担本项实现与结果选择） |
| inputs and outputs inspected | OUTPUT_MANIFEST.yaml 路径＋hash；抽查 adapter/counts/坐标/GT 隔离/split 行 |
| computation reproduced | 命令见 OUTPUT_MANIFEST；第三方优先复算小判定表，不默认重跑全矩阵 |
| original requirements unresolved | 无（限制见 EVIDENCE，不属未完成） |
| deviations reviewed | S3-D01 待 review |
| claim decisions | 技术回放（非生物学 claim） |
| stage closure / downstream permission | D3 PASS_WITH_LIMITATIONS；下游见 NEXT |
| remaining actions | 无；Stage 4 扩展资产沿用本合同 |
