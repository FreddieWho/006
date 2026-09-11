# Stage 1 决策记录

## 必需任务清单

| requirement_id | original_requirement | mandatory | execution_status | result | evidence_locator | deviation_id | downstream_impact |
|---|---|---|---|---|---|---|---|
| S1-001 | 跨仓库 logical-unit registry（含身份/模态/治疗/重复字段） | 是 | COMPLETE | 11,191 logical rows，主键闭合 | `multimodal_logical_units.tsv`、`OUTPUT_MANIFEST.yaml` | — | Stage 4–9 候选来源的身份依据 |
| S1-002 | 空间平台资格分层 | 是 | COMPLETE | Visium/ST/Xenium/imaging/h5ad 等按平台分层 | `spatial_physical_units.tsv`、`REGISTRY_REPORT.md` | — | Stage 3 adapter 输入 |
| S1-003 | 数据角色与 claim 隔离（D1） | 是 | COMPLETE | 184 条 claim×source，主角色＋选择暴露＋泄漏组 | `data_role_assignment.tsv` | S1-D01 | 用过的验证数据同 claim 下不恢复独立身份 |
| S1-004 | 重复与泄漏审计 | 是 | COMPLETE | 8,731 条血缘边；unknown 不进独立验证计数 | `duplicate_lineage.tsv` | — | 拆分与计数依据 |
| S1-005 | Stage 4–9 实际候选来源的身份链/时间点/治疗/endpoint/重复/暴露核对；库存/可运行/可识别分列 | 是 | COMPLETE | 见证据表“资格三列”；不可识别项单列 | `EVIDENCE_AND_CONFLICTS.md` | — | 下游按三列消费，不得混用 |
| S1-006 | 临床/纵向可识别性判定 | 是 | COMPLETE | 直接纵向空间 rewiring=NOT_IDENTIFIABLE；资格审计可结束、科学问题仍不可识别 | `EVIDENCE_AND_CONFLICTS.md` | — | 空间疗效/修复分支阻塞（见 NEXT） |

结果规则：资格审计 COMPLETE；科学问题按 SUPPORTED/NOT_SUPPORTED/INCONCLUSIVE/NOT_IDENTIFIABLE 单列（本阶段多为 NOT_IDENTIFIABLE 的资格边界，非生物学否定）。
范围变更：无。
reviewer/date：待第三方实际审查时填写（见文末模板）。

## 偏移表

| deviation_id | original | revised | reason | benefit_evidence | claim_loss | exposure_change | decision | decision_by | date |
|---|---|---|---|---|---|---|---|---|---|
| S1-D01 | 给整个项目永久贴数据角色标签 | 按具体 claim 管理角色；用过的验证数据同 claim 下不恢复独立身份 | 同一数据在不同 claim 下角色不同；防止选择后验证 | 预期：角色精确；证据：`data_role_assignment.tsv:selection_exposure` 全量 `not_run_v7_stage1_metadata_audit`＋role_source 分列 | 无；secondary role 仅 support | 选择暴露逐 claim 记录 | ACCEPTED | repair-audit | 2026-09-11 |

## 第三方审查记录模板

| 项目 | 内容 |
|---|---|
| scope_version / run_id | plan_v7.0_2026-08-31 / v7-stage1-registry-1.0 |
| reviewer / date / independence | 待填写（须未承担本项实现与结果选择） |
| inputs and outputs inspected | OUTPUT_MANIFEST.yaml 全部路径＋hash；重点抽查 logical/physical/crosswalk/completeness/lineage/role 各表主键行 |
| computation reproduced | 命令见 OUTPUT_MANIFEST；确定性排序→byte-identical；第三方优先复算小判定表，不默认重跑全账本 |
| original requirements unresolved | S1-006 的 NOT_IDENTIFIABLE 项（原因＋结论影响见 EVIDENCE） |
| deviations reviewed | S1-D01 待 review |
| claim decisions | 资格角色（非生物学 claim）；科学 claim 为 NOT_RUN |
| stage closure / downstream permission | 资格审计 CLOSED；科学问题按 EVIDENCE 单列；下游见 NEXT_PHASE_READINESS |
| remaining actions | GSE211956 语义、Mendeley provenance、配对空间资产寻找（责任人/完成条件见 NEXT） |
