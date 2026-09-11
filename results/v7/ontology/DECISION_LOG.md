# Stage 2 决策记录

## 必需任务清单

| requirement_id | original_requirement | mandatory | execution_status | result | evidence_locator | deviation_id | downstream_impact |
|---|---|---|---|---|---|---|---|
| S2-001 | 经验程序基线（8 FM coverage/carrier/ambient/cohort dependence；matched null；leave-cohort-out） | 是 | COMPLETE | 44 对象/88 分片；null 28 成功、11 不可估计；不补零 | `diagnostics/`、`measurement_reliability.tsv` | — | D2 分级依据 |
| S2-002 | 机制轴补齐（≥10 轴；curated/经验/组合显式；不按 response 挑成员） | 是 | COMPLETE | 10 父轴/31 组件；response-blind | `module_dictionary.tsv`、`module_membership.tsv` | S2-D01 | 父轴按多组件向量解释 |
| S2-003 | Cell-state ontology（coarse/mid/fine；保留原始 label/confidence/unknown/mixed/low-quality） | 是 | COMPLETE | coarse/mid/fine 映射完成 | `cell_state_ontology.tsv`、`cell_state_label_mapping.tsv` | — | 跨数据映射依据 |
| S2-004 | 跨模态投影（cell/spot/bin/region＋patient 分数；coverage/platform/resolution/uncertainty） | 是 | COMPLETE | 跨模态映射＋三类输入；空间真实投影按合同进 Stage 3 | `cross_modal_mapping.tsv`、`STAGE3_MEASUREMENT_CONTRACT.yaml` | — | 下游按三类消费 |
| S2-005 | 每个主分析特征连接覆盖率/carrier/null/跨研究稳定性/score 语义 | 是 | COMPLETE | 逐特征行已连接（见证据表） | `measurement_reliability.tsv` | — | 失败链条定向修复依据 |
| S2-006 | 核对 Stage 7 pre/post 可比尺度且未回归掉治疗信号 | 是 | COMPLETE | Stage 7 使用同一 Stage 2 分数 schema；仅去技术 nuisance，保留 arm/time/arm×time/pairing（见 Stage 7 DECISION_LOG） | `EVIDENCE_AND_CONFLICTS.md:Stage7核对` | — | Stage 7 输入合法 |
| S2-007 | 具体失败链条定向修复 | 是 | COMPLETE | 保留 39 特征覆盖审计；失败特征不删除，只降级为注释/敏感性 | `measurement_reliability.tsv:D2_class` | S2-D01 | 完成度不得靠删特征改善 |

结果规则：逐特征×模态判定 可测/有限/拒绝；无“scRNA D2 自动推广”。
范围变更：无。
reviewer/date：待第三方实际审查时填写（见文末模板）。

## 偏移表

| deviation_id | original | revised | reason | benefit_evidence | claim_loss | exposure_change | decision | decision_by | date |
|---|---|---|---|---|---|---|---|---|---|
| S2-D01 | 单值机制分数/删失败特征 | 父轴多组件＋平台子空间＋abstention；保留 39 特征审计 | 方向不同的生物学不应压成单值；删特征美化完成度 | 预期：语义保真；证据：`measurement_reliability.tsv` 32 uncertainty＋11 null 不可估计均保留 | 无科学 claim 损失；部分特征降级为敏感性/注释 | 无数据暴露变化 | ACCEPTED | repair-audit | 2026-09-11 |

## 第三方审查记录模板

| 项目 | 内容 |
|---|---|
| scope_version / run_id | plan_v7.0 / STAGE2_RUN_MANIFEST.yaml |
| reviewer / date / independence | 待填写（须未承担本项实现与结果选择） |
| inputs and outputs inspected | OUTPUT_MANIFEST.yaml 路径＋hash；抽查 reliability 行的 coverage/null/稳定性/score 语义 |
| computation reproduced | 命令见 OUTPUT_MANIFEST；第三方优先复算小判定表，不默认重跑全量计量 |
| original requirements unresolved | 空间真实投影（按合同进 Stage 3，非本阶段缺口） |
| deviations reviewed | S2-D01 待 review |
| claim decisions | 逐特征×模态 D2（非疗效/机制 claim） |
| stage closure / downstream permission | 词汇接口 CLOSED；模态桥未校准部分单列未完成；下游见 NEXT |
| remaining actions | 无（失败链条修复已纳入保留策略） |
