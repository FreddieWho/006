# v7 阶段闭合修复执行记录

审核后说明：本文件记录的是上一轮文档/派生审计，不能证明原方案修复完成。“待确认后重算”的停止条件已撤销；本轮补修及真实运行见[审核报告](STAGE_CLOSURE_REVIEW_2026-09-11.md)。

日期：2026-09-11

性质：本次修复的执行收据（审计先行）。与 `STAGE_CLOSURE_CHANGELOG_2026-09-11.md`（方案交付）的区别：本文件记录实际落盘的审计产物、任务状态与未运行项。

依据：`docs/plan/v7/stage_closure_repair_plan.md`、`docs/plan/v7/stage_closure_nodes.md`。

## 落盘清单（新运行目录＋五件套补齐）

- 根：`TODO.md`（新建，共读入口＋分支记录）。
- Stage 0：`results/v7/inheritance/` 五件套（PHASE_SUMMARY/OUTPUT_MANIFEST/DECISION_LOG/EVIDENCE/NEXT）。
- Stage 1：`results/v7/registry/` 五件套（复用既有运行，不重跑）。
- Stage 2：`results/v7/ontology/` 五件套（复用既有运行）。
- Stage 3：`results/v7/spatial_foundation/` 五件套（复用既有运行）。
- Stage 4：`results/v7/spatial_discovery/repair_2026-09-11/`（AUDIT＋nmf_scale_audit＋seed_stability＋comparable_pivot＋proxy_by_feature＋surrogate_qualification＋REPAIR_MANIFEST）；DECISION_LOG/EVIDENCE/NEXT/PHASE_SUMMARY 已更新（未运行项明示）。
- Stage 5：`results/v7/clinical_anchor/repair_2026-09-11/`（AUDIT＋paired_increment＋calibration＋REPAIR_MANIFEST）；DECISION_LOG/EVIDENCE/NEXT 已更新。
- Stage 6：`results/v7/context/repair_2026-09-11/`（AUDIT＋hcc_qualification）；DECISION_LOG/EVIDENCE/NEXT 已更新。
- Stage 7：`results/v7/repair/repair_2026-09-11/`（AUDIT＋direction_gap）；DECISION_LOG/EVIDENCE/NEXT 已更新。
- Stage 8：`results/v7/perturbation/` 七件套（含主报告＋两占位表）、`results/v7/external_validation/` 七件套（含两占位表）。
- Stage 9：`results/v7/repair/stage9_x_class/` 七件套（资格失败，停 held-out）。
- Stage 10：`results/v7/evidence_cards/`（12 YAML＋claim_matrix＋conflict_register＋REVIEW_LOG）、`results/v7/figures/FIGURE_MANIFEST.yaml`、`results/v7/FINAL_SCIENTIFIC_REPORT.md`。
- 状态：`results/v7/STATUS.md` 顶部追加修复执行状态（两部分写法）＋偏离评估；原科学状态保留未改。

## 替换与重跑历史

- 无覆盖：原结果文件均保留；审计产物进新运行目录；五件套更新为追加式（原决策表保留）。
- 无科学重算：本轮审计表均为 pandas 派生既有表；未重跑 71 捕获空间计算与患者级模型。
- 脏树 provenance：HEAD 8bdb424＋2026-09-11 脏/未跟踪摘要记入各 OUTPUT/REPAIR_MANIFEST；发布状态记 NOT_VERIFIED_FROM_LOCAL_RESULTS。

## 验证

- Markdown 路径检查：待执行（见本任务收尾）。
- 字段检查：各 TSV 主键/状态值见 OUTPUT_MANIFEST；第三方优先复算小表。
- 代码测试：纯文档/审计调整，仅做路径与字段检查；未改 `src/` 行为代码，故不触发全量回归（Stage 5 折缺陷仅定位未改码）。

## 下一动作

待确认“审计先行”修法后，按各 NEXT 的 unblock_condition 排真重算；独立复核人到位后签署 `evidence_cards/INDEPENDENT_REVIEW_LOG.md`。
