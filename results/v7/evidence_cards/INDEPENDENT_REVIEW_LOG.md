# 独立复核记录（Stage 10）

## 2026-09-11 本轮审核记录

- 原提交审核者：本轮Codex会话；对上一轮修复产物进行代码、表、来源与收据核对。结论：原修复未完成，不通过全方案验收。
- 查出遗漏：混合治疗环境、不同患者折、全样本患病率基线、逐特征配对置换、NMF误差尺度及留出缺口；纠正GSE238264无response的错误断言。
- 补修与实际重算：见docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md及新运行RUN_RECEIPT。
- 当前claim判定入口改为review_claim_matrix.tsv；原claim_matrix.tsv和C*.yaml保留为上一轮审计历史，冲突时不得优先于本轮审核结果。
- 独立性边界：审核者实施了本轮补修；新增代码仅声明针对性自测、真实重算和输出检查。对新实现的第三方独立复核仍NOT_RUN，不能签署“独立复核通过”。
- 完整方案科学闭合：OPEN；文档和计算补修不消除原方案的独立结构/临床/组合及空间纵向缺口。

## 原模板（历史，待新第三方审查者补充）

状态：模板已立，无署名复核（作者自查不替代独立审查）。

| 项目 | 内容 |
|---|---|
| scope_version | plan_v7.0_2026-08-31 |
| reviewer / date / independence | 待填写（须未承担该项实现与结果选择） |
| inputs and outputs inspected | `evidence_cards/claim_matrix.tsv` 全表＋各 `C*.yaml` 证据定位路径＋hash（见各阶段 OUTPUT_MANIFEST） |
| computation reproduced | 第三方优先复算：NMF 可比 RMSE 透视表、单环境配对增量表、单候选实现表（三选一即可）；命令见各阶段 DECISION_LOG |
| original requirements unresolved | C03/C04/C09 NOT_IDENTIFIABLE；C10/C11 NOT_RUN；C12 资格失败；C05/C07/C08 仅探索性/描述性 |
| deviations reviewed | S0-D01/D02、S1-D01、S2-D01、S3-D01、S4-D01/D02、S5-D01、S6-D01、S7-D01、S8-D01、S9-D01（接受/拒绝/证据待补及理由待填） |
| claim decisions | 见 claim_matrix.tsv（支持／不支持／不确定／不可识别及证据定位） |
| stage closure / downstream permission | 修订范围内文档/审计 CLOSED；原科学范围 Stage 4/5/6/7 部分 OPEN、Stage 8/9/10 科学 OPEN；下游按各 NEXT 执行 |
| remaining actions | 各 NEXT 的 unblock_condition；责任人待指定 |

已签署的 review 到达后追加行记录，不得以签署掩盖缺失计算。
