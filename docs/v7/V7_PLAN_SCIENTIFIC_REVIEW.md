# v7 plan 科学一致性复核

日期：2026-08-31
结论：`PASS`

复核重点与闭合结果：

1. GSE291246 已纠正为 BCC Xenium；Mendeley 为 11 个 HCC object；ST_CRC_CMS 为 7 patient、14 section、12 个可审计 source。
2. patient 是唯一外层独立单位；block、section、region、spot/cell 均为嵌套观测。
3. 数据角色按 claim 管理；参与 ontology、候选、阈值或模型选择的数据不得再承担同一 claim 的独立验证。
4. Stage 4B 在读取 response 前形成版本化 shared architecture，并以整患者留出验证 architecture surrogate。
5. Stage 6 只注释 Stage 4B 已定义的结构，不用 response 反向重定义 shared architecture。
6. 当前没有闭合的同患者纵向 spatial PD1+X 链；没有真实纵向 spatial 或通过独立验证的 surrogate 时，`architecture_rewiring=NOT_IDENTIFIABLE`，Branch C 阻塞。
7. R-04 继承路径、状态证据与 SHA-256 已核对；当前只能继承接口和诊断，不能继承有效 K 或潜在场结论。

该复核确认 v7 没有降低宏观目标，同时避免伪重复、选择后验证、跨模态拼接式过度主张和将 post-only spatial 误写为纵向 repair。
