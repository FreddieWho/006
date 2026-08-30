# v7 scientific status

状态：`V7_STAGE1_COMPLETE_STAGE2_NEXT`

- v7 plan/roadmap：完成，科学一致性复核 `PASS`；
- v6/013 继承审计与仓库迁移：完成 Stage 0 范围；
- v7 registry：`STAGE1_COMPLETE`（跨仓库逻辑/物理单位、患者—区块—切片、治疗/疗效完整性、重复血缘和 claim-specific role 已生成并校验）；
- v7 ontology/measurement：`NOT_RUN`；
- v7 spatial foundation/discovery：`NOT_RUN`；
- v7 clinical anchor/context：`NOT_RUN`；
- v7 PD1+X repair：`NOT_RUN`；
- direct longitudinal spatial PD1+X rewiring：`NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`；
- v7 perturbation/external validation：`NOT_RUN`；
- v7 X-class generalization：`NOT_RUN`；
- v7 biological claims：无。

Stage 1 账本规模：1,670 个 `/006` 患者主表行、10,166 个样本主表行、1,002 个 `/013_spatial` physical rows（其中 167 个 resolved included candidate），以及 1,044 个合并空间物理单位行。Stage 1 仅是元数据审计，不产生生物学结论。

关键边界：Lambrecht HCC 的旧摘要（25 pre/22 paired）与有效主表（44 patients/112 samples）不一致，已在报告中保留为冲突；GSE238264 为 7 名 HCC post-only PD1+cabozantinib 空间样本；GSE291246 为 35 个 BCC Xenium sections、17 个 QC 可重放记录且无 response。直接纵向 spatial PD1+X rewiring 仍为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

下一节点：Stage 2 共同生物学词汇与跨模态测量。
