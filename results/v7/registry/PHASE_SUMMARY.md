# Stage 1 阶段摘要：身份、数据角色和可识别性

- 科学问题：哪些资产支持哪个 claim，真正独立样本数是多少？（节点计划 Stage 1）
- 范围版本：plan_v7.0_2026-08-31；本轮为元数据审计范围（未打开表达矩阵/图像/h5ad 内容，未运行科学模型）。
- 现状：`REGISTRY_REPORT.md`＋7 张账本表＋`STAGE1_RUN_MANIFEST.json`已存在并校验；本次补齐五件套，不重跑全量注册计算。
- 方法：跨 `/006` 与只读 `/013_spatial` 的 logical/physical registry；患者—区块—切片 crosswalk；治疗/response 完整性；重复血缘；claim×source 角色表（D1）。只针对 Stage 4–9 实际候选来源核对身份链、counts/坐标/模态、时间点、治疗/endpoint、重复关系和选择暴露；将“库存存在”“本地可运行”“临床/纵向可识别”分列。
- 主要结果：1,670 患者主表行、10,166 样本主表行、1,002 `/013` physical rows（167 RESOLVED_INCLUDED_CANDIDATE）、1,044 合并空间物理单位行、184 条 claim×source 角色记录；每个 config source 均有 dataset-level logical unit，角色外键闭合。
- 负结果/限制：直接纵向空间 PD1+X rewiring 为 NOT_IDENTIFIABLE（无同患者空间 pre/on/post 链）；GSE211956 response 语义待核验；Mendeley provenance 未闭合；GSE291246 无 response 只作 support；Lambrecht 旧摘要冲突保留。
- 判定：每个候选来源已获 discovery/validation/support/reject 角色；未知身份不增加独立患者。缺数据的资格审计可结束，但对应科学问题仍 NOT_IDENTIFIABLE（见 DECISION_LOG 任务表）。
- 下游用途：分别交付空间发现、结构验证、临床、配对、扰动、X-class 的合格资产清单与缺口；代理资格不额外要求 response（结构代理训练/验证 response-blind，临床关联才需要 response）。
