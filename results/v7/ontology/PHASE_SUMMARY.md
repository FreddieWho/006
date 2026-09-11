# Stage 2 阶段摘要：测量与共同词汇

- 科学问题：哪些功能对象在什么模态和载体中可可靠测量？（节点计划 Stage 2）
- 范围版本：plan_v7.0；D2 仅 scRNA 测量证据，不外推到 bulk/spatial/perturbation。
- 现状：44 对象、88 count pseudobulk 分片、39 可执行特征；7 scRNA measurable、32 uncertainty；10 机制父轴中 1 全可测；exact-mask null 28 计算成功、11 不可估计。既有 `ONTOLOGY_AND_MEASUREMENT_REPORT.md`＋reliability/bridge contract；本次补齐五件套，不重跑全量计量。
- 方法：response-blind 词汇（8 legacy FM＋10 机制父轴/31 组件）；raw counts 重建 coarse/mid pseudobulk；GSE207422 scRNA–bulk 桥、GSE193736 扰动技术桥；exact-mask null、leave-study-family-out、coverage、ambient proxy、carrier、method agreement、confounding 诊断定 D2。
- 主要结果：逐特征×模态判定可测/有限/拒绝；共同词汇接口可完成；Stage 7 pre/post 使用可比尺度且未回归掉治疗信号（见证据表 Stage7 核对行）。
- 负结果/限制：bulk 仅 P05/P08 精确配对（描述性，39 特征均为 `paired_descriptive_insufficient_n`）；扰动仅重复性接口（median Spearman 0.881），不作方向结论；真实空间投影 `NOT_RUN_PENDING_STAGE3`（按合同进入 Stage 3）；ambient proxy 仅观察性；单 scRNA platform 下 platform 效应不可估计。
- 判定：共同词汇接口 COMPLETE；模态桥按三类输入交付（主要/敏感性/仅注释），无法校准的桥单列未完成；不存在“scRNA D2 自动推广到 bulk/spatial”。
- 下游用途：Stage 4–7 按三类输入消费；39 特征覆盖审计保留，不得静默删除失败特征。
