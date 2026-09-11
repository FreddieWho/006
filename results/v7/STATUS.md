# v7 scientific status

## 当前审核结论（2026-09-11）

状态：`V7_REVIEW_REMEDIATION_PARTIAL_SCIENTIFIC_OPEN`。

上一轮“修复完成”审核不通过：完成的是文档审计，核心重算多项未执行。本轮已实际补修，完整科学方案仍OPEN。

主入口：[审核报告](../../docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md)。

- Stage4：患者平衡原语、NMF共同RMSE、整队列留出、3初始化与子空间比较已运行；71捕获/2769项置换完成，其中2752项各999次、17项不可估计；真实topology增量与独立结构验证未闭合。
- Stage5：125患者、5队列拆6治疗环境；共同患者折、训练折基线、24组配对增量已重算，未见稳定正增量；独立验证仍缺。
- 空间response修正：GSE238264原始GEO phenotype核验7人4R/3NR；已补5原语的post-only横断面探索，不能称基线预测、独立验证或纵向repair。
- Stage6：独立HCC复现及可识别混杂分解未完成。
- Stage7：47配对患者的整向量置换、fraction、组成调整和cross-fit分子方向已补；可信barrier、独立组合及空间repair未闭合。
- Stage8：GSE193736 24样本、156个真实条件对比已执行；PD1X靶点/载体兼容链和外部验证未闭合。
- Stage9–10：多机制类外推、独立第三方复核和全方案闭合仍未完成。

新结果位于各阶段review/direction/post_spatial运行目录，以本轮RUN_RECEIPT与节点五件套为准；具体路径见审核报告。常规补修已获授权，不等待重复确认。

## 以下是上一轮历史记录，已被本轮审核结论替代

下文的“修复执行完成”“待确认”“偏离低”和旧科学状态只作provenance，不能作为当前执行依据。

## 阶段闭合修复执行状态（2026-09-11，给人读）

已完成：按修复方案把 Stage 0–3 补成完整五件套；给 Stage 4/5/6/7 各建了修复审计目录并写清哪些算完、哪些没算；把以前 tablas 里容易误读的地方改了名字（比如空间增量那一列现在明确只是组成调整前后的参考值，不是真正的增量证据）；Stage 8/9/10 的架子（证据卡、冲突表、最终报告草稿）已搭好。

正在做：等你确认这套“先把缺口摆到台面上”的修法是否接受；确认后才排真正的重算（共同患者折重跑、高精度置换、held-out 比较等）。

卡在哪里：真正的空间增量比较、NMF 的留出验证、高精度置换、共同折重跑、独立肝癌复现、修复方向、扰动方向和外部验证都还没跑；X 类别也不够做跨类预测。不是不想跑，是这些要么需要新开运行目录重算，要么缺合格的配对数据。

准备怎么解决：缺口都写进了每个阶段的任务表和“解除条件”里；有配对数据就补独立验证，没有就保持“不可判定”并守住结论边界，不硬凑阳性结果。

```
2026/9/11
ROADMAP  [#######---] 7/10 节点有可运行产物（0–7），8–10 为资格审计＋骨架
本周投入  科学问题 ████░░░░░░ 40%   基础设施 ██████░░░░ 60%

偏离程度  低
偏离位置  修法本身：原计划隐含“补齐＝重跑”，本轮以审计＋显式 NOT_RUN 代替部分重跑；
         已在各任务表单列未运行项，未改写历史，未放宽证据标准。
建议     接受“审计先行、重算排期”的节奏，按 NEXT 中的解除条件逐个销项。
```

## 给 agent 的接手信息（十行内）

- 活跃节点：修复执行完成待确认；科学主线仍停在 S7 分子完成／空间阻塞，原状态未改。
- 核心文件：`results/v7/evidence_cards/claim_matrix.tsv`、`results/v7/FINAL_SCIENTIFIC_REPORT.md`、各 `repair_2026-09-11/` 审计目录。
- 复现审计表：见各 DECISION_LOG 末尾复核命令（纯 pandas 派生，不重跑大矩阵）。
- 最近 DECISIONS：S4-D01 接受（Moran 差值改名）、S5 折缺陷定位、S9 资格失败停 held-out（详见各 DECISION_LOG）。
- 下一步：确认修法→按解除条件排重算；找配对空间＋response 资产；约独立复核人。

---

以下为原科学状态（本轮未改科学结论，只加修复审计层）：

状态：`V7_STAGE7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED`

> 2026-09-01 已修复上述问题：matched-null 保留 exact observation mask，并将 D2 限定为 scRNA 测量证据；bulk 和 perturbation 仅保留各自的描述性/重复性证据。

- v7 plan/roadmap：完成，科学一致性复核 `PASS`；
- v6/013 继承审计与仓库迁移：完成 Stage 0 范围；
- v7 registry：`STAGE1_COMPLETE`（跨仓库逻辑/物理单位、患者—区块—切片、治疗/疗效完整性、重复血缘和 claim-specific role 已生成并校验）；
- v7 ontology/measurement：`STAGE2_COMPLETE`（44 个对象、88 个 count pseudobulk 分片、39 个可执行特征；scRNA D2 中 7 个 `measurable`、32 个 `joint_representation_with_uncertainty`；10 个机制父轴中 1/9；exact-mask null 28 个计算成功、11 个明确不可估计）；
- v7 spatial foundation：`STAGE3_COMPLETE_WITH_LIMITATIONS`（8 个逻辑 pilot 单元、11 个物理捕获、7 个患者身份包；Visium/Xenium；D3 `PASS_WITH_LIMITATIONS`）；
- v7 spatial discovery：`S4_COMPLETE_WITH_LIMITATIONS`（71 个物理捕获、43 位患者、3 个数据集；response-blind；5 个跨数据集空间结构候选；2,752/2,769 个 topology 记录可估计；NMF 324 条低秩拟合可估计）；
- v7 clinical anchor/context：`S5_BASELINE_COMPLETE_WITH_LIMITATIONS`（125 位患者、5 个队列；患者级 pre 基线、bootstrap 和 response permutation 已完成）；
- v7 context/HCC decomposition：`S6_COMPLETE_WITH_LIMITATIONS`（5 个 response-blind 锁定候选；15 个 dataset/platform 实现、5 个 HCC residual 描述和 25 个临床注释组合）；
- v7 repair：`S7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED`（TASK01/TASK02 共 47 位配对患者、3,432 条患者内位移记录；response/pairing permutation 已完成；空间重排与 surrogate 不可识别）；
- v7 PD1+X repair：`S7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED`（分子/细胞位移已完成；空间重排未识别）；
- direct longitudinal spatial PD1+X rewiring：`NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`；
- v7 perturbation/external validation：`NOT_RUN`；
- v7 X-class generalization：`NOT_RUN`；
- v7 biological claims：支持 response-blind 词汇/测量层结论、跨数据集空间邻接候选和环境内分子/细胞临床锚；不把 marker proxy 当真实 composition，不把 bulk 两个 exact pairs 当成特征验证，不把 perturbation 重复性当成方向证据；尚不支持空间屏障、PD-1 failure 或 PD1+X repair 结论。

Stage 1 账本规模：1,670 个 `/006` 患者主表行、10,166 个样本主表行、1,002 个 `/013_spatial` physical rows（其中 167 个 resolved included candidate），以及 1,044 个合并空间物理单位行。Stage 1 仅是元数据审计，不产生生物学结论。

Stage 2 计量规模：32,869 个表达单元、1,281,891 条表达单元×特征记录、1,876 个患者—时间点；8 个 legacy FM 通过独立的 topic mass、native score、top-gene sensitivity 和 library-size 公式重放。TASK01、GSE207422_sc、GSE301741 的全细胞技术面板共 234 项，cell-mean 与 count-pseudobulk rank concordance 中位数为 0.758。GSE207422 仅 P05/P08 精确配对；39 个特征的患者级配对状态均为 `paired_descriptive_insufficient_n`；GSE193736 返回 24 个扰动样本的重复性接口，中位 Spearman 为 0.881，不作靶点方向结论。

关键边界：Lambrecht HCC 的旧摘要（25 pre/22 paired）与有效主表（44 patients/112 samples）不一致，已在报告中保留为冲突；GSE238264 为 7 名 HCC post-only PD1+cabozantinib 空间样本；GSE291246 为 35 个 BCC Xenium sections、17 个 QC 可重放记录且无 response。直接纵向 spatial PD1+X rewiring 仍为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

下一节点：启动 Stage 8 真实扰动、target/carrier reachability 与外部患者证据；空间→临床直接桥接仍为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`，需继续寻找 response-linked paired spatial 资产。直接空间疗效关联、PD1+X spatial repair 与外部验证仍为 `NOT_RUN`。
