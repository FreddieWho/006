# v7 Stage 2：共同生物学词汇与跨模态测量

## 技术路线

- 固定继承 8 个 legacy FM，并并行建立 10 个机制父轴、31 个显式组件。
- 对 44 个 count-compatible 单细胞对象从 raw counts 重建患者/样本×时间点×coarse/mid state pseudobulk。
- 在 GSE207422 建立 response-blind scRNA–bulk 桥，在 GSE193736 建立扰动 bulk 技术桥。
- 在 TASK01、GSE207422_sc、GSE301741 对全部合格细胞执行逐细胞技术计分，并与 count pseudobulk 的排序一致性对照；不导出全细胞分数表。
- 通过 exact-observation-mask matched-expression null、leave-study-family-out、coverage、native ambient proxy、carrier、method agreement 和 confounding 诊断确定 scRNA 测量等级。null 候选不足时保留 `NOT_TESTABLE`，不把未测基因当作 0。
- 空间 spot/bin/region 仅输出 Stage 3 接口，本阶段未执行真实空间投影。

## 结果 / 证据

- 完整处理对象：44；pseudobulk 分片：88；计数守恒全部通过。
- 可执行测量特征：39 个（8 FM + 31 机制组件）；另有 10 个机制父轴按组件向量解释。
- 39 个可执行特征的 scRNA D2 分布：{"joint_representation_with_uncertainty": 32, "measurable": 7}；10 个机制父轴：{"joint_representation_with_uncertainty": 9, "measurable": 1}。D2 仅表示 scRNA 测量证据，不外推到其他模态。
- exact-mask null 诊断状态：{"NOT_TESTABLE_INSUFFICIENT_EXACT_MASK_NULLS": 11, "PASS_COMPUTED": 28}；具体不可估计的特征保留为不确定性，不强行使用。
- 8 个 legacy FM 均通过当前输入上的冻结公式重放。旧 v6 行级数值因上游 cell-to-context 分配更新而不再要求逐行相等，已作为 historical context drift 明示记录。
- 逐细胞技术面板完成队列：3；234 个“队列×层级×特征”比较的 rank concordance 中位数为 0.758，低或负一致性的组件保留为 D2 uncertainty。
- GSE207422 精确患者—时间点配对仅 P05、P08；P07 为时间点不一致的支持性桥接。P05/P08 的 39 维路线 Spearman 分别为 0.729 和 0.282，只作接口描述；39 个特征的患者级配对校准都是 `paired_descriptive_insufficient_n`。
- GSE193736 的 24 个扰动条件/重复样本均可使用同一测量接口；39 个特征的重复谱 median Spearman 为 0.881，但不据此作扰动方向或疗效结论。
- ambient 诊断使用 low-quality/ambient 标签作观察性 proxy，不等同于已证明污染因果；platform 效应在本阶段单一 scRNA platform 下不可估计。

## 结论

Stage 2 已建立可追溯、response-blind 的共同生物学词汇。scRNA 的 D2 测量基础已可供后续空间与患者级分析；bulk 仅有两位精确配对的描述性支持，perturbation 仅有重复性接口证据，不是已验证的等价模态。空间屏障、PD-1 失败和 PD1+X 修复仍需后续阶段的数据与分析。

## 风险与限制

- 机制父轴是多组件向量，不把方向不同的生物学强行压成单值。
- D2 uncertainty 是测量证据的一部分；不能把带不确定性的程序写成已确认机制。matched-null 的区间是算法稳定/置换包络，不是生物效应置信区间。
- 逐细胞先归一化与 counts 先汇总再归一化是不同估计量；技术面板只检验排序稳健性，不要求数值相等。
- GSE207422 的精确配对患者数很少，只支持接口可行性。
- 真实空间投影状态为 `NOT_RUN_PENDING_STAGE3`。
