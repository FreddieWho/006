# v6.2.1 realignment 执行状态

> 当前阶段已于 2026-08-27 闭合：GSE301741 身份交接通过，但 Phase4A mid-level coverage 仍阻塞。请先阅读 `results/v6_2/scientific_engineering_realignment_v1/gse301741_handoff_repair_v1/PHASE_REPORT_20260827.md`。本文其余内容是早期 realignment 历史记录；其中预注册式 gate、response lock 和冻结状态不再约束新的探索分支。

更新时间：2026-08-21（Asia/Shanghai）

## 已完成

1. **G0 冻结**：D01–D08 已按导师用户同意的推荐选项批准；合同、response lock、Phase8 `NOT RUN` 和 Phase6 membership hash 均已冻结。
2. **WP5 response-blind 表征必要性评估**：比较 direct module、linear latent、hierarchical module 和线性 SRB Stage-A 候选。SRB 没有达到预注册的额外收益门槛，主干保留为透明的 `direct_module`。
3. **WP6a metadata-only 锚点评估**：审计 68 个 cohort 的癌种、治疗情境、endpoint/timepoint schema、label quality、环境角色和可用性元数据。未读取任何 response/outcome 数值。

## 当前科学阻塞

WP6a 的结论不是“无效”，而是**在 response lock 下无法识别共享 response 方向**。治疗情境和癌种组合不完全可交换，endpoint/timepoint schema 也不完全统一；因此技术上缺少一个预先定义、可比较的数值 estimand。

这会阻塞：

- WP6b 的数值 anchor 对比；
- shared-direction/transport 结论；
- 后续 joint response/context model。

## 需要导师决定的事项

推荐选项是 `context_stratified_estimand_v1`：在治疗情境×癌种层内比较，跨情境只作 support-only，不把非等价 endpoint 池化。另需明确批准 WP6b 的一次性 response 解锁；WP6b 只做预注册数值对比，不重训表征。

备选是：

- `high_quality_endpoint_estimand_v1`：只纳入 endpoint/timepoint 元数据兼容且质量高的 cohort，换取更干净但更小的样本；
- `no_shared_direction_estimand_v1`：不估计共享方向，只报告各 cohort 探索性结果，科学边界最强但不能形成 pooled anchor claim。

在决定前，WP6b、joint response/context、Phase8 和任何疗效/barrier 结论保持未启动。

## 当前更新（WP6b-0 已执行）

导师已批准一次性 WP6b response 解锁。随后完成了限定范围的数值可行性审计，但因主环境复制不足停止：754 个 baseline 聚合行、报告为 58 个分析环境、2 个主锚候选、仅 1 个主环境达标、仅 HCC 一个主癌种、没有重复 endpoint×treatment 轴。没有拟合 `βe`、`μ` 或 `δe`，没有启动 WP7/Phase8。审计后已修正逐行资格过滤；该修正只可能减少资格，精确重跑需新的 response 解锁。

当前状态以 [WP6b 完成清单](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_completion_manifest.yaml) 和 [数据缺口规格](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_data_gap_specification.yaml) 为准；G0 冻结文件中的“解锁前”状态保留为历史快照，不作当前运行状态解读。
