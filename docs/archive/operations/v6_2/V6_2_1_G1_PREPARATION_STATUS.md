# G1 / WP6b 当前状态：数值可行性审计已停止

更新时间：2026-08-22（Asia/Shanghai）

## 已完成

- G1 共享 estimand 与一次性 WP6b response 解锁已由导师分别批准；解锁验证器通过，随后已消费该一次性授权。
- 在固定的 8-FM `direct_module` 表征、固定 join key 和预注册环境定义下，完成 WP6b-0 数值覆盖审计。
- 本阶段只读取允许的 response label/质量/环境字段用于类别覆盖和环境可行性判断；没有拟合环境效应 `βe`、共享向量 `μ` 或情境偏移 `δe`。

## 审计结果

`WP6B_STOPPED_FEASIBILITY`：754 个 baseline 聚合行、58 个分析环境中，报告为 2 个主锚候选，且仅 1 个达到主环境门槛（TASK01，25 例，13 responder / 12 failure）。另一个 clean 候选 GSE301741 仅 2 例（1/1）。符合主门槛的癌种只有 HCC，未形成重复的 endpoint×treatment 轴。该计数来自首轮环境汇总；随后已修正逐行资格过滤，修正只会减少候选资格，不能把本轮 STOP 变成 PASS；精确重跑需新的 response 解锁。

因此不能进入 WP7，也不能报告 shared direction、context modulation、barrier attribution 或 Phase8 结果。该停止是数据可比性/复制不足，不是模型拟合失败。

## 当前阻塞与技术含义

需要至少再有一个独立、可比较的 clean PD1/PD-L1 baseline 环境（样本数 ≥15、每类 ≥5），并最好来自第二癌种或同一 endpoint×treatment estimand 的非混杂复制。support-only、双免疫联合、高混杂环境不能为达到门槛而并入 clean anchor，也不能跨 endpoint 直接池化。

另外，TASK01 的绑定记录标为 `temporal_sensitivity_only`；在任何 WP7 效应模型前还需明确该时间点边界是否可以作为正式 anchor。

## 证据文件

- [WP6b 完成清单](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_completion_manifest.yaml)
- [数值可行性摘要](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_numeric_feasibility/wp6b_feasibility_summary.yaml)
- [环境覆盖表](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_numeric_feasibility/anchor_environment_feasibility.csv)
- [数据缺口规格](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_data_gap_specification.yaml)
- [一次性解锁合同（已消费）](../../results/v6_2/scientific_engineering_realignment_v1/wp6b_response_unlock_contract_v1.yaml)
- [当前执行状态](../../results/v6_2/scientific_engineering_realignment_v1/current_execution_state_v1.yaml)

`g1_preparation_manifest.yaml` 保留为解锁前授权快照；当前真值以 WP6b 完成清单和数据缺口规格为准。
