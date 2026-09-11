# v7 Stage 5 患者级临床基线报告

## 状态

**S5_BASELINE_COMPLETE_WITH_LIMITATIONS**。本轮将疗效字段仅用于临床锚定，按治疗环境和癌种分开分析；没有读取空间结构候选，也没有把细胞当作独立患者。

## 数据与方法

- 纳入 125 位患者，队列为 GSE286827、GSE301741、LAMBRECHT_HCC、TASK01、TASK02。
- 预治疗样本按 T0/baseline/pre 口径选择；预治疗缺失、重复或冲突在患者表中保留。
- 主要分子模型使用 D2 `measurable` 特征；39 特征模型作为预先定义的敏感性分析。
- 细胞组成使用 scRNA pseudobulk 的 coarse-state 细胞数比例，按无 response 的总体丰度选择最多12个状态；这是患者级组成特征，不是疗效标签。
- 每个环境内部做患者分层交叉验证，报告 AUROC、PR-AUC 和 Brier；同时运行 response permutation 负对照。

## 结论边界

输出可以回答哪些分子/细胞状态在各治疗环境内与疗效共同变化，以及空间结果是否值得继续桥接。它不能单独证明 PD-1 因果机制，也不能证明 X 的增量疗效。TASK01/TASK02 的配对位移仅作描述；直接纵向空间重排仍是 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \
  python scripts/v7/clinical/run_stage5.py
```
