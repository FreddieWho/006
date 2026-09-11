# v7 Stage 7 PD1+X repair 基线报告

## 状态

**S7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED**。本轮在 TASK01/TASK02 中完成 3432 条患者内分子位移记录，覆盖 47 位患者和 T1/T2 两个目标时间点；空间重排分支保持不可识别。

## 结果 / 证据

- 以患者为单位计算 T1−T0、T2−T0 的 39 个 Stage 2 特征位移，并按 response 分层描述；没有把分子位移改写为空间屏障重排。
- `39` 个特征形成 response-stratified displacement 表；`repair_direction` 保持未指定，因为没有给每个分子特征预设 X-class 机制方向。
- 运行 response permutation 和 pairing permutation 负对照；结果用于判断配对/标签结构是否足以解释观察到的差异。
- `observed_spatial_rewiring.tsv` 与 `surrogate_inferred_architecture_rewiring.tsv` 均明确为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

## 结论边界

当前可以继续研究 PD-1 单药与 PD1+lenvatinib 的患者内分子/细胞状态位移，不能声称 X 的因果增量，也不能声称空间屏障被修复。要升级 spatial repair，需要同一患者的纵向空间数据，或已经在独立配对患者上验证通过的 Stage 4B surrogate。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \
  python scripts/v7/repair/run_stage7.py --response-permutations 199 --config config/v7/stage7.yaml
```
