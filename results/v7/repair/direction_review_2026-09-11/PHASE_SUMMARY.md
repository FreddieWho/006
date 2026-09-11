# Stage 7｜训练内分子方向与cross-fit评估

执行状态：COMPUTATION_COMPLETE（见 RUN_RECEIPT.json）；阶段科学闭合：OPEN。

7个既定D2主要特征；训练患者基线定义方向，患者时间点同折；199次标签置换重新拟合方向。

基线分子区分能力弱，结果是敏感性分析；不能称已验证repair。

这是2026-09-11审核后的新运行，历史根目录结果保留。输入来源、代码、输出hash见本目录RUN_RECEIPT.json；第三方先读此摘要，再读EVIDENCE_AND_CONFLICTS.md和DECISION_LOG.md。

复算命令（新目录，禁止覆盖本运行）：

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/v7/review/run_direction_review.py --output <new_output_directory>
```

完整审核报告：docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md。
