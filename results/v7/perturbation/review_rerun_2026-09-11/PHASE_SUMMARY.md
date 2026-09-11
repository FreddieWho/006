# Stage 8｜真实扰动条件测量对比

执行状态：COMPUTATION_COMPLETE（见 RUN_RECEIPT.json）；阶段科学闭合：OPEN。

GSE193736的24个样本、156个LTBR-labelled对tNGFR-labelled条件对比，保持lineage与rest/stim分层。

样本replicate不自动等于独立donor；无PD1X靶点兼容链和外部验证。

这是2026-09-11审核后的新运行，历史根目录结果保留。输入来源、代码、输出hash见本目录RUN_RECEIPT.json；第三方先读此摘要，再读EVIDENCE_AND_CONFLICTS.md和DECISION_LOG.md。

复算命令（新目录，禁止覆盖本运行）：

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/v7/review/run_perturbation_review.py --output <new_output_directory>
```

完整审核报告：docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md。
