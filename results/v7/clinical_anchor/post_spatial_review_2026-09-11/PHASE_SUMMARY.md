# Stage 5｜治疗后空间response探索

执行状态：COMPUTATION_COMPLETE（见 RUN_RECEIPT.json）；阶段科学闭合：OPEN。

原始GEO phenotype验证7位空间患者4R/3NR，分析5个事先定义的proxy，保留全部结果。

post-only开发支持；无基线/纵向空间，也不证明composition之外的topology增量。

这是2026-09-11审核后的新运行，历史根目录结果保留。输入来源、代码、输出hash见本目录RUN_RECEIPT.json；第三方先读此摘要，再读EVIDENCE_AND_CONFLICTS.md和DECISION_LOG.md。

复算命令（新目录，禁止覆盖本运行）：

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/v7/review/run_spatial_response_review.py --output <new_output_directory>
```

完整审核报告：docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md。
