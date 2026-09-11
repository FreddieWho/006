# Stage 4｜空间模块低秩留出比较

执行状态：COMPUTATION_COMPLETE（见 RUN_RECEIPT.json）；阶段科学闭合：OPEN。

71捕获资格审计；54个Visium捕获进入39特征完整子空间；整队列留出、3随机初始化。K=2均收敛，高阶模型仍不稳定。

这检验模块向量可压缩性；不等于开放空间场、独立GT验证或topology增量。

这是2026-09-11审核后的新运行，历史根目录结果保留。输入来源、代码、输出hash见本目录RUN_RECEIPT.json；第三方先读此摘要，再读EVIDENCE_AND_CONFLICTS.md和DECISION_LOG.md。

复算命令（新目录，禁止覆盖本运行）：

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/v7/review/run_closure_remediation.py spatial_factors --output <new_output_directory>
```

完整审核报告：docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md。
