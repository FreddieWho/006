# Stage 7｜患者向量置换及细胞组成敏感性

执行状态：COMPUTATION_COMPLETE（见 RUN_RECEIPT.json）；阶段科学闭合：OPEN。

47位配对患者、3432条分子位移；199次标签/完整向量配对置换；156项分子调整、40项细胞比例对比及区间。

独立组合验证、可信barrier资格、空间重排仍未完成；非因果解释。

这是2026-09-11审核后的新运行，历史根目录结果保留。输入来源、代码、输出hash见本目录RUN_RECEIPT.json；第三方先读此摘要，再读EVIDENCE_AND_CONFLICTS.md和DECISION_LOG.md。

复算命令（新目录，禁止覆盖本运行）：

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/v7/review/run_closure_remediation.py repair --output <new_output_directory>
```

完整审核报告：docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md。
