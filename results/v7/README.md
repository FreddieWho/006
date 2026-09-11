# v7 results

本目录是v7的唯一结果入口。当前以[STATUS](STATUS.md)和[修复审核报告](../../docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md)为准：基础层与多项补修计算已执行，完整科学闭合仍OPEN。Stage6/7已有首轮与补修结果，Stage8已补真实条件对比，不能继续统称“Stage6及以后未运行”。

审核后重算入口位于`scripts/v7/review/`；各新运行目录的RUN_RECEIPT及五件套记录精确输入、代码、结果和未完成任务。`review_2026-09-11/artifact_verification.json`是本轮收据核对入口，`evidence_cards/review_claim_matrix.tsv`是当前claim矩阵。下列原阶段命令保留为历史可重算入口，不能覆盖或替代新运行。

未来子目录按 `docs/plan/v7/roadmap.md` 创建：

- `registry/`
- `ontology/`
- `spatial_foundation/`
- `spatial_discovery/`
- `clinical_anchor/`
- `context/`
- `repair/`
- `perturbation/`
- `external_validation/`
- `evidence_cards/`
- `figures/`

大型矩阵、图像和 checkpoint 不进入普通 Git；以 manifest、SHA-256 和可重算入口引用。

Stage 1 重放：`python scripts/v7/registry/build_stage1_registry.py`。该入口只读取元数据和控制面文件，不读取表达矩阵或图像。

当前阶段入口：

- Stage 4：`python scripts/v7/spatial/run_stage4.py --permutations 19 --nmf-sample-limit 1500 --config config/v7/stage4.yaml`；结果在 `spatial_discovery/`。
- Stage 5：`python scripts/v7/clinical/run_stage5.py`；结果在 `clinical_anchor/`。
- Stage 6：`python scripts/v7/context/run_stage6.py`；结果在 `context/`。
- Stage 7：`python scripts/v7/repair/run_stage7.py`；结果在 `repair/`。

四阶段均已保留 `PHASE_SUMMARY.md`、`OUTPUT_MANIFEST.yaml`、`DECISION_LOG.md`、`EVIDENCE_AND_CONFLICTS.md` 和 `NEXT_PHASE_READINESS.yaml`。空间结果为 response-blind 结构候选，临床结果为环境内患者级探索性基线，Stage 6 为 context 描述和桥接审计，Stage 7 为分子/细胞 repair 基线；空间分支仍不可识别，所有阶段都不单独支持因果或 PD1+X 空间修复结论。
