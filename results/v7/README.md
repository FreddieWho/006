# v7 results

本目录是 v7 的唯一结果入口。Stage 1 已完成元数据 registry；后续科学计算仍未运行。

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
