# v7 scripts

本目录只放可直接重放的 v7 命令入口。可复用实现进入 `src/v7/`，测试进入 `tests/v7/`。

计划子目录：

- `registry/`
- `ontology/`
- `spatial/`
- `clinical/`
- `repair/`
- `perturbation/`
- `evaluation/`

Stage 1 入口：`python scripts/v7/registry/build_stage1_registry.py`。它生成 `results/v7/registry/` 下的七张账本、报告和运行 manifest；当前仍未运行 v7 生物学模型。
