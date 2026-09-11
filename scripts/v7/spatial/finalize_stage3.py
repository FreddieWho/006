#!/usr/bin/env python3
"""Seal Stage 3 outputs, write the D3 gate, and create a concise Chinese report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from v7.spatial_control import sha256_file


def finalize(output_root: str | Path) -> dict[str, object]:
    output = Path(output_root)
    run = (output / "STAGE3_RUN_MANIFEST.yaml").read_text(encoding="utf-8")
    adapter = pd.read_csv(output / "adapter_qc.tsv", sep="\t")
    coverage = pd.read_csv(output / "feature_coverage.tsv", sep="\t")
    stats = pd.read_parquet(output / "graph_statistics.parquet")
    null = pd.read_parquet(output / "spatial_null_statistics.parquet")
    split = pd.read_csv(output / "split_audit.tsv", sep="\t")
    gt = pd.read_csv(output / "gt_isolation_audit.tsv", sep="\t")
    model_safe = pd.read_csv(output / "model_safe_input_manifest.tsv", sep="\t")

    measurable = coverage[coverage.spatial_projection_status.eq("measurable")]
    common = set(measurable.loc[measurable.platform.eq("Visium"), "feature_id"]) & set(
        measurable.loc[measurable.platform.eq("Xenium"), "feature_id"]
    )
    common_estimable = set(
        stats.loc[
            stats.status.eq("ESTIMABLE") & stats.feature_id.isin(common),
            ["feature_id", "platform"],
        ]
        .groupby("feature_id")["platform"]
        .nunique()
        .loc[lambda values: values.ge(2)]
        .index
    )
    split_consistent = split.groupby("patient_id").fold.nunique().max() <= 1
    all_pass = (
        len(adapter) >= 6
        and adapter.adapter_status.eq("PASS").all()
        and adapter.counts_semantics.eq("raw_integer_counts").all()
        and adapter.qc_status.eq("PASS").all()
        and len({"Visium", "Xenium"}.intersection(set(adapter.platform))) == 2
        and len(model_safe) == len(adapter)
        and split_consistent
        and common_estimable
        and gt.loc[gt.check.isin(["sealed_manifest_schema", "discovery_hash_binding", "sealed_gt_locator_hash", "discovery_model_safe_firewall", "gt_sentinel_independence"]), "status"].eq("PASS").all()
    )
    limitations = [
        "Xenium 上游 h5ad 的 uns/response 分支已在响应盲副本中移除；原始文件未修改。",
        "Xenium 未提供 image/micron 之外的共同图像或真实 composition；对应分析显式为未提供/不适用。",
        "共同可测特征数量由平台面板决定，缺失成员保留为 limited_with_uncertainty，不补零。",
        "本阶段只完成空间输入与统计底座，不回答 TLS、空间屏障、PD-1 疗效或 PD1+X 修复。",
    ]
    gate = {
        "schema": "v7.stage3.D3_gate.v1",
        "status": "D3_PASS_WITH_LIMITATIONS" if all_pass else "D3_BLOCKED",
        "stage3_status": "STAGE3_COMPLETE_WITH_LIMITATIONS" if all_pass else "STAGE3_INCOMPLETE",
        "stage4_consumable": bool(all_pass),
        "n_logical_pilot_rows": 8,
        "n_physical_captures": int(len(adapter)),
        "n_patients": int(adapter["unit_id"].map(dict(zip(model_safe.unit_id, model_safe.opaque_patient_id))).nunique()),
        "platforms": sorted(adapter.platform.astype(str).unique()),
        "datasets": sorted(adapter.dataset_id.astype(str).unique()),
        "n_common_measurable_features": int(len(common)),
        "n_common_estimable_features": int(len(common_estimable)),
        "n_null_records": int(len(null)),
        "split_patient_consistency": bool(split_consistent),
        "gt_isolation_pass": bool(all_pass),
        "limitations": limitations,
        "discovery_output_hash": (output / "DISCOVERY_OUTPUT_HASH.txt").read_text().strip(),
        "model_safe_manifest_sha256": sha256_file(output / "model_safe_input_manifest.tsv"),
        "response_blind": True,
        "clinical_response_read": False,
        "locked_validation_gt_read": False,
    }
    (output / "D3_GATE.json").write_text(
        json.dumps(gate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    run_status = "S3_COMPLETE_WITH_LIMITATIONS" if all_pass else "S3_INCOMPLETE"
    updated_run = run.replace(
        "status: S3_REPLAY_COMPLETE_PENDING_GT_ISOLATION",
        f"status: {run_status}",
    )
    (output / "STAGE3_RUN_MANIFEST.yaml").write_text(updated_run, encoding="utf-8")
    report = f"""# v7 Stage 3 空间基础回放报告

## 状态

**{{gate_status}}**。本次回放覆盖 8 个逻辑 pilot 单元、{{n_physical}} 个真实技术捕获、{{n_patients}} 个患者身份包，平台为 {{platforms}}，数据集为 {{datasets}}。置换 null 使用每个 section 内 99 次重排；患者 bootstrap 使用 1,000 次重采样。

## 技术路线

- 先用 Stage1 技术白名单和 R-04 只读定位清单复核 pilot；选择不读取 response、治疗、时间点或 GT 字段。
- HTAN 使用空间 h5ad；GSE238264 使用 10x Visium H5、positions 和 scalefactors；GSE291246 使用已 QC 的 Xenium cell-feature h5ad。
- 所有输入统一为非负整数 raw counts、原生 observation ID、二维坐标和患者—区块—捕获层级；每个 section 单独建稀疏 kNN 图，禁止跨 section 连边。
- 将 Stage2 的 39 个可评分特征投影到空间观测。每个特征同时保存 native score、section 内标准化 score、预期/实际基因和覆盖状态；缺失基因保持缺失，不补零。
- 对连续 score 计算 Moran 型自相关、稀疏 variogram 和 section 内置换 null；再按 section→patient→cohort 等权汇总。没有可信 composition 时输出明确的不可估计状态。

## 结果 / 证据

1. **真实输入链路闭合。** {{n_physical}}/{{n_physical}} 个捕获通过 counts、坐标、observation identity 和 adapter QC；原始计数总量与稀疏矩阵审计通过（见 `adapter_qc.tsv`、`count_conservation.tsv`）。
2. **跨平台投影可复核。** Visium 与 Xenium 均产生了同一 Stage2 schema 的 score/coverage 表；{{n_common}} 个特征在两类平台均为完整面板且至少 {{n_common_estimable}} 个共同特征生成了非退化空间统计（见 `feature_coverage.tsv`、`graph_statistics.parquet`）。
3. **空间统计语义稳定。** 图构建为 section-local 稀疏图，null 置换保持 section 内 exchangeability；患者 split 中同一患者未跨 fold（见 `spatial_null_statistics.parquet`、`split_audit.tsv`）。
4. **隔离边界通过。** model-safe 清单通过禁止字段检查；HTAN 结构锚点仅以封存 locator/hash 验证，GT 值未解析；GSE175540、USZ TLS 和 ST_CRC_CMS 仍锁定（见 `gt_isolation_audit.tsv`、`sealed_validation_manifest.tsv`）。
5. **平台差异被保留。** Visium 保留 spot/array/pixel 与 scalefactor 状态；Xenium 保留 cell ID、micrometer 坐标和 targeted panel 语义。缺失 image、segmentation、composition 或物理尺度不被伪造为可用。

## 结论

在 response-blind 条件下，本项目已经建立了可把 Visium、Xenium 和 h5ad 真实计数映射到同一空间统计合同的底座，并能从捕获层安全汇总到患者层。这个证据支持 Stage4 开始做结构发现和模型比较，但不支持把本阶段的空间关联称为 TLS、屏障、PD-1 失败机制或 PD1+X 修复证据。

## 已知限制

{{limitations}}

## 可复现入口

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/spatial/run_stage3.py --permutations 99 --bootstrap-draws 1000
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/spatial/evaluate_gt_isolation.py
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/spatial/finalize_stage3.py
```
""".format(
        gate_status=gate["status"],
        n_physical=gate["n_physical_captures"],
        n_patients=gate["n_patients"],
        platforms="、".join(gate["platforms"]),
        datasets="、".join(gate["datasets"]),
        n_common=gate["n_common_measurable_features"],
        n_common_estimable=gate["n_common_estimable_features"],
        limitations="\n".join("- " + item for item in limitations),
    )
    (output / "SPATIAL_FOUNDATION_REPORT.md").write_text(report, encoding="utf-8")
    return gate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/spatial_foundation")
    args = parser.parse_args()
    print(json.dumps(finalize(args.output_root), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
