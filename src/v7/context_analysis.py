"""Stage 6 context decomposition and spatial-to-clinical bridge audit.

This stage consumes the response-blind Stage 4 architecture table and the
patient-level Stage 5 clinical anchor.  Response is used only to annotate
already locked spatial candidates; it never re-selects, re-clusters, or
renames an architecture.  Without paired spatial/response patients the bridge
is recorded as not identifiable rather than inferred from separate cohorts.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .spatial_control import sha256_file


HCC_SPATIAL_DATASETS = frozenset({"GSE238264"})
HCC_CLINICAL_DATASETS = frozenset({"LAMBRECHT_HCC"})
CONTEXT_LABELS = {
    "GSE238264": "HCC_spatial_post_treatment",
    "GSE291246": "BCC_Xenium",
    "HTAN_VANDERBILT_CRC": "CRC_Visium",
}


def _safe_quantile(values: Iterable[object], quantile: float) -> float:
    array = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(float)
    array = array[np.isfinite(array)]
    return float(np.quantile(array, quantile)) if len(array) else np.nan


def _sign_label(values: Iterable[object]) -> tuple[str, float, float]:
    array = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(float)
    array = array[np.isfinite(array)]
    if not len(array):
        return "UNRESOLVED_NO_ESTIMATES", np.nan, np.nan
    positive = float(np.mean(array > 0))
    negative = float(np.mean(array < 0))
    if positive >= 0.8:
        label = "POSITIVE_DOMINANT"
    elif negative >= 0.8:
        label = "NEGATIVE_DOMINANT"
    elif max(positive, negative) >= 0.6:
        label = "CONTEXT_MODULATED"
    else:
        label = "CONFLICTING_OR_NEAR_ZERO"
    return label, positive, negative


def _read_stage4(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(root / "barrier_architecture_candidates.tsv", sep="\t")
    primitives = pd.read_csv(root / "niche_primitives.tsv", sep="\t")
    increments = pd.read_csv(root / "topology_increment.tsv", sep="\t")
    required = {"primitive_id", "dataset_id", "patient_id", "estimate", "status"}
    missing = required - set(primitives.columns)
    if missing:
        raise ValueError(f"Stage 4 primitive table missing columns: {sorted(missing)}")
    return candidates, primitives, increments


def _context_implementations(primitives: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (primitive_id, dataset_id, platform), frame in primitives.groupby(
        ["primitive_id", "dataset_id", "platform"], dropna=False, sort=True
    ):
        values = pd.to_numeric(frame.estimate, errors="coerce")
        label, positive, negative = _sign_label(values)
        rows.append(
            {
                "primitive_id": primitive_id,
                "dataset_id": dataset_id,
                "context_label": CONTEXT_LABELS.get(str(dataset_id), "UNMAPPED_CONTEXT"),
                "platform": platform,
                "n_units": int(frame.unit_id.nunique()),
                "n_patients": int(frame.patient_id.nunique()),
                "n_estimable": int(frame.status.eq("ESTIMABLE").sum()),
                "median_estimate": float(values.median()),
                "q25_estimate": _safe_quantile(values, 0.25),
                "q75_estimate": _safe_quantile(values, 0.75),
                "positive_fraction": positive,
                "negative_fraction": negative,
                "direction_label": label,
                "response_blind": True,
                "claim_scope": "context_implementation_descriptive",
            }
        )
    return pd.DataFrame(rows)


def _context_heterogeneity(implementations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for primitive_id, frame in implementations.groupby("primitive_id", sort=True):
        medians = pd.to_numeric(frame.median_estimate, errors="coerce").dropna().to_numpy(float)
        label, positive, negative = _sign_label(medians)
        rows.append(
            {
                "primitive_id": primitive_id,
                "n_contexts": int(frame.dataset_id.nunique()),
                "n_platforms": int(frame.platform.nunique()),
                "context_median_min": float(np.min(medians)) if len(medians) else np.nan,
                "context_median_max": float(np.max(medians)) if len(medians) else np.nan,
                "context_median_range": float(np.ptp(medians)) if len(medians) else np.nan,
                "context_direction_label": label,
                "context_positive_fraction": positive,
                "context_negative_fraction": negative,
                "status": "DESCRIPTIVE_NO_INFERENCE",
                "claim_scope": "context_modulation_not_causal",
            }
        )
    return pd.DataFrame(rows)


def _hcc_residual_architectures(implementations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for primitive_id, frame in implementations.groupby("primitive_id", sort=True):
        hcc = frame[frame.dataset_id.isin(HCC_SPATIAL_DATASETS)]
        non_hcc = frame[~frame.dataset_id.isin(HCC_SPATIAL_DATASETS)]
        hcc_values = pd.to_numeric(hcc.median_estimate, errors="coerce").dropna()
        non_values = pd.to_numeric(non_hcc.median_estimate, errors="coerce").dropna()
        hcc_median = float(hcc_values.median()) if len(hcc_values) else np.nan
        non_median = float(non_values.median()) if len(non_values) else np.nan
        rows.append(
            {
                "primitive_id": primitive_id,
                "hcc_spatial_datasets": ";".join(sorted(hcc.dataset_id.astype(str).unique())),
                "non_hcc_spatial_datasets": ";".join(sorted(non_hcc.dataset_id.astype(str).unique())),
                "hcc_median_estimate": hcc_median,
                "non_hcc_median_estimate": non_median,
                "hcc_minus_non_hcc": hcc_median - non_median
                if np.isfinite(hcc_median) and np.isfinite(non_median)
                else np.nan,
                "status": "HCC_CONTEXT_DESCRIPTIVE"
                if np.isfinite(hcc_median) and np.isfinite(non_median)
                else "NOT_ESTIMABLE_MISSING_CONTEXT",
                "hcc_spatial_response_status": "NOT_IDENTIFIABLE_POST_ONLY_NO_RESPONSE_LINK",
                "hcc_clinical_anchor_dataset": "LAMBRECHT_HCC",
                "hcc_claim_scope": "HCC_context_not_pan_cancer",
                "response_blind_spatial_input": True,
            }
        )
    return pd.DataFrame(rows)


def _clinical_annotations(
    candidates: pd.DataFrame,
    effects: pd.DataFrame,
    ecology: pd.DataFrame,
) -> pd.DataFrame:
    environments = sorted(effects.environment.astype(str).unique())
    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.iterrows():
        primitive_id = str(candidate.primitive_id)
        for environment in environments:
            env_effect = effects[effects.environment.eq(environment)]
            primary = env_effect[env_effect.model.eq("abundance_primary")]
            composition = env_effect[env_effect.model.eq("composition")]
            joint = env_effect[env_effect.model.eq("joint_primary")]
            env_ecology = ecology[ecology.environment.eq(environment)]
            rows.append(
                {
                    "primitive_id": primitive_id,
                    "clinical_environment": environment,
                    "architecture_status_at_clinical_read": candidate.architecture_status,
                    "n_patients": int(primary.n_patients.iloc[0]) if len(primary) else np.nan,
                    "abundance_primary_auroc": float(primary.auroc.iloc[0]) if len(primary) else np.nan,
                    "composition_auroc": float(composition.auroc.iloc[0]) if len(composition) else np.nan,
                    "joint_primary_auroc": float(joint.auroc.iloc[0]) if len(joint) else np.nan,
                    "responder_compatible_feature_count": int(len(env_ecology)),
                    "spatial_response_bridge_status": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
                    "response_used_for_annotation_only": True,
                    "claim_scope": "clinical_anchor_annotation_not_spatial_association",
                }
            )
    return pd.DataFrame(rows)


def run_context_analysis(
    *,
    output_root: str | Path = "results/v7/context",
    stage4_root: str | Path = "results/v7/spatial_discovery",
    stage5_root: str | Path = "results/v7/clinical_anchor",
) -> dict[str, Any]:
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    spatial = Path(stage4_root)
    clinical = Path(stage5_root)
    candidates, primitives, increments = _read_stage4(spatial)
    effects = pd.read_csv(clinical / "environment_effects.tsv", sep="\t")
    ecology = pd.read_csv(clinical / "responder_compatible_ecology.tsv", sep="\t")

    # Candidate names and direction are inherited exactly from the response-
    # blind Stage 4 table.  Stage 6 cannot promote an unresolved primitive.
    locked = candidates[candidates.architecture_status.ne("UNRESOLVED_TOO_FEW_ESTIMABLE_UNITS")].copy()
    implementations = _context_implementations(
        primitives[primitives.primitive_id.isin(locked.primitive_id)]
    )
    heterogeneity = _context_heterogeneity(implementations)
    hcc = _hcc_residual_architectures(implementations)
    annotations = _clinical_annotations(locked, effects, ecology)

    annotations.to_csv(output / "shared_architecture_clinical_annotations.tsv", sep="\t", index=False)
    implementations.to_csv(output / "context_implementations.tsv", sep="\t", index=False)
    hcc.to_csv(output / "hcc_residual_architectures.tsv", sep="\t", index=False)
    heterogeneity.to_csv(output / "context_heterogeneity.tsv", sep="\t", index=False)

    # A compact input audit records that Stage 4 candidates were frozen before
    # this response-annotated stage.
    input_audit = pd.DataFrame(
        [
            {
                "input": "barrier_architecture_candidates.tsv",
                "sha256": sha256_file(spatial / "barrier_architecture_candidates.tsv"),
                "role": "locked_response_blind_candidate",
            },
            {
                "input": "niche_primitives.tsv",
                "sha256": sha256_file(spatial / "niche_primitives.tsv"),
                "role": "response_blind_context_measurement",
            },
            {
                "input": "environment_effects.tsv",
                "sha256": sha256_file(clinical / "environment_effects.tsv"),
                "role": "clinical_annotation_only",
            },
        ]
    )
    input_audit.to_csv(output / "context_input_audit.tsv", sep="\t", index=False)

    manifest = {
        "schema": "v7.stage6.context_analysis.run.v1",
        "status": "S6_COMPLETE_WITH_LIMITATIONS",
        "n_locked_candidates": int(len(locked)),
        "n_context_implementation_rows": int(len(implementations)),
        "n_hcc_rows": int(len(hcc)),
        "n_clinical_annotation_rows": int(len(annotations)),
        "spatial_datasets": sorted(primitives.dataset_id.astype(str).unique()),
        "clinical_environments": sorted(effects.environment.astype(str).unique()),
        "candidate_locked_before_response_read": True,
        "response_used_only_for_annotation": True,
        "spatial_response_bridge_status": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
        "hcc_spatial_response_status": "NOT_IDENTIFIABLE_POST_ONLY_NO_RESPONSE_LINK",
        "stage4_manifest_sha256": sha256_file(spatial / "STAGE4_RUN_MANIFEST.json"),
        "stage5_manifest_sha256": sha256_file(clinical / "STAGE5_RUN_MANIFEST.json"),
    }
    (output / "STAGE6_RUN_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    hcc_text = hcc[hcc.status.eq("HCC_CONTEXT_DESCRIPTIVE")]
    largest = hcc_text.iloc[int(np.argmax(np.abs(hcc_text.hcc_minus_non_hcc.to_numpy(float))))] if len(hcc_text) else None
    largest_line = (
        f"最大绝对 HCC—非 HCC 中位差出现在 `{largest.primitive_id}`，差值为 "
        f"{largest.hcc_minus_non_hcc:.3f}；这只是 context 描述，不是疗效或因果效应。"
        if largest is not None
        else "当前没有同时具备 HCC 与非 HCC 空间估计的原语。"
    )
    report = f"""# v7 Stage 6 context 与桥接报告

## 状态

**S6_COMPLETE_WITH_LIMITATIONS**。本阶段使用 Stage 4 已在 response 读取前锁定的 {len(locked)} 个空间候选，比较其在 {len(manifest['spatial_datasets'])} 个空间数据集中的实现，并把 Stage 5 的患者级临床基线作为注释层。没有把临床 response 重新用于空间候选选择。

## 结果 / 证据

- 形成 `{len(implementations)}` 行 dataset/platform context 实现表和 `{len(heterogeneity)}` 行 context 异质性表；候选的 carrier、平台和方向变化均保留。
- HCC 空间结果来自 GSE238264 的 post-only 样本；{largest_line}
- Stage 5 的 5 个临床环境被逐一注释，但空间—response 直接桥接为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。
- HCC 的临床锚来自 LAMBRECHT_HCC，与 GSE238264 不是同一患者或同一空间时间点，不能拼成纵向 HCC 空间疗效证据。

## 结论边界

当前支持“同一功能原语可由不同数据集/平台以不同量级实现”的 context 描述，并支持后续 HCC 与 PD1+X 设计优先检查髓系、CAF/血管和 T-cell 邻接。当前不支持泛癌共享疗效机制、HCC 空间疗效关联、PD1+X 修复或因果解释。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/context/run_stage6.py
```
"""
    (output / "SHARED_AND_HCC_CONTEXT_REPORT.md").write_text(report, encoding="utf-8")
    return {**manifest, "output_root": str(output.resolve())}


__all__ = ["run_context_analysis"]
