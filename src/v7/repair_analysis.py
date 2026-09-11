"""Stage 7 molecular/cell-state repair baseline with a fail-closed spatial branch.

The molecular branch uses audited TASK01/TASK02 patient time points.  It
computes paired T1/T2 minus T0 displacements and response-stratified
descriptions, while keeping response permutation and pairing permutation nulls
separate.  Spatial rewiring is not inferred from these molecular changes.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .clinical_baseline import _read_abundance_scores, _read_response_table
from .spatial_control import sha256_file, stable_hash


TARGET_ENVIRONMENTS = ("TASK01", "TASK02")
TARGET_TIMEPOINTS = ("T1", "T2")


def _paired_rows(scores: pd.DataFrame, response: pd.DataFrame) -> pd.DataFrame:
    score = scores[
        scores.analysis_dataset.isin(TARGET_ENVIRONMENTS)
        & scores.timepoint.isin(("T0", *TARGET_TIMEPOINTS))
    ].copy()
    response_cols = [
        "analysis_dataset",
        "patient_key",
        "response_binary",
        "treatment_values",
        "endpoint_values",
        "complete_for_response",
        "conflict_flag",
    ]
    score = score.merge(response[response_cols], on=["analysis_dataset", "patient_key"], how="inner")
    rows: list[pd.DataFrame] = []
    for environment in TARGET_ENVIRONMENTS:
        frame = score[score.analysis_dataset.eq(environment)]
        pre = frame[frame.timepoint.eq("T0")]
        for target in TARGET_TIMEPOINTS:
            post = frame[frame.timepoint.eq(target)]
            if pre.empty or post.empty:
                continue
            merged = pre.merge(
                post,
                on=["analysis_dataset", "patient_key", "feature_id"],
                suffixes=("_pre", "_post"),
            )
            if merged.empty:
                continue
            merged["delta"] = merged.score_standardized_post - merged.score_standardized_pre
            merged["target_timepoint"] = target
            rows.append(
                merged[
                    [
                        "analysis_dataset",
                        "patient_key",
                        "feature_id",
                        "target_timepoint",
                        "response_binary_pre",
                        "treatment_values_pre",
                        "endpoint_values_pre",
                        "complete_for_response_pre",
                        "conflict_flag_pre",
                        "coverage_pre",
                        "coverage_post",
                        "n_cells_used_pre",
                        "n_cells_used_post",
                        "score_standardized_pre",
                        "score_standardized_post",
                        "delta",
                    ]
                ].rename(columns={"analysis_dataset": "environment", "response_binary_pre": "response_binary"})
            )
    if not rows:
        return pd.DataFrame()
    result = pd.concat(rows, ignore_index=True)
    result["status"] = np.where(
        result.groupby(["environment", "patient_key", "target_timepoint"]).delta.transform("count") > 0,
        "PAIRED_T0_TO_TARGET",
        "NOT_ESTIMABLE",
    )
    return result


def _repair_compatibility(paired: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (environment, target, feature), frame in paired.groupby(
        ["environment", "target_timepoint", "feature_id"], sort=True
    ):
        responders = frame.loc[frame.response_binary.eq(1), "delta"].dropna()
        nonresponders = frame.loc[frame.response_binary.eq(0), "delta"].dropna()
        r_median = float(responders.median()) if len(responders) else np.nan
        n_median = float(nonresponders.median()) if len(nonresponders) else np.nan
        contrast = r_median - n_median if np.isfinite(r_median) and np.isfinite(n_median) else np.nan
        if len(responders) < 3 or len(nonresponders) < 3:
            status = "NOT_ESTIMABLE_TOO_FEW_PATIENTS"
        else:
            status = "DESCRIPTIVE_RESPONSE_STRATIFIED_DISPLACEMENT"
        rows.append(
            {
                "environment": environment,
                "target_timepoint": target,
                "feature_id": feature,
                "n_responder": int(len(responders)),
                "n_non_responder": int(len(nonresponders)),
                "responder_median_delta": r_median,
                "non_responder_median_delta": n_median,
                "responder_minus_non_responder": contrast,
                "status": status,
                "repair_direction": "NOT_ASSIGNED_WITHOUT_MECHANISTIC_TARGET_DIRECTION",
                "causal_claim": "forbidden",
            }
        )
    return pd.DataFrame(rows)


def _contrast(frame: pd.DataFrame, labels: pd.Series) -> float:
    values = pd.to_numeric(frame.delta, errors="coerce")
    r = values[labels.to_numpy() == 1].dropna()
    n = values[labels.to_numpy() == 0].dropna()
    if len(r) < 3 or len(n) < 3:
        return np.nan
    return float(np.median(r) - np.median(n))


def _permute_patient_vectors(frame: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Break pairing using one donor mapping for the complete post vector.

    Called within one environment/timepoint. Donor missingness travels with
    the vector; do not independently shuffle features or impute missing scores.
    """
    keys = ["patient_key", "feature_id"]
    if frame.duplicated(keys).any():
        raise ValueError("duplicate patient-feature rows in paired permutation")
    patients = np.sort(frame.patient_key.unique())
    donors = dict(zip(patients, rng.permutation(patients), strict=True))
    post = frame.set_index(keys).score_standardized_post
    donor_index = pd.MultiIndex.from_arrays([
        frame.patient_key.map(donors), frame.feature_id,
    ], names=keys)
    shuffled = frame.copy()
    shuffled["delta"] = post.reindex(donor_index).to_numpy(float) - frame.score_standardized_pre.to_numpy(float)
    return shuffled


def _negative_controls(paired: pd.DataFrame, permutations: int, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (environment, target), frame in paired.groupby(["environment", "target_timepoint"], sort=True):
        patients = frame[["patient_key", "response_binary"]].drop_duplicates("patient_key")
        if patients.response_binary.nunique() < 2 or len(patients) < 6:
            rows.append(
                {
                    "environment": environment,
                    "target_timepoint": target,
                    "control": "response_permutation",
                    "status": "NOT_TESTABLE",
                    "n_permutations": permutations,
                }
            )
            continue
        observed = []
        for _, feature_frame in frame.groupby("feature_id", sort=True):
            labels = feature_frame.response_binary
            value = _contrast(feature_frame, labels)
            if np.isfinite(value):
                observed.append(abs(value))
        rng = np.random.default_rng(seed ^ int(stable_hash({"environment": environment, "target": target})[:8], 16))
        response_null: list[float] = []
        pairing_null: list[float] = []
        response_values = patients.response_binary.to_numpy(int)
        patient_order = patients.patient_key.to_numpy(str)
        for _ in range(permutations):
            shuffled_values = rng.permutation(response_values)
            mapping = dict(zip(patient_order, shuffled_values, strict=True))
            shuffled_labels = frame.patient_key.map(mapping).astype(int)
            per_feature = []
            for _, feature_frame in frame.groupby("feature_id", sort=True):
                local_labels = shuffled_labels.loc[feature_frame.index]
                value = _contrast(feature_frame, local_labels)
                if np.isfinite(value):
                    per_feature.append(abs(value))
            if per_feature:
                response_null.append(float(np.median(per_feature)))

            # Break the T0-to-target pairing while preserving each feature's
            # target distribution and the patient-level response labels.
            shuffled = _permute_patient_vectors(frame, rng)
            per_feature = []
            for _, feature_frame in shuffled.groupby("feature_id", sort=True):
                value = _contrast(feature_frame, feature_frame.response_binary)
                if np.isfinite(value):
                    per_feature.append(abs(value))
            if per_feature:
                pairing_null.append(float(np.median(per_feature)))

        observed_value = float(np.median(observed)) if observed else np.nan
        for control, values in (("response_permutation", response_null), ("pairing_permutation", pairing_null)):
            rows.append(
                {
                    "environment": environment,
                    "target_timepoint": target,
                    "control": control,
                    "status": "ESTIMABLE" if values else "NOT_TESTABLE",
                    "n_permutations": permutations,
                    "n_successful": len(values),
                    "observed_median_abs_contrast": observed_value,
                    "null_median": float(np.median(values)) if values else np.nan,
                    "null_q025": float(np.quantile(values, 0.025)) if values else np.nan,
                    "null_q975": float(np.quantile(values, 0.975)) if values else np.nan,
                }
            )
    return pd.DataFrame(rows)


def run_repair_analysis(
    *,
    output_root: str | Path = "results/v7/repair",
    response_table_path: str | Path = "results/v7/registry/treatment_response_completeness.tsv",
    patient_scores_path: str | Path = "results/v7/ontology/scores/patient_timepoint_scores.parquet",
    config_path: str | Path = "config/v7/stage7.yaml",
    permutations: int = 50,
    seed: int = 20260911,
) -> dict[str, Any]:
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    response = _read_response_table(response_table_path)
    scores = _read_abundance_scores(patient_scores_path)
    paired = _paired_rows(scores, response)
    compatibility = _repair_compatibility(paired)
    controls = _negative_controls(paired, permutations, seed)

    paired.to_csv(output / "paired_displacements.tsv", sep="\t", index=False)
    compatibility.to_csv(output / "repair_compatibility.tsv", sep="\t", index=False)
    controls.to_csv(output / "negative_combination_controls.tsv", sep="\t", index=False)
    spatial = pd.DataFrame(
        [
            {
                "branch": "observed_spatial_rewiring",
                "status": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
                "reason": "no_same_patient_pre_on_post_spatial_response_asset",
                "claim_allowed": "none",
            }
        ]
    )
    surrogate = pd.DataFrame(
        [
            {
                "branch": "surrogate_inferred_architecture_rewiring",
                "status": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
                "reason": "Stage4B_surrogate_not_validated_on_independent_paired_patients",
                "claim_allowed": "none",
            }
        ]
    )
    spatial.to_csv(output / "observed_spatial_rewiring.tsv", sep="\t", index=False)
    surrogate.to_csv(output / "surrogate_inferred_architecture_rewiring.tsv", sep="\t", index=False)
    input_audit = pd.DataFrame(
        [
            {"input": "treatment_response_completeness.tsv", "sha256": sha256_file(response_table_path), "role": "response_label_audit"},
            {"input": "patient_timepoint_scores.parquet", "sha256": sha256_file(patient_scores_path), "role": "paired_molecular_scores"},
        ]
    )
    input_audit.to_csv(output / "repair_input_audit.tsv", sep="\t", index=False)

    manifest = {
        "schema": "v7.stage7.repair_analysis.run.v1",
        "status": "S7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED",
        "molecular_branch": "COMPLETE_WITH_LIMITATIONS",
        "spatial_branch": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
        "n_paired_rows": int(len(paired)),
        "n_paired_patients": int(paired.patient_key.nunique()) if not paired.empty else 0,
        "n_paired_patients_by_environment": paired.groupby("environment").patient_key.nunique().to_dict()
        if not paired.empty
        else {},
        "target_timepoints": sorted(paired.target_timepoint.unique()) if not paired.empty else [],
        "response_permutations": permutations,
        "response_read": True,
        "spatial_response_read": False,
        "response_table_sha256": sha256_file(response_table_path),
        "patient_scores_sha256": sha256_file(patient_scores_path),
        "config_sha256": sha256_file(config_path),
    }
    (output / "STAGE7_RUN_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    response_features = int(compatibility.feature_id.nunique()) if not compatibility.empty else 0
    report = f"""# v7 Stage 7 PD1+X repair 基线报告

## 状态

**S7_MOLECULAR_COMPLETE_SPATIAL_BLOCKED**。本轮在 TASK01/TASK02 中完成 {manifest['n_paired_rows']} 条患者内分子位移记录，覆盖 {manifest['n_paired_patients']} 位患者和 T1/T2 两个目标时间点；空间重排分支保持不可识别。

## 结果 / 证据

- 以患者为单位计算 T1−T0、T2−T0 的 39 个 Stage 2 特征位移，并按 response 分层描述；没有把分子位移改写为空间屏障重排。
- `{response_features}` 个特征形成 response-stratified displacement 表；`repair_direction` 保持未指定，因为没有给每个分子特征预设 X-class 机制方向。
- 运行 response permutation 和 pairing permutation 负对照；结果用于判断配对/标签结构是否足以解释观察到的差异。
- `observed_spatial_rewiring.tsv` 与 `surrogate_inferred_architecture_rewiring.tsv` 均明确为 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

## 结论边界

当前可以继续研究 PD-1 单药与 PD1+lenvatinib 的患者内分子/细胞状态位移，不能声称 X 的因果增量，也不能声称空间屏障被修复。要升级 spatial repair，需要同一患者的纵向空间数据，或已经在独立配对患者上验证通过的 Stage 4B surrogate。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/repair/run_stage7.py --response-permutations {permutations} --config config/v7/stage7.yaml
```
"""
    (output / "REPAIR_REPORT.md").write_text(report, encoding="utf-8")
    return {**manifest, "output_root": str(output.resolve())}


__all__ = ["run_repair_analysis"]
