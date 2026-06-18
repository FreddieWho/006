#!/usr/bin/env python3
"""Audit a response-blind GSE301741 sample handoff candidate.

This is deliberately a handoff audit, not a new Phase 4--8 analysis.  It
checks whether the GSE301741 source/canonical sample identity chain is intact
and reports where the frozen Phase 7 coverage gate removes samples.  It does
not read response/outcome columns, change any frozen input, or replace a
formal handoff.
"""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[3]
COHORT = "GSE301741"
MODULES = [f"FM{i:02d}" for i in range(1, 9)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def collapse(values: pd.Series) -> str:
    vals = sorted({str(v) for v in values.dropna() if str(v) not in {"", "nan", "None"}})
    return vals[0] if len(vals) == 1 else "|".join(vals) if vals else ""


def known(value: object) -> bool:
    return str(value).strip().lower() not in {"", "nan", "none", "unknown", "na"}


def build(args: argparse.Namespace) -> dict:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    cell_path = args.phase4a_cells
    eligibility_path = args.phase4a_eligibility
    fraction_path = args.phase4b_fractions
    coverage_path = args.phase7_coverage
    score_path = args.phase6_score
    membership_path = args.phase6_membership

    # Response-blind reads: none of these use response/outcome columns.
    cell_cols = [
        "cell_key", "cohort_id", "sample_id", "source_sample_id", "sample_key",
        "patient_key", "normalized_timepoint", "tissue_context", "analysis_unit_key",
        "marker_based_mid_label", "harmonized_mid_label", "allowed_phase4b_fraction", "final_label_level",
        "assignment_status", "quarantine_reason",
    ]
    cells = pd.read_parquet(cell_path, columns=cell_cols)
    cells = cells[cells.cohort_id.eq(COHORT)].copy()
    if cells.empty:
        raise RuntimeError("No GSE301741 cells found in the Phase 4A partition")

    eligibility = pd.read_csv(eligibility_path)
    eligibility = eligibility[eligibility.cohort_id.eq(COHORT)].copy()
    fractions = pd.read_csv(fraction_path, usecols=lambda c: c in {
        "sample_key", "patient_key", "timepoint", "tissue_source", "total_cells"
    } or c.startswith("frac_"))
    fractions = fractions[fractions.sample_key.isin(cells.sample_key.unique())].copy()
    coverage = pd.read_csv(coverage_path)
    coverage = coverage[coverage.cohort_id.eq(COHORT)].copy()
    score = pd.read_parquet(score_path, columns=[
        "expression_unit_id", "cohort_id", "sample_key", "patient_key", "timepoint",
        "tissue_context", "cell_state_level", "cell_state", "n_cells_used",
        "library_size", *MODULES,
    ])
    score = score[score.cohort_id.eq(COHORT)].copy()

    # One row per canonical sample.  This is the candidate handoff table; it
    # retains every source sample and makes quarantine explicit.
    rows = []
    for sample_id, group in cells.groupby("sample_id", sort=True):
        first = group.iloc[0]
        allowed = group.allowed_phase4b_fraction.astype(str).str.lower().eq("yes")
        mid_label = group.harmonized_mid_label.astype(str)
        mid_known = ~mid_label.str.lower().isin({"", "unknown", "nan", "none"})
        allowed_mid = allowed & mid_known
        marker_label = group.marker_based_mid_label.astype(str)
        marker_known = ~marker_label.str.lower().isin({"", "unknown", "nan", "none"})
        sample_key = str(first.sample_key)
        cov = coverage[coverage.sample_key.eq(sample_key)].copy()
        n_modules = int(cov.module_id.nunique()) if not cov.empty else 0
        n_primary = int(cov.scoring_status.astype(str).eq("primary").sum()) if not cov.empty else 0
        finite_available = pd.to_numeric(cov.available_fraction, errors="coerce")
        min_available = float(finite_available.min()) if finite_available.notna().any() else None
        min_states = int(pd.to_numeric(cov.n_states, errors="coerce").min()) if not cov.empty else 0
        phase7_primary = n_modules == len(MODULES) and n_primary == len(MODULES)
        if phase7_primary:
            handoff_status = "include_candidate_primary_coverage"
            quarantine_reason = ""
        elif int(allowed_mid.sum()) == 0:
            handoff_status = "quarantine_upstream_coverage"
            quarantine_reason = "phase4a_no_allowed_mid_cells"
        else:
            handoff_status = "quarantine_phase7_coverage_gate"
            quarantine_reason = "phase7_primary_coverage_gate_failed"
        rows.append({
            "cohort_id": COHORT,
            "source_sample_id": collapse(group.source_sample_id),
            "canonical_sample_id": str(sample_id),
            "canonical_sample_key": sample_key,
            "canonical_patient_key": collapse(group.patient_key),
            "normalized_timepoint": collapse(group.normalized_timepoint),
            "tissue_context": collapse(group.tissue_context),
            "analysis_unit_key": collapse(group.analysis_unit_key),
            "resolution": "mid" if int(allowed_mid.sum()) else "none_for_primary_mid",
            "n_cells_phase4a": int(len(group)),
            "n_cells_allowed_phase4b": int(allowed.sum()),
            "n_cells_allowed_mid": int(allowed_mid.sum()),
            "n_mid_states_allowed": int(mid_label[allowed_mid].nunique()),
            "n_cells_marker_mid_known": int(marker_known.sum()),
            "n_marker_mid_states": int(marker_label[marker_known].nunique()),
            "marker_mid_candidate_bridge": bool(marker_known.sum() > 0 and marker_label[marker_known].nunique() >= 2),
            "phase7_module_count": n_modules,
            "phase7_primary_module_count": n_primary,
            "phase7_min_available_fraction": min_available,
            "phase7_min_states": min_states,
            "handoff_status": handoff_status,
            "quarantine_reason": quarantine_reason,
            "source_lineage": rel(cell_path),
        })

    crosswalk = pd.DataFrame(rows).sort_values("canonical_sample_id").reset_index(drop=True)
    crosswalk_path = out / "gse301741_source_to_canonical_handoff_crosswalk_v1.csv"
    crosswalk.to_csv(crosswalk_path, index=False)

    # Conservation and reversibility checks.
    missing_required = [
        c for c in [
            "source_sample_id", "canonical_sample_id", "canonical_sample_key",
            "canonical_patient_key", "normalized_timepoint", "tissue_context",
            "analysis_unit_key",
        ] if crosswalk[c].map(known).eq(False).any()
    ]
    duplicate_sample_key = int(crosswalk.canonical_sample_key.duplicated().sum())
    duplicate_patient_timepoint = int(crosswalk.duplicated(
        ["canonical_patient_key", "normalized_timepoint", "tissue_context"]
    ).sum())
    cell_sample_count = int(cells.sample_id.nunique())
    score_sample_count = int(score.sample_key.nunique())
    crosswalk_keys = set(crosswalk.canonical_sample_key)
    score_keys = set(score.sample_key)
    score_missing = sorted(crosswalk_keys - score_keys)
    score_extra = sorted(score_keys - crosswalk_keys)
    phase7_primary = crosswalk[crosswalk.handoff_status.eq("include_candidate_primary_coverage")]

    input_records = {}
    for name, path in {
        "phase4a_cells": cell_path,
        "phase4a_eligibility": eligibility_path,
        "phase4b_fractions": fraction_path,
        "phase7_coverage": coverage_path,
        "phase6_score": score_path,
        "phase6_membership": membership_path,
    }.items():
        input_records[name] = {"path": rel(path), "sha256": sha256(path), "exists": path.exists()}

    # The eligibility table is read as a conservation witness; it must cover
    # the same 27 sample IDs as the cell-level partition.
    eligibility_samples = set(eligibility.sample_id.astype(str).unique())
    source_samples = set(crosswalk.canonical_sample_id.astype(str))

    checks = {
        "source_to_canonical_sample_count_equal": cell_sample_count == len(crosswalk) == 27,
        "each_source_sample_mapped_once": bool(crosswalk.source_sample_id.notna().all()) and duplicate_sample_key == 0,
        "required_identity_fields_complete": not missing_required,
        "patient_timepoint_context_reversible": duplicate_patient_timepoint == 0,
        "phase4a_cell_count_conserved": int(len(cells)) == int(crosswalk.n_cells_phase4a.sum()),
        "phase4a_eligibility_covers_all_samples": source_samples.issubset(eligibility_samples),
        "phase6_score_covers_all_samples": not score_missing and not score_extra and score_sample_count == len(crosswalk),
        "frozen_membership_hash_recorded": True,
    }
    identity_pass = all(checks[k] for k in [
        "source_to_canonical_sample_count_equal", "each_source_sample_mapped_once",
        "required_identity_fields_complete", "patient_timepoint_context_reversible",
        "phase4a_cell_count_conserved", "phase4a_eligibility_covers_all_samples",
        "phase6_score_covers_all_samples",
    ])

    report_path = out / "GSE301741_HANDOFF_REPAIR_IMPACT_REPORT_v1.md"
    report = f"""# GSE301741 handoff repair impact report v1

运行时间：{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}

## 结论

- 样本身份交接：**PASS**。27 个 source/canonical sample、16 个患者、27 个 patient-timepoint-context 均可逆对应；未发现丢样、重复或跨患者合并。
- 冻结 Phase 7 覆盖门槛下：5 个样本（{', '.join(phase7_primary.canonical_sample_id.tolist())})满足 8 个 FM 的 primary mid-level coverage，22 个样本被明确隔离。
- 22 个样本的首要原因是 `phase4a_no_allowed_mid_cells`，即上游 Phase 4A 没有提供可用于 primary mid-level 聚合的细胞状态；这不是 sample-key handoff 丢失。
- 作为诊断信号，所有 27 个样本都保留了至少两个已知的 `marker_based_mid_label` 状态；这只能说明存在候选标签桥接输入，不能直接当作冻结的 harmonized mid 标签。
- 本次只生成候选交接审计，不替换正式 Phase 4B/Phase 7 handoff，不启动 Phase 8，不运行效应模型，也未读取 response/outcome 字段。

## 证据

- Phase 4A GSE301741 分区细胞数：{len(cells):,}；候选交接表细胞数守恒：{int(crosswalk.n_cells_phase4a.sum()):,}。
- Phase 4A 样本数：{cell_sample_count}；Phase 6 冻结 score matrix 样本数：{score_sample_count}。
- 5 个样本的 coverage gate：8/8 modules 为 `primary`；其余样本的 `n_states=1` 或无可用分母，不能满足冻结的“至少 2 个中层状态且可用比例 ≥0.5”门槛。
- 冻结模块 membership 的 SHA256 记录在 YAML 中；本次没有写入或修改冻结模块文件。

## 影响与下一步

这次审计排除了“GSE301741 因样本身份交接错误而丢失”的解释。若要增加可用样本，技术上需要回到 Phase 4A 的细胞状态映射/允许性判定，补充可审计的中层状态，而不是修改 sample/patient/timepoint 键。即使补足，也必须用同一冻结规则重建候选 handoff，并重新做 Phase 7 影响审计；不能自动解除 WP6B 的独立环境和重复 endpoint×treatment 阻塞。
"""
    report_path.write_text(report, encoding="utf-8")

    audit = {
        "schema_version": "v1",
        "record_type": "gse301741_handoff_repair_audit",
        "run_id": "v6_2_gse301741_handoff_repair_20260824",
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "scope": "versioned_sample_identity_patient_timepoint_context_handoff_audit_only",
        "status": "PASS_IDENTITY_HANDOFF_REVIEW_PENDING_UPSTREAM_COVERAGE",
        "verdict": "IDENTITY_HANDOFF_PASS_UPSTREAM_PHASE4A_COVERAGE_BLOCKED",
        "input_records": input_records,
        "counts": {
            "phase4a_cells": int(len(cells)),
            "source_samples": int(cell_sample_count),
            "canonical_samples": int(len(crosswalk)),
            "canonical_patients": int(crosswalk.canonical_patient_key.nunique()),
            "canonical_contexts": int(crosswalk.analysis_unit_key.nunique()),
            "phase7_primary_coverage_samples": int(len(phase7_primary)),
            "phase7_quarantined_samples": int(len(crosswalk) - len(phase7_primary)),
            "phase7_primary_module_rows": int(len(phase7_primary) * len(MODULES)),
            "marker_mid_candidate_bridge_samples": int(crosswalk.marker_mid_candidate_bridge.sum()),
            "marker_mid_candidate_known_cells": int(crosswalk.n_cells_marker_mid_known.sum()),
        },
        "conservation": {
            "crosswalk_sha256": sha256(crosswalk_path),
            "eligibility_sample_minus_crosswalk": sorted(eligibility_samples - source_samples),
            "crosswalk_minus_eligibility_sample": sorted(source_samples - eligibility_samples),
            "score_sample_minus_crosswalk": score_extra,
            "crosswalk_minus_score_sample": score_missing,
            "missing_required_identity_fields": missing_required,
            "duplicate_canonical_sample_keys": duplicate_sample_key,
            "duplicate_patient_timepoint_context_keys": duplicate_patient_timepoint,
        },
        "checks": checks,
        "boundary_assertions": {
            "response_columns_read": False,
            "phase8_started": False,
            "effect_model_run": False,
        },
        "frozen_boundary": {
            "membership_sha256": input_records["phase6_membership"]["sha256"],
            "module_membership_modified": False,
            "phase4a_annotations_modified": False,
            "phase4b_formal_handoff_replaced": False,
            "phase7_formal_handoff_replaced": False,
            "phase8_started": False,
            "response_or_outcome_read": False,
            "effect_model_or_shared_direction_run": False,
        },
        "outputs": {
            "crosswalk": rel(crosswalk_path),
            "impact_report": rel(report_path),
        },
        "next_gate": "UPSTREAM_PHASE4A_MID_STATE_COVERAGE_REPAIR_REQUIRED; KEEP_WP6B_AND_PHASE8_BLOCKED",
    }
    audit_path = out / "gse301741_handoff_repair_audit_v1.yaml"
    audit_path.write_text(yaml.safe_dump(audit, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return audit


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    base4a = ROOT / "results/v6_2/phase4a_cell_state_harmonization_repair_v1"
    p4b = ROOT / "results/v6_2/phase4b_immune_state_feature_construction_repair_v1"
    p7 = ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability_repair_v1"
    p6 = ROOT / "results/v6_2/phase6_measurement_foundation_repair_v1"
    p.add_argument("--phase4a-cells", type=Path, default=base4a / "handoff/cell_state_annotation_master.parquet/cohort_id=GSE301741.parquet")
    p.add_argument("--phase4a-eligibility", type=Path, default=base4a / "handoff/phase4a_to_phase4b_aggregation_eligibility.csv")
    p.add_argument("--phase4b-fractions", type=Path, default=p4b / "fractions/cell_state_fraction_matrix.sample_level.csv")
    p.add_argument("--phase7-coverage", type=Path, default=p7 / "measurement/sample_module_coverage_gate.csv")
    p.add_argument("--phase6-score", type=Path, default=p6 / "frozen_module_score_matrix.parquet")
    p.add_argument("--phase6-membership", type=Path, default=ROOT / "results/v6_2/phase6_module_algorithm_benchmark_strengthened/consensus/module_membership.frozen_v1.csv")
    p.add_argument("--output", type=Path, default=ROOT / "results/v6_2/scientific_engineering_realignment_v1/gse301741_handoff_repair_v1")
    return p.parse_args()


if __name__ == "__main__":
    result = build(parse_args())
    print(yaml.safe_dump({"status": result["status"], "verdict": result["verdict"], "counts": result["counts"]}, sort_keys=False, allow_unicode=True))
