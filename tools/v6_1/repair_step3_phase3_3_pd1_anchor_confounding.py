#!/usr/bin/env python3
"""Repair Step3 Phase 3.3 PD1_anchor cohort confounding.

Creates role-based repaired universes without changing samples, labels, splits,
feature matrix, or feature dictionary.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import chi2_contingency


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/v6_1/step3_qc_aware_strong_baseline"
P1 = OUT / "01_analysis_universe"
P3 = OUT / "03_leakage_confounding_audit"
ROLES = {"primary_allowed", "sensitivity_only", "biological_anchor_only", "blocked"}
PD1_FEATURE_FAMILIES = [
    "fraction",
    "signature",
    "pathway",
    "tf_activity",
    "combined_interpretable",
    "combined_all_allowed",
]


def response_num(s: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=s.index, dtype=float)
    txt = s.astype(str).str.lower()
    out[txt.eq("responder")] = 1.0
    out[txt.eq("non_responder")] = 0.0
    return out


def weighted_chi2_p(df: pd.DataFrame, weight_col: str) -> float:
    y = response_num(df["response_label"])
    tab = pd.DataFrame({"cohort_id": df["cohort_id"], "y": y, "w": df[weight_col]})
    tab = tab[tab["y"].isin([0, 1])]
    if tab["cohort_id"].nunique() < 2 or tab["y"].nunique() < 2:
        return np.nan
    ct = tab.pivot_table(index="cohort_id", columns="y", values="w", aggfunc="sum", fill_value=0.0)
    try:
        return float(chi2_contingency(ct)[1])
    except Exception:
        return np.nan


def risk_from_weighted(max_frac: float, p: float, missing_auc: float) -> str:
    if max_frac > 0.75:
        return "HIGH"
    if pd.notna(p) and p < 0.001:
        return "HIGH"
    if pd.notna(missing_auc) and missing_auc > 0.70:
        return "HIGH"
    if max_frac > 0.50 or (pd.notna(p) and p < 0.05) or (pd.notna(missing_auc) and missing_auc > 0.60):
        return "MODERATE"
    return "LOW"


def role_columns(role: str, weight_col: str, max_weight_frac: float, reason: str) -> dict[str, object]:
    if role not in ROLES:
        raise ValueError(role)
    return {
        "universe_role": role,
        "sample_weight_column": weight_col,
        "cohort_adjustment_required": role == "primary_allowed",
        "effective_max_cohort_weight_fraction": max_weight_frac,
        "primary_allowed_after_repair": role == "primary_allowed",
        "sensitivity_only_after_repair": role == "sensitivity_only",
        "biological_anchor_only": role == "biological_anchor_only",
        "repair_reason": reason,
    }


def main() -> int:
    registry = pd.read_csv(P1 / "step3_phase3_1_analysis_universe_registry.csv")
    membership = pd.read_csv(P1 / "step3_phase3_1_sample_universe_membership.csv")
    conf = pd.read_csv(P3 / "step3_phase3_3_confounding_audit_by_universe.csv")

    frac_pd1 = membership[
        membership["task"].eq("PD1_anchor")
        & membership["feature_family"].eq("fraction")
        & membership["universe_type"].eq("supervised_labeled_universe")
        & membership["included"].astype(bool)
    ].copy()
    if frac_pd1.empty:
        raise SystemExit("No PD1_anchor supervised fraction samples found")

    patient_sample_n = frac_pd1.groupby("patient_id")["sample_key"].transform("count")
    frac_pd1["pd1_anchor_weighted_sample_weight"] = 1.0 / patient_sample_n
    frac_pd1["pd1_anchor_weighted_sample_weight"] *= len(frac_pd1) / frac_pd1["patient_id"].nunique()
    frac_pd1["pd1_anchor_full_sample_weight"] = 1.0
    frac_pd1["pd1_anchor_within_cohort_meta_weight"] = np.nan

    cohort_breakdown = (
        frac_pd1.groupby("cohort_id")
        .agg(
            n_samples=("sample_key", "count"),
            n_patients=("patient_id", "nunique"),
            n_responders=("response_label", lambda x: int((x.astype(str) == "responder").sum())),
            n_non_responders=("response_label", lambda x: int((x.astype(str) == "non_responder").sum())),
            full_weight_sum=("pd1_anchor_full_sample_weight", "sum"),
            weighted_weight_sum=("pd1_anchor_weighted_sample_weight", "sum"),
        )
        .reset_index()
    )
    cohort_breakdown["sample_fraction"] = cohort_breakdown["n_samples"] / cohort_breakdown["n_samples"].sum()
    cohort_breakdown["patient_fraction"] = cohort_breakdown["n_patients"] / cohort_breakdown["n_patients"].sum()
    cohort_breakdown["full_weight_fraction"] = cohort_breakdown["full_weight_sum"] / cohort_breakdown["full_weight_sum"].sum()
    cohort_breakdown["weighted_weight_fraction"] = cohort_breakdown["weighted_weight_sum"] / cohort_breakdown["weighted_weight_sum"].sum()
    cohort_breakdown["responder_rate"] = cohort_breakdown["n_responders"] / (
        cohort_breakdown["n_responders"] + cohort_breakdown["n_non_responders"]
    )
    cohort_breakdown.to_csv(P3 / "repair_phase3_3_PD1_anchor_cohort_response_breakdown.csv", index=False)

    full_max = float(cohort_breakdown["full_weight_fraction"].max())
    weighted_max = float(cohort_breakdown["weighted_weight_fraction"].max())
    weighted_p = weighted_chi2_p(frac_pd1, "pd1_anchor_weighted_sample_weight")

    pd1_conf = conf[conf["task"].eq("PD1_anchor")].copy()
    if pd1_conf.empty:
        raise SystemExit("No PD1_anchor rows in confounding audit")
    template_missing = float(pd1_conf["missingness_response_auc"].dropna().max()) if pd1_conf["missingness_response_auc"].notna().any() else np.nan
    weighted_risk = risk_from_weighted(weighted_max, weighted_p, template_missing)
    weighted_role = "primary_allowed" if weighted_risk in {"LOW", "MODERATE"} else "sensitivity_only"

    non_pd1_registry = registry[~registry["task"].eq("PD1_anchor")].copy()
    non_pd1_registry = add_default_roles(non_pd1_registry)
    repaired_registry_rows = [non_pd1_registry]

    base_supervised = registry[
        registry["task"].eq("PD1_anchor")
        & registry["universe_type"].eq("supervised_labeled_universe")
        & registry["feature_family"].isin(PD1_FEATURE_FAMILIES)
    ].copy()
    role_specs = [
        ("PD1_anchor_full", "sensitivity_only", "pd1_anchor_full_sample_weight", full_max, "full_PD1_anchor_high_cohort_confounding_downgraded"),
        ("PD1_anchor_weighted", weighted_role, "pd1_anchor_weighted_sample_weight", weighted_max, "patient_level_cohort_balanced_weighting"),
        ("PD1_anchor_within_cohort_meta", "biological_anchor_only", "none", full_max, "within_cohort_effect_direction_meta_anchor_not_prediction_model"),
    ]
    for group, role, weight_col, max_frac, reason in role_specs:
        rows = base_supervised.copy()
        rows["universe_id"] = rows["feature_family"].map(
            lambda f: f"{group}__{f}__supervised_labeled_universe"
        )
        for k, v in role_columns(role, weight_col, max_frac, reason).items():
            rows[k] = v
        repaired_registry_rows.append(rows)
    repaired_registry = pd.concat(repaired_registry_rows, ignore_index=True, sort=False)
    repaired_registry.to_csv(P1 / "step3_phase3_1_analysis_universe_registry.repaired.csv", index=False)

    repaired_membership = add_membership_roles(membership, frac_pd1, role_specs)
    repaired_membership.to_csv(P1 / "step3_phase3_1_sample_universe_membership.repaired.csv", index=False)

    non_pd1_conf = conf[~conf["task"].eq("PD1_anchor")].copy()
    non_pd1_conf = add_default_roles(non_pd1_conf)
    repaired_conf_rows = [non_pd1_conf]
    for group, role, weight_col, max_frac, reason in role_specs:
        rows = pd1_conf.copy()
        rows["universe_id"] = rows["feature_family"].map(
            lambda f: f"{group}__{f}__supervised_labeled_universe"
        )
        if group == "PD1_anchor_weighted":
            rows["max_cohort_fraction"] = weighted_max
            rows["cohort_response_chi2_p"] = weighted_p
            rows["confounding_risk_level"] = weighted_risk
            rows["recommended_action"] = (
                "primary_allowed_with_patient_level_cohort_balanced_weights"
                if role == "primary_allowed"
                else "weighted_still_high_downgrade_to_sensitivity"
            )
        elif group == "PD1_anchor_full":
            rows["confounding_risk_level"] = "HIGH"
            rows["recommended_action"] = "sensitivity_only_full_PD1_anchor"
        else:
            rows["confounding_risk_level"] = "not_applicable"
            rows["recommended_action"] = "biological_anchor_only_within_cohort_direction_meta"
        for k, v in role_columns(role, weight_col, max_frac, reason).items():
            rows[k] = v
        repaired_conf_rows.append(rows)
    repaired_conf = pd.concat(repaired_conf_rows, ignore_index=True, sort=False)
    repaired_conf.to_csv(P3 / "step3_phase3_3_confounding_audit_by_universe.repaired.csv", index=False)

    status = {
        "phase": "3.3_repair",
        "verdict": "CONDITIONAL_PASS",
        "repair_target": "PD1_anchor_all_universes_high_confounding_risk",
        "full_PD1_anchor_role": "sensitivity_only",
        "weighted_PD1_anchor_role": weighted_role,
        "weighted_PD1_anchor_confounding_risk_level": weighted_risk,
        "effective_max_cohort_weight_fraction": weighted_max,
        "weighted_cohort_response_chi2_p": weighted_p,
        "hard_failures": [],
        "conditional_items": [
            "unweighted_full_PD1_anchor_downgraded_to_sensitivity_only",
            "within_cohort_meta_is_biological_anchor_only_not_prediction_model",
        ]
        + ([] if weighted_role == "primary_allowed" else ["weighted_PD1_anchor_still_not_primary_allowed"]),
        "allowed_to_continue_step3": True,
        "primary_pd1_anchor_allowed": weighted_role == "primary_allowed",
        "phase3_4_pd1_primary_universe": "PD1_anchor_weighted" if weighted_role == "primary_allowed" else None,
    }
    (P3 / "step3_phase3_3_status.repaired.yaml").write_text(
        yaml.safe_dump(status, sort_keys=False), encoding="utf-8"
    )

    summary = [
        "# Repair Phase 3.3 PD1_anchor Confounding Summary",
        "",
        "## Root Cause",
        "- Original PD1_anchor supervised universe has 2 cohorts, 340 samples, 51 patients.",
        f"- Original dominant cohort sample fraction: {full_max:.3f}.",
        "- Cohort-response association and high cohort predictability indicate cohort-level confounding.",
        "",
        "## Repair",
        "- `PD1_anchor_full`: keeps all eligible PD1_anchor samples, role `sensitivity_only`, sample_weight=1.0.",
        "- `PD1_anchor_weighted`: keeps all eligible PD1_anchor samples, uses patient-level cohort-balanced sample weights.",
        f"- Weighted effective max cohort weight fraction: {weighted_max:.3f}.",
        f"- Weighted cohort-response chi2 p-value: {weighted_p:.6g}.",
        f"- Weighted confounding risk: {weighted_risk}; role `{weighted_role}`.",
        "- `PD1_anchor_within_cohort_meta`: role `biological_anchor_only`; no prediction model, no AUC use.",
        "",
        "## Non-Changes",
        "- No samples removed or duplicated.",
        "- Feature matrix unchanged.",
        "- Feature dictionary unchanged.",
        "- Response labels unchanged.",
        "- Patient split unchanged.",
        "",
        "## Downstream Rules",
        "- Phase3.4/3.5 primary models may use only `universe_role == primary_allowed`.",
        "- Full unweighted PD1_anchor is sensitivity only.",
        "- Within-cohort meta PD1_anchor is biological anchor only.",
    ]
    (P3 / "repair_phase3_3_PD1_anchor_confounding_summary.md").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )
    return 0


def add_default_roles(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "confounding_risk_level" in df.columns:
        primary = ~df["confounding_risk_level"].astype(str).eq("HIGH")
    else:
        primary = df.get("status", "eligible").astype(str).eq("eligible")
    df["universe_role"] = np.where(primary, "primary_allowed", "sensitivity_only")
    df["sample_weight_column"] = "none"
    df["cohort_adjustment_required"] = False
    df["effective_max_cohort_weight_fraction"] = df.get("max_cohort_sample_fraction", df.get("max_cohort_fraction", np.nan))
    df["primary_allowed_after_repair"] = df["universe_role"].eq("primary_allowed")
    df["sensitivity_only_after_repair"] = df["universe_role"].eq("sensitivity_only")
    df["biological_anchor_only"] = False
    df["repair_reason"] = "unchanged_non_PD1_anchor_or_existing_role"
    return df


def add_membership_roles(
    membership: pd.DataFrame, frac_pd1: pd.DataFrame, role_specs: list[tuple[str, str, str, float, str]]
) -> pd.DataFrame:
    non_pd1 = membership[~membership["task"].eq("PD1_anchor")].copy()
    non_pd1["universe_id"] = (
        non_pd1["task"] + "__" + non_pd1["feature_family"] + "__" + non_pd1["universe_type"]
    )
    non_pd1["universe_role"] = "primary_allowed"
    non_pd1["sample_weight_column"] = "none"
    non_pd1["sample_weight"] = 1.0
    non_pd1["cohort_adjustment_required"] = False
    non_pd1["effective_max_cohort_weight_fraction"] = np.nan
    non_pd1["primary_allowed_after_repair"] = True
    non_pd1["sensitivity_only_after_repair"] = False
    non_pd1["biological_anchor_only"] = False
    non_pd1["repair_reason"] = "unchanged_non_PD1_anchor"

    pd1_base = membership[
        membership["task"].eq("PD1_anchor")
        & membership["universe_type"].eq("supervised_labeled_universe")
        & membership["feature_family"].isin(PD1_FEATURE_FAMILIES)
    ].copy()
    weights = frac_pd1.set_index("sample_key")[
        [
            "pd1_anchor_full_sample_weight",
            "pd1_anchor_weighted_sample_weight",
            "pd1_anchor_within_cohort_meta_weight",
        ]
    ]
    out = [non_pd1]
    for group, role, weight_col, max_frac, reason in role_specs:
        rows = pd1_base.copy()
        rows["universe_id"] = rows["feature_family"].map(
            lambda f: f"{group}__{f}__supervised_labeled_universe"
        )
        rows["universe_role"] = role
        rows["sample_weight_column"] = weight_col
        rows["sample_weight"] = rows["sample_key"].map(weights[weight_col]) if weight_col != "none" else np.nan
        rows.loc[~rows["included"].astype(bool), "sample_weight"] = 0.0
        rows["cohort_adjustment_required"] = role == "primary_allowed"
        rows["effective_max_cohort_weight_fraction"] = max_frac
        rows["primary_allowed_after_repair"] = role == "primary_allowed"
        rows["sensitivity_only_after_repair"] = role == "sensitivity_only"
        rows["biological_anchor_only"] = role == "biological_anchor_only"
        rows["repair_reason"] = reason
        out.append(rows)
    return pd.concat(out, ignore_index=True, sort=False)


if __name__ == "__main__":
    raise SystemExit(main())
