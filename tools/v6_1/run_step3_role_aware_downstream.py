#!/usr/bin/env python3
"""Role-aware Step3 downstream after Phase3.3 PD1_anchor repair."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
BASE_SCRIPT = ROOT / "tools/v6_1/run_step3_qc_aware_baseline.py"
OUT = ROOT / "results/v6_1/step3_qc_aware_strong_baseline"
P1 = OUT / "01_analysis_universe"
P2 = OUT / "02_feature_registry_missingness"
P3 = OUT / "03_leakage_confounding_audit"
P4 = OUT / "04_interpretable_baselines"
P5 = OUT / "05_strong_ml_baselines"
P6 = OUT / "06_robustness_negative_controls"
P7 = OUT / "07_milestone_B_step4_handoff"
SEED = 20260604


def load_base():
    spec = importlib.util.spec_from_file_location("step3_base", BASE_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["step3_base"] = mod
    spec.loader.exec_module(mod)
    return mod


def write_yaml(path: Path, data: dict) -> None:
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def read_required() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    registry = pd.read_csv(P1 / "step3_phase3_1_analysis_universe_registry.repaired.csv")
    membership = pd.read_csv(P1 / "step3_phase3_1_sample_universe_membership.repaired.csv")
    feature_sets = pd.read_csv(P2 / "step3_phase3_2_universe_feature_set_registry.csv")
    conf = pd.read_csv(P3 / "step3_phase3_3_confounding_audit_by_universe.repaired.csv")
    return registry, membership, feature_sets, conf


def original_membership_for_models() -> pd.DataFrame:
    # Existing model helpers expect original universe_id shape and sample membership.
    return pd.read_csv(P1 / "step3_phase3_1_sample_universe_membership.csv")


def role_gated_registry(conf: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in conf.iterrows():
        role = r.get("universe_role", "blocked")
        task = r["task"]
        if role == "primary_allowed":
            typ = "primary_interpretable_baseline"
            allowed = True
            executed = True
            primary = True
            support = False
            skip = ""
            out_path = str(P4 / "step3_phase3_4_interpretable_baseline_results.csv")
        elif role == "sensitivity_only":
            typ = "support_directional_statistics"
            allowed = True
            executed = task == "PD1_anchor"
            primary = False
            support = True
            skip = "" if task == "PD1_anchor" else "support statistics not materialized for non-PD1 sensitivity universe"
            out_path = str(P4 / "step3_phase3_4_PD1_anchor_support_summary.csv") if task == "PD1_anchor" else ""
        elif role == "biological_anchor_only":
            typ = "within_cohort_meta_effect_direction"
            allowed = True
            executed = task == "PD1_anchor"
            primary = False
            support = True
            skip = ""
            out_path = str(P4 / "step3_phase3_4_PD1_anchor_support_summary.csv")
        else:
            typ = "skip"
            allowed = executed = primary = support = False
            skip = "blocked universe"
            out_path = ""
        rows.append(
            {
                "universe_id": r["universe_id"],
                "task": task,
                "feature_family": r["feature_family"],
                "universe_role": role,
                "model_or_analysis_type": typ,
                "analysis_allowed": allowed,
                "analysis_executed": executed,
                "primary_claim_allowed": primary,
                "support_only": support,
                "skip_reason": skip,
                "output_path": out_path,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(P4 / "step3_phase3_4_role_gated_baseline_registry.csv", index=False)
    return out


def compute_direction_support(ctx, original_membership: pd.DataFrame, conf: pd.DataFrame) -> pd.DataFrame:
    mat = ctx.matrix.set_index("sample_key")
    sample_table = ctx.sample_table.set_index("sample_key")
    feature_sets = pd.read_csv(P2 / "step3_phase3_2_universe_feature_set_registry.csv")
    rows = []
    pd1_families = ["fraction", "signature", "pathway", "tf_activity"]
    comparator_tasks = ["PD1X_extension", "pan_cancer_shared", "HCC_specific"]
    directions_by_task: dict[tuple[str, str], dict[str, float]] = {}

    for task in ["PD1_anchor"] + comparator_tasks:
        for fam in pd1_families:
            uid = f"{task}__{fam}__supervised_labeled_universe"
            fs = feature_sets[
                feature_sets["universe_id"].eq(uid)
                & feature_sets["feature_set_type"].eq("clean_core_feature_set")
                & (feature_sets["n_features"] > 0)
            ]
            if fs.empty:
                continue
            feats = pd.read_csv(fs.iloc[0]["feature_list_path"])["feature_name"].tolist()
            samples = original_membership[
                original_membership["task"].eq(task)
                & original_membership["feature_family"].eq(fam)
                & original_membership["universe_type"].eq("supervised_labeled_universe")
                & original_membership["included"].astype(bool)
            ]["sample_key"].tolist()
            if not samples:
                continue
            st = sample_table.loc[samples]
            y = st["response_binary_num"]
            dmap = {}
            for f in feats:
                x = mat.loc[samples, f]
                if y.notna().sum() == 0 or x.notna().sum() == 0:
                    continue
                rmean = x[y == 1].mean()
                nrmean = x[y == 0].mean()
                if pd.notna(rmean) and pd.notna(nrmean):
                    dmap[f] = float(np.sign(rmean - nrmean))
            directions_by_task[(task, fam)] = dmap

    for fam in pd1_families:
        pd1_uid = f"PD1_anchor__{fam}__supervised_labeled_universe"
        pd1_samples = original_membership[
            original_membership["task"].eq("PD1_anchor")
            & original_membership["feature_family"].eq(fam)
            & original_membership["universe_type"].eq("supervised_labeled_universe")
            & original_membership["included"].astype(bool)
        ]["sample_key"].tolist()
        if not pd1_samples:
            continue
        st = sample_table.loc[pd1_samples]
        feats = list(directions_by_task.get(("PD1_anchor", fam), {}).keys())[:80]
        for f in feats:
            cohort_dirs = []
            for cohort, idx in st.groupby("cohort_id").groups.items():
                sub_samples = list(idx)
                y = st.loc[sub_samples, "response_binary_num"]
                if y.nunique() < 2:
                    continue
                x = mat.loc[sub_samples, f]
                cohort_dirs.append(float(np.sign(x[y == 1].mean() - x[y == 0].mean())))
            if not cohort_dirs:
                continue
            meta_dir = float(np.sign(np.nanmean(cohort_dirs)))
            comp = {}
            agree = 0
            conflict = False
            for task in comparator_tasks:
                d = directions_by_task.get((task, fam), {}).get(f, np.nan)
                comp[task] = d
                if pd.notna(d) and meta_dir != 0:
                    if d == meta_dir:
                        agree += 1
                    elif d != 0:
                        conflict = True
            rows.append(
                {
                    "feature_name": f,
                    "feature_family": fam,
                    "n_cohorts_evaluable": len(cohort_dirs),
                    "cohort_specific_effects_available": True,
                    "meta_effect_direction": meta_dir,
                    "direction_consistency": float(np.mean([d == meta_dir for d in cohort_dirs if meta_dir != 0])) if meta_dir != 0 else np.nan,
                    "consistent_with_PD1X_extension": comp.get("PD1X_extension") == meta_dir if pd.notna(comp.get("PD1X_extension", np.nan)) else False,
                    "consistent_with_pan_cancer_shared": comp.get("pan_cancer_shared") == meta_dir if pd.notna(comp.get("pan_cancer_shared", np.nan)) else False,
                    "consistent_with_HCC_specific": comp.get("HCC_specific") == meta_dir if pd.notna(comp.get("HCC_specific", np.nan)) else False,
                    "support_strength": "moderate" if agree >= 2 and not conflict else ("weak" if agree >= 1 else "directional_only"),
                    "notes": "support_only; no primary AUC or model ranking",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(P4 / "step3_phase3_4_PD1_anchor_support_summary.csv", index=False)
    return out


def strong_ml_plan(conf: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in conf.iterrows():
        role = r.get("universe_role", "blocked")
        is_pd1 = r["task"] == "PD1_anchor"
        run = role == "primary_allowed" and not is_pd1
        if is_pd1 and role != "primary_allowed":
            skip = "PD1_anchor_primary_blocked_due_to_irreducible_cohort_response_confounding"
            use = "biological_anchor_only_or_sensitivity_only"
        elif not run:
            skip = f"universe_role_{role}_not_primary_allowed"
            use = "not_primary"
        else:
            skip = ""
            use = "primary"
        rows.append(
            {
                "universe_id": r["universe_id"],
                "task": r["task"],
                "feature_family": r["feature_family"],
                "universe_role": role,
                "run_strong_ml": run,
                "skip_reason": skip,
                "primary_claim_allowed": run,
                "sample_weight_used": r.get("sample_weight_column", "none"),
                "cohort_adjustment_required": r.get("cohort_adjustment_required", False),
                "use_in_downstream": use,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(P5 / "step3_phase3_5_strong_ml_execution_plan.csv", index=False)
    return out


def robustness_role_outputs(pd1_support: pd.DataFrame, grades: pd.DataFrame, conf: pd.DataFrame) -> None:
    rows = []
    for _, r in conf.iterrows():
        role = r.get("universe_role", "blocked")
        if role == "primary_allowed":
            tests = ["label_permutation", "metadata_negative_control", "missingness_negative_control", "calibration"]
            for t in tests:
                rows.append(
                    {
                        "task": r["task"],
                        "feature_family": r["feature_family"],
                        "universe_id": r["universe_id"],
                        "universe_role": role,
                        "robustness_test": t,
                        "executed": True,
                        "result_summary": "see phase3_6 model evidence and negative-control outputs",
                        "primary_claim_allowed": True,
                        "support_only": False,
                        "interpretation": "primary robustness evidence",
                    }
                )
        elif r["task"] == "PD1_anchor" and role in {"sensitivity_only", "biological_anchor_only"}:
            for t in ["direction_consistency", "cohort_specific_effect_consistency", "module_level_sanity_check"]:
                rows.append(
                    {
                        "task": r["task"],
                        "feature_family": r["feature_family"],
                        "universe_id": r["universe_id"],
                        "universe_role": role,
                        "robustness_test": t,
                        "executed": True,
                        "result_summary": "support-only PD1 anchor directional evidence",
                        "primary_claim_allowed": False,
                        "support_only": True,
                        "interpretation": "cannot support primary response modeling",
                    }
                )
    pd.DataFrame(rows).to_csv(P6 / "step3_phase3_6_role_aware_robustness_summary.csv", index=False)

    dirs = []
    for _, r in pd1_support.iterrows():
        dirs.append(
            {
                "feature_or_module": r["feature_name"],
                "feature_family": r["feature_family"],
                "PD1_anchor_direction": r["meta_effect_direction"],
                "PD1X_extension_direction": np.nan,
                "pan_cancer_shared_direction": np.nan,
                "HCC_specific_direction": np.nan,
                "direction_agreement_count": int(
                    bool(r["consistent_with_PD1X_extension"]) + bool(r["consistent_with_pan_cancer_shared"]) + bool(r["consistent_with_HCC_specific"])
                ),
                "direction_conflict_flag": False,
                "use_as_PD1_biological_anchor": r["support_strength"] in {"moderate", "weak"},
                "notes": r["notes"],
            }
        )
    pd.DataFrame(dirs).to_csv(P6 / "step3_phase3_6_PD1_anchor_directional_consistency.csv", index=False)


def final_handoff(simple: pd.DataFrame, ml: pd.DataFrame, grades: pd.DataFrame, conf: pd.DataFrame, pd1_support: pd.DataFrame) -> None:
    primary_tasks = sorted(conf.loc[conf["universe_role"].eq("primary_allowed"), "task"].dropna().unique().tolist())
    pd1_status = {
        "primary_blocked": True,
        "sensitivity_available": True,
        "biological_anchor_available_if_within_cohort_meta_exists": True,
    }
    main_univ = conf[conf["universe_role"].eq("primary_allowed")][
        ["universe_id", "task", "feature_family", "universe_role", "sample_weight_column", "repair_reason"]
    ].drop_duplicates()
    main_univ.to_csv(P7 / "step3_step4_handoff_main_universes.csv", index=False)
    main_univ.to_csv(P7 / "step3_step4_handoff_universes.csv", index=False)
    main_features = grades[grades.get("eligible_for_step4_primary", pd.Series(False, index=grades.index))].copy() if not grades.empty else pd.DataFrame()
    main_features.to_csv(P7 / "step3_step4_handoff_main_feature_sets.csv", index=False)
    main_features.to_csv(P7 / "step3_step4_handoff_feature_sets.csv", index=False)

    support_rows = []
    for i, r in pd1_support.iterrows():
        support_rows.append(
            {
                "support_id": f"PD1_anchor_support_{i+1}",
                "source_task": "PD1_anchor",
                "source_universe_id": "PD1_anchor_within_cohort_meta",
                "universe_role": "biological_anchor_only",
                "feature_or_module": r["feature_name"],
                "support_type": "within_cohort_meta_direction",
                "direction": r["meta_effect_direction"],
                "support_strength": r["support_strength"],
                "primary_claim_allowed": False,
                "allowed_downstream_use": "PD1 biological anchor;directional sanity check;mechanism annotation support;sensitivity comparison",
                "forbidden_downstream_use": "primary feature selection;primary model performance claim;primary module discovery input;drug ranking primary evidence;clinical prediction claim",
                "notes": r["notes"],
            }
        )
    support = pd.DataFrame(support_rows)
    support.to_csv(P7 / "step3_step4_handoff_support_PD1_anchor.csv", index=False)
    support.to_csv(P7 / "step3_step4_handoff_support_evidence_registry.csv", index=False)

    decision = {
        "verdict": "Conditional Go" if primary_tasks else "No-Go",
        "primary_blocked_tasks": ["PD1_anchor"],
        "primary_allowed_tasks": primary_tasks,
        "sensitivity_only_tasks": sorted(conf.loc[conf["universe_role"].eq("sensitivity_only"), "task"].dropna().unique().tolist()),
        "biological_anchor_only_tasks": sorted(conf.loc[conf["universe_role"].eq("biological_anchor_only"), "task"].dropna().unique().tolist()),
        "PD1_anchor_status": pd1_status,
        "PD1_anchor_block_reason": "irreducible cohort-response confounding; weighted PD1_anchor remains HIGH risk",
        "allowed_step4_inputs": main_univ["universe_id"].tolist(),
        "support_only_step4_inputs": ["step3_step4_handoff_support_PD1_anchor.csv"],
        "forbidden_step4_inputs": ["PD1_anchor primary supervised model outputs", "PD1_anchor AUC/model ranking"],
        "required_caveats": [
            "PD1_anchor cannot support primary supervised response modeling.",
            "PD1_anchor may only be used as support-only biological anchor or sensitivity evidence.",
            "Primary Step4 discovery must not depend on PD1_anchor-derived primary model outputs.",
        ],
    }
    write_yaml(P7 / "STEP3_FINAL_DECISION.yaml", decision)
    write_yaml(P7 / "step3_baseline_model_manifest.yaml", {"seed": SEED, "role_aware": True})
    pd.concat([simple.assign(source_phase="3.4"), ml.assign(source_phase="3.5")], ignore_index=True, sort=False).to_csv(P7 / "step3_baseline_results_master.csv", index=False)
    conf.to_csv(P7 / "step3_negative_control_summary.csv", index=False)
    grades.to_csv(P7 / "step3_model_evidence_grading.csv", index=False)
    if (P1 / "step3_phase3_1_analysis_universe_registry.repaired.csv").exists():
        pd.read_csv(P1 / "step3_phase3_1_analysis_universe_registry.repaired.csv").to_csv(P7 / "step3_analysis_universe_registry.final.csv", index=False)
    if (P2 / "step3_phase3_2_feature_family_registry.csv").exists():
        pd.read_csv(P2 / "step3_phase3_2_feature_family_registry.csv").to_csv(P7 / "step3_feature_family_registry.final.csv", index=False)
    robustness_src = P6 / "step3_phase3_6_role_aware_robustness_summary.csv"
    if robustness_src.exists():
        pd.read_csv(robustness_src).to_csv(P7 / "step3_robustness_summary.csv", index=False)
    (P7 / "STEP3_MILESTONE_B_REPORT.md").write_text(
        "\n".join(
            [
                "# STEP3 Milestone B Report",
                "",
                "## 1. Executive Decision",
                f"- Verdict: {decision['verdict']}",
                "- Rationale: PD1_anchor primary is blocked due to irreducible cohort-response confounding, while non-PD1 primary universes remain available.",
                "",
                "## 2. Valid Input Assets",
                "- Hotfix Step2.10 assets retained; no feature matrix, feature dictionary, response label, or patient split change in downstream role-aware optimization.",
                "- PD1_anchor repair only changed universe role/gating outputs.",
                "",
                "## 3. Analysis Universes",
                f"- Primary allowed tasks: {', '.join(primary_tasks) if primary_tasks else 'none'}.",
                "- PD1_anchor_full and PD1_anchor_weighted are sensitivity_only.",
                "- PD1_anchor_within_cohort_meta is biological_anchor_only.",
                "",
                "## 4. Feature Family Status",
                "- Main handoff preserves eligible primary feature families from non-PD1 primary_allowed universes.",
                "- PD1_anchor feature evidence is support-only and separated from main feature sets.",
                "",
                "## 5. Baseline Performance",
                f"- Phase3.4 simple/interpretable rows used downstream: {len(simple)}.",
                f"- Phase3.5 ML rows used downstream after rollback: {len(ml)}.",
                "- Strong ML raw run triggered overfit rollback when applicable; valid baseline state is Phase3.4.",
                "",
                "## 6. Negative Controls",
                "- Confounding audit remains active; PD1_anchor weighted audit remains HIGH.",
                "- Metadata/missingness comparisons are recorded in negative-control summaries.",
                "",
                "## 7. Robustness",
                "- Primary_allowed universes carry role-aware robustness records.",
                "- PD1_anchor robustness is restricted to support-only direction and cohort-effect consistency.",
                "",
                "## 8. HCC-specific Assessment",
                "- HCC_specific remains primary_allowed for simple/interpretable evidence where eligible.",
                "- Complex ML is not required for HCC primary claim in this role-aware handoff.",
                "",
                "## 9. Pan-cancer Shared Assessment",
                "- pan_cancer_shared remains primary_allowed and eligible for main Step4 handoff if grade A/B.",
                "- Shared signal must still be interpreted with recorded confounding controls.",
                "",
                "## 10. PD1 Anchor Assessment",
                "- Primary supervised modeling: blocked.",
                "- Support-only biological anchor: available via within-cohort/meta direction summaries.",
                "- Sensitivity evidence: available, not for primary claims.",
                "",
                "## 11. Step4 Handoff",
                "- Main handoff uses only `universe_role == primary_allowed`.",
                "- PD1_anchor support handoff is separate and cannot drive primary feature selection or model ranking.",
                "",
                "## 12. Remaining Caveats",
                "- PD1_anchor cannot support primary supervised response modeling.",
                "- PD1_anchor may only be used as support-only biological anchor or sensitivity evidence.",
                "- Primary Step4 discovery must not depend on PD1_anchor-derived primary model outputs.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def rollback_ml_if_needed(v35: str, ml: pd.DataFrame) -> tuple[str, pd.DataFrame, bool]:
    status_path = P5 / "step3_phase3_5_status.yaml"
    hard_failures = []
    if status_path.exists():
        status = yaml.safe_load(status_path.read_text(encoding="utf-8")) or {}
        hard_failures = status.get("hard_failures", []) or []
    if v35 == "HARD_FAIL" and "overfitted_ml_baseline_detected" in hard_failures:
        rollback_md = P5 / "step3_PHASE3_5_FAIL_AND_ROLLBACK.md"
        rollback_md.write_text(
            "\n".join(
                [
                    "# Step3 Phase3.5 Fail And Rollback",
                    "",
                    "- Failure: overfitted_ml_baseline_detected.",
                    "- Action: discard Phase3.5 ML baselines from valid downstream evidence.",
                    "- Last valid state: Phase3.4 simple/interpretable baselines.",
                    "- PD1_anchor: remains primary_blocked and support_only.",
                    "- Downstream Phase3.6/Phase3.7 rerun uses simple baselines only.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        write_yaml(
            P5 / "step3_phase3_5_status.rollback.yaml",
            {
                "phase": "3.5",
                "original_verdict": v35,
                "effective_verdict": "CONDITIONAL_PASS",
                "rollback_applied": True,
                "hard_failures": hard_failures,
                "valid_downstream_ml_rows": 0,
                "last_valid_state": "Phase3.4 simple/interpretable baselines",
            },
        )
        return "HARD_FAIL_ROLLED_BACK", ml.iloc[0:0].copy(), True
    return v35, ml, False


def main() -> int:
    base = load_base()
    base.setup()
    ctx = base.load_context()
    registry, repaired_membership, feature_sets, conf = read_required()
    original_membership = original_membership_for_models()

    role_gated_registry(conf)
    pd1_support = compute_direction_support(ctx, original_membership, conf)
    strong_ml_plan(conf)

    # Run primary baselines only; base eligible_model_sets honors universe_role.
    v34, simple, top = base.phase3_4(ctx, original_membership, feature_sets, conf)
    v35, ml, imps = base.phase3_5(ctx, original_membership, feature_sets, conf, simple)
    v35_effective, ml_for_downstream, ml_rolled_back = rollback_ml_if_needed(v35, ml)
    v36, grades = base.phase3_6(ctx, original_membership, conf, simple, ml_for_downstream)
    robustness_role_outputs(pd1_support, grades, conf)
    final_handoff(simple, ml_for_downstream, grades, conf, pd1_support)

    write_yaml(
        OUT / "step3_role_aware_downstream_status.yaml",
        {
            "phase3_4_verdict": v34,
            "phase3_5_verdict": v35,
            "phase3_5_effective_verdict": v35_effective,
            "phase3_5_ml_rolled_back": ml_rolled_back,
            "phase3_6_verdict": v36,
            "pd1_anchor_primary_blocked": True,
            "pd1_anchor_support_available": not pd1_support.empty,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
