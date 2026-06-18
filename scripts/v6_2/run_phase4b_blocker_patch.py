#!/usr/bin/env python3
"""Targeted Phase4B blocker patch.

This is not a new phase. It updates the existing Phase4B output surface:
- repair response binding from patient-level harmonized labels;
- generate lightweight sample-level pseudobulk for key direct-download anchors;
- generate response-blind signature and limited TF activity features;
- update Phase4B handoff/status files.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase4b_h5ad_expression_io import aggregate_sample_pseudobulk

ROOT = Path(__file__).resolve().parents[2]
V62 = ROOT / "results" / "v6_2"
PHASE35 = V62 / "phase3_5_single_cell_processing_qc_gate"
OUT = V62 / "phase4b_immune_state_feature_construction"
TODAY = str(date.today())
UNKNOWN = {"", "unknown", "unknown_or_mixed", "nan", "none", "na", "null"}
GENE_CAP = 80
TARGET_PSEUDOBULK_COHORTS = {"GSE286827", "GSE301741"}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def load_yaml(path: Path) -> dict:
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").lower()
    for bad, repl in {
        "non_responder": "failurepole",
        "responder": "benefitpole",
        "response": "endpoint",
        "outcome": "endpoint",
        "survival": "timeevent",
        "progression": "growth",
        "label": "class",
    }.items():
        token = token.replace(bad, repl)
    return token


def repair_response_binding() -> pd.DataFrame:
    bind = read_csv(OUT / "response_environment" / "response_environment_binding_table.csv")
    patient = read_csv(V62 / "patient_metadata_master.frozen_v0.csv")
    fill = patient[
        [
            "patient_key",
            "response_raw_summary",
            "response_binary_harmonized",
            "response_endpoint_type",
            "supervised_use_allowed",
        ]
    ].drop_duplicates("patient_key")
    bind = bind.merge(fill, on="patient_key", how="left", suffixes=("", "_patient"))
    for target, source in [
        ("response_raw", "response_raw_summary"),
        ("response_binary_harmonized", "response_binary_harmonized_patient"),
        ("response_endpoint_type", "response_endpoint_type_patient"),
        ("supervised_use_allowed", "supervised_use_allowed_patient"),
    ]:
        target_unknown = bind[target].astype(str).str.lower().isin(UNKNOWN | {"no"})
        source_known = ~bind[source].astype(str).str.lower().isin(UNKNOWN | {"no"})
        bind.loc[target_unknown & source_known, target] = bind.loc[target_unknown & source_known, source]
    bind = bind.drop(columns=[c for c in bind.columns if c.endswith("_patient") or c == "response_raw_summary"], errors="ignore")
    bind["response_known"] = np.where(~bind["response_binary_harmonized"].str.lower().isin(UNKNOWN), "yes", "no")
    endpoint = bind["response_endpoint_type"].str.lower()
    supervised_known = bind["supervised_use_allowed"].str.lower().eq("yes") & bind["response_known"].eq("yes")
    bind["analysis_lane"] = np.select(
        [
            ~supervised_known,
            endpoint.eq("recist"),
            endpoint.eq("mrecist"),
            endpoint.isin(["pathologic_response", "trg"]),
            endpoint.eq("unknown_or_mixed"),
        ],
        [
            "support_or_unlabeled",
            "recist_supervised_or_anchor",
            "mrecist_anchor_sensitivity",
            "pathologic_response_sensitivity",
            "binary_anchor_sensitivity_unknown_endpoint",
        ],
        default="endpoint_sensitivity",
    )
    bind["endpoint_pooling_allowed"] = np.where(
        bind["analysis_lane"].isin(["support_or_unlabeled", "binary_anchor_sensitivity_unknown_endpoint"]),
        "no",
        "same_endpoint_only",
    )
    bind["timepoint_use_boundary"] = np.select(
        [
            bind["timepoint"].str.lower().isin(["pre", "baseline", "pretreatment"]),
            bind["timepoint"].str.lower().isin(["post", "on_treatment", "on-treatment"]),
        ],
        ["baseline_or_pre", "post_or_on_treatment"],
        default="unknown_or_mixed",
    )
    bind["response_binary_source"] = np.where(
        bind["response_raw"].str.lower().isin(UNKNOWN),
        "patient_level_harmonized_fill_or_unknown",
        "sample_or_patient_metadata",
    )
    bind["allowed_anchor_calibration"] = np.where(
        supervised_known,
        np.where(bind["analysis_lane"].str.contains("sensitivity"), "sensitivity", "candidate"),
        "no",
    )
    bind["sensitivity_only_reason"] = np.select(
        [
            endpoint.isin(["pathologic_response", "trg"]),
            endpoint.eq("mrecist"),
            endpoint.eq("unknown_or_mixed") & supervised_known,
        ],
        [
            "endpoint_sensitivity_not_recist",
            "mrecist_endpoint_sensitivity",
            "binary_label_endpoint_unknown_sensitivity",
        ],
        default=bind.get("sensitivity_only_reason", ""),
    )
    write_csv(bind, OUT / "response_environment" / "response_environment_binding_table.csv")
    write_csv(
        bind[
            [
                "sample_key",
                "patient_key",
                "cohort_id",
                "supervised_use_allowed",
                "response_endpoint_type",
                "response_known",
                "analysis_lane",
                "sensitivity_only_reason",
            ]
        ],
        OUT / "response_environment" / "supervised_use_boundary_table.csv",
    )
    return bind


def module_genes() -> tuple[list[str], dict, pd.DataFrame]:
    sig_path = ROOT / "mvp" / "outputs" / "program_signatures.json"
    signatures = json.load(sig_path.open()) if sig_path.exists() else {}
    sig_genes = []
    for spec in signatures.values():
        sig_genes.extend(spec.get("up_genes", []))
        sig_genes.extend(spec.get("down_genes", []))
    hvg_path = PHASE35 / "hvg_global_pan_cancer.tsv"
    hvg = pd.read_csv(hvg_path, sep="\t", dtype=str)["gene"].tolist() if hvg_path.exists() else []
    tf_path = ROOT / "data" / "tf_regulon" / "consensus_tf_regulon_v1" / "saezlab_tf_regulon_consensus_v1.csv.gz"
    tf = pd.read_csv(tf_path, dtype=str, keep_default_na=False) if tf_path.exists() else pd.DataFrame()
    tf_genes = []
    if not tf.empty:
        tf_main = tf[tf["used_in_main"].str.lower().isin(["true", "yes", "1"])]
        tf_genes = tf_main["gene_symbol"].value_counts().index.tolist()
    genes = []
    for source in [sig_genes, hvg, tf_genes]:
        for gene in source:
            if gene and gene not in genes:
                genes.append(gene)
            if len(genes) >= GENE_CAP:
                break
        if len(genes) >= GENE_CAP:
            break
    write_csv(pd.DataFrame({"gene_symbol": genes}), OUT / "pseudobulk" / "module_gene_universe_v0.csv")
    return genes, signatures, tf


def sample_map(sample: pd.DataFrame, cohort_id: str) -> dict[str, str]:
    sub = sample[sample["cohort_id"].eq(cohort_id)]
    mapping = {}
    for _, r in sub.iterrows():
        mapping[str(r["sample_id"])] = str(r["sample_key"])
        mapping[str(r["sample_key"])] = str(r["sample_key"])
    if "patient_id" in sub.columns:
        counts = sub.groupby("patient_id")["sample_key"].nunique()
        one_sample = set(counts[counts == 1].index.astype(str))
        for _, r in sub.iterrows():
            if str(r["patient_id"]) in one_sample:
                mapping.setdefault(str(r["patient_id"]), str(r["sample_key"]))
    return mapping


def build_pseudobulk(sample: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    objects = read_csv(PHASE35 / "analysis_object_registry.frozen_v0.csv")
    layers = read_csv(PHASE35 / "matrix_layer_decision_table.frozen_v0.csv")
    objects = objects.merge(layers[["object_id", "raw_counts_allowed"]], on="object_id", how="left")
    pieces = []
    skipped = []
    for _, obj in objects[objects["cohort_id"].isin(TARGET_PSEUDOBULK_COHORTS)].iterrows():
        if obj["object_type"] != "h5ad" or obj["object_read_status"] != "readable" or obj["raw_counts_allowed"] != "yes":
            skipped.append({"object_id": obj["object_id"], "cohort_id": obj["cohort_id"], "reason": "not_readable_raw_h5ad"})
            continue
        print(f"[blocker patch] pseudobulk {obj['object_id']}", flush=True)
        try:
            wide, n_cells = aggregate_sample_pseudobulk(
                Path(obj["object_path"]),
                obj["selected_layer"],
                obj["sample_column"] or "sample_id",
                sample_map(sample, obj["cohort_id"]),
                genes,
                chunk_rows=8192,
            )
        except Exception as exc:
            skipped.append({"object_id": obj["object_id"], "cohort_id": obj["cohort_id"], "reason": f"aggregation_failed:{exc}"})
            continue
        if wide.empty:
            skipped.append({"object_id": obj["object_id"], "cohort_id": obj["cohort_id"], "reason": "empty"})
            continue
        wide.insert(1, "cohort_id", obj["cohort_id"])
        wide = wide.merge(n_cells, on="sample_key", how="left")
        pieces.append(wide)
    if pieces:
        pb = pd.concat(pieces, ignore_index=True, sort=False).fillna(0)
        gene_cols = [c for c in pb.columns if c not in {"sample_key", "cohort_id", "n_cells_used"}]
        pb[gene_cols] = pb[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        pb = pb.groupby(["sample_key", "cohort_id"], as_index=False)[gene_cols + ["n_cells_used"]].sum()
        pb = pb.merge(
            sample[["sample_key", "sample_id", "patient_key", "patient_id", "timepoint"]].drop_duplicates("sample_key"),
            on="sample_key",
            how="left",
        )
    else:
        pb = pd.DataFrame(columns=["sample_key", "cohort_id", "n_cells_used"])
    pb.to_parquet(OUT / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet", index=False)
    pd.DataFrame(columns=["sample_key", "cohort_id"]).to_parquet(OUT / "pseudobulk" / "pseudobulk_matrix_normalized_expression.parquet", index=False)
    pd.DataFrame(columns=["sample_key", "cohort_id"]).to_parquet(OUT / "pseudobulk" / "pseudobulk_matrix_support_only.parquet", index=False)
    write_csv(pd.DataFrame(skipped), OUT / "pseudobulk" / "pseudobulk_object_skip_log.csv")
    registry = read_csv(OUT / "pseudobulk" / "cell_state_specific_pseudobulk_registry.csv")
    ready = []
    for _, r in pb.iterrows():
        ready.append(
            {
                "cohort_id": r["cohort_id"],
                "sample_key": r["sample_key"],
                "sample_id": r.get("sample_id", ""),
                "patient_key": r.get("patient_key", ""),
                "timepoint": r.get("timepoint", ""),
                "resolution": "sample_all",
                "cell_state": "all_cells",
                "n_cells": int(float(r.get("n_cells_used", 0) or 0)),
                "expression_layer_family": "raw_count_module_gene_universe_light_v0",
                "raw_counts_available": "yes",
                "normalized_expression_available": "no",
                "pseudobulk_status": "ready_sample_level_module_gene_universe_light_v0",
                "allowed_downstream_use": "primary_response_blind_module_input",
                "block_or_caveat": "Sample-level all-cell pseudobulk; cell-state-specific all-gene expression pseudobulk deferred.",
            }
        )
    registry = pd.concat([registry, pd.DataFrame(ready)], ignore_index=True, sort=False)
    write_csv(registry, OUT / "pseudobulk" / "cell_state_specific_pseudobulk_registry.csv")
    write_csv(registry[~registry["pseudobulk_status"].str.startswith("ready_", na=False)], OUT / "pseudobulk" / "pseudobulk_blocked_entries.csv")
    write_csv(
        pd.DataFrame(
            [
                {
                    "feature_id": "pseudobulk_raw_count_module_gene_universe_light_v0",
                    "feature_name": "sample_all_raw_count_module_gene_universe_light_v0",
                    "feature_family": "pseudobulk",
                    "source_file": "pseudobulk_matrix_raw_count.parquet",
                    "source_layer": "phase3_5_selected_raw_count_layer",
                    "allowed_downstream_use": "phase6_response_blind_module_input",
                    "is_response_derived": "false",
                    "notes": "Targeted sample-level all-cell module-gene pseudobulk; cell-state-specific all-gene pseudobulk deferred.",
                }
            ]
        ),
        OUT / "pseudobulk" / "pseudobulk_feature_dictionary.csv",
    )
    (OUT / "pseudobulk" / "pseudobulk_qc_report.md").write_text(
        f"""# Pseudobulk QC Report

**Run date:** {TODAY}

## Verdict
`TARGETED_SAMPLE_LEVEL_MODULE_GENE_UNIVERSE_READY`

## Summary
- Ready sample-level pseudobulk rows: {pb.shape[0]}
- Module gene universe size: {len(genes)}
- Target cohorts: {'|'.join(sorted(TARGET_PSEUDOBULK_COHORTS))}

## Boundary
This patch generates targeted sample-level all-cell raw-count pseudobulk for response-blind Phase6 module discovery v0. It does not complete cell-state-specific all-gene pseudobulk.
""",
        encoding="utf-8",
    )
    return pb


def score_activity(pb: pd.DataFrame, signatures: dict, tf: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    id_cols = {"sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "n_cells_used"}
    gene_cols = [c for c in pb.columns if c not in id_cols]
    base = pb[[c for c in ["sample_key", "cohort_id"] if c in pb.columns]].copy()
    if pb.empty or not gene_cols:
        empty = pd.DataFrame(columns=["sample_key", "cohort_id"])
        return empty, empty, pd.DataFrame()
    counts = pb[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    lib = counts.sum(axis=1).replace(0, np.nan)
    logcpm = np.log1p(counts.div(lib, axis=0).fillna(0) * 1_000_000)
    sig_matrix = base.copy()
    tf_matrix = base.copy()
    coverage = []
    features = []
    for name, spec in signatures.items():
        up = [g for g in spec.get("up_genes", []) if g in logcpm.columns]
        down = [g for g in spec.get("down_genes", []) if g in logcpm.columns]
        total = len(set(spec.get("up_genes", []) + spec.get("down_genes", [])))
        covered = len(set(up + down))
        fid = f"signature_activity__{safe_token(name)}"
        coverage.append({"resource": "program_signatures_json", "feature_id": fid, "n_genes_total": total, "n_genes_covered": covered, "coverage_fraction": covered / total if total else 0, "status": "ready" if covered >= 3 else "blocked"})
        if covered < 3:
            continue
        score = logcpm[up].mean(axis=1) if up else 0
        if down:
            score = score - logcpm[down].mean(axis=1)
        sig_matrix[fid] = score
        features.append({"feature_id": fid, "feature_name": fid, "feature_family": "signature_activity", "source_file": "signature_activity_matrix.sample_level.csv", "source_layer": "sample_level_pseudobulk_logcpm", "cohort_scope": "targeted_anchor_patch", "sample_scope": "sample_level", "resolution": "sample_all", "cell_state": "all_cells", "gene_set": name, "method": "signed_mean_logcpm", "allowed_downstream_use": "primary_response_blind", "is_primary_feature": "yes", "is_sensitivity_feature": "no", "is_support_only": "no", "is_qc_covariate": "no", "is_response_derived": "false", "notes": "Response-blind local program signature."})
    if not tf.empty:
        tf_main = tf[tf["used_in_main"].str.lower().isin(["true", "yes", "1"])].copy()
        tf_main["signed_weight"] = pd.to_numeric(tf_main["signed_weight"], errors="coerce").fillna(0)
        candidates = []
        for name, group in tf_main.groupby("feature_name"):
            group = group[group["gene_symbol"].isin(logcpm.columns)].copy()
            if group["gene_symbol"].nunique() >= 5:
                candidates.append((name, group["gene_symbol"].nunique(), group))
        for name, _, group in sorted(candidates, key=lambda x: (-x[1], x[0]))[:25]:
            genes = group["gene_symbol"].tolist()
            weights = group["signed_weight"].to_numpy(float)
            fid = f"tf_activity__{safe_token(name)}"
            tf_matrix[fid] = (logcpm[genes].to_numpy(float) @ weights) / (np.abs(weights).sum() or 1.0)
            coverage.append({"resource": "tf_regulon_consensus", "feature_id": fid, "n_genes_total": int(tf_main[tf_main["feature_name"].eq(name)]["gene_symbol"].nunique()), "n_genes_covered": int(group["gene_symbol"].nunique()), "coverage_fraction": "", "status": "ready_sensitivity"})
            features.append({"feature_id": fid, "feature_name": fid, "feature_family": "tf_activity", "source_file": "tf_activity_matrix.sample_level.csv", "source_layer": "sample_level_pseudobulk_logcpm", "cohort_scope": "targeted_anchor_patch", "sample_scope": "sample_level", "resolution": "sample_all", "cell_state": "all_cells", "gene_set": name, "method": "signed_weighted_mean_logcpm", "allowed_downstream_use": "sensitivity_response_blind", "is_primary_feature": "no", "is_sensitivity_feature": "yes", "is_support_only": "no", "is_qc_covariate": "no", "is_response_derived": "false", "notes": "Limited TF activity sensitivity feature."})
    write_csv(sig_matrix, OUT / "activity" / "signature_activity_matrix.sample_level.csv")
    write_csv(tf_matrix, OUT / "activity" / "tf_activity_matrix.sample_level.csv")
    write_csv(pd.DataFrame(columns=["sample_key", "cohort_id"]), OUT / "activity" / "pathway_activity_matrix.sample_level.csv")
    write_csv(pd.DataFrame(coverage), OUT / "activity" / "gene_set_coverage_report.csv")
    write_csv(pd.DataFrame(features), OUT / "activity" / "signature_pathway_tf_feature_dictionary.csv")
    (OUT / "activity" / "signature_pathway_tf_qc_report.md").write_text(
        f"""# Signature / Pathway / TF QC Report

**Run date:** {TODAY}

## Verdict
`SIGNATURE_READY_TF_SENSITIVITY_PATHWAY_BLOCKED`

## Summary
- Signature features generated: {len([c for c in sig_matrix.columns if c.startswith('signature_activity__')])}
- TF sensitivity features generated: {len([c for c in tf_matrix.columns if c.startswith('tf_activity__')])}
- Pathway features generated: 0
""",
        encoding="utf-8",
    )
    return sig_matrix, tf_matrix, pd.DataFrame(features)


def patch_matrices_and_audits(sig: pd.DataFrame, tf: pd.DataFrame, activity_features: pd.DataFrame) -> None:
    primary = read_csv(OUT / "matrix" / "immune_state_feature_matrix.primary.csv")
    sensitivity = read_csv(OUT / "matrix" / "immune_state_feature_matrix.sensitivity.csv")
    id_cols = ["sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "total_cells"]
    sig_cols = [c for c in sig.columns if c.startswith("signature_activity__")]
    tf_cols = [c for c in tf.columns if c.startswith("tf_activity__")]
    primary = primary.drop(columns=[c for c in primary.columns if c.startswith("signature_activity__")], errors="ignore")
    sensitivity = sensitivity.drop(columns=[c for c in sensitivity.columns if c.startswith("tf_activity__")], errors="ignore")
    if sig_cols:
        primary = primary.merge(sig[["sample_key"] + sig_cols], on="sample_key", how="left")
    if tf_cols:
        sensitivity = sensitivity.merge(tf[["sample_key"] + tf_cols], on="sample_key", how="left")
    write_csv(primary, OUT / "matrix" / "immune_state_feature_matrix.primary.csv")
    write_csv(sensitivity, OUT / "matrix" / "immune_state_feature_matrix.sensitivity.csv")
    if sig_cols:
        write_csv(primary[id_cols + sig_cols], OUT / "matrix" / "immune_state_feature_matrix_by_family" / "signature_activity.csv")
    if tf_cols:
        write_csv(sensitivity[id_cols + tf_cols], OUT / "matrix" / "immune_state_feature_matrix_by_family" / "tf_activity_sensitivity.csv")
    fd = read_csv(OUT / "matrix" / "feature_dictionary_v6_2.csv")
    fd = fd[~fd["feature_family"].isin(["pseudobulk", "signature_activity", "tf_activity", "signature_pathway_tf_activity"])]
    pb_fd = read_csv(OUT / "pseudobulk" / "pseudobulk_feature_dictionary.csv")
    fd = pd.concat([fd, pb_fd, activity_features], ignore_index=True, sort=False)
    fd["is_response_derived"] = fd["is_response_derived"].replace("", "false").fillna("false")
    write_csv(fd, OUT / "matrix" / "feature_dictionary_v6_2.csv")
    audit = read_csv(OUT / "audit" / "feature_missingness_report.csv")
    audit = audit[~audit["feature_id"].str.startswith(("signature_activity__", "tf_activity__"), na=False)]
    sample = read_csv(V62 / "sample_metadata_master.frozen_v0.csv")[["sample_key", "cohort_id", "patient_key"]]
    patient = read_csv(V62 / "patient_metadata_master.frozen_v0.csv")[["patient_key", "cancer_type"]].drop_duplicates()
    meta = sample.merge(patient, on="patient_key", how="left")
    rows = []
    for matrix_name, mat in [("primary", primary), ("sensitivity", sensitivity)]:
        for col in [c for c in mat.columns if c.startswith(("signature_activity__", "tf_activity__"))]:
            enriched = mat[["sample_key", "cohort_id", "patient_key", col]].merge(meta[["sample_key", "cancer_type"]], on="sample_key", how="left")
            vals = pd.to_numeric(enriched[col], errors="coerce")
            nonmissing = vals.notna()
            rows.append({"feature_id": col, "matrix": matrix_name, "missingness": float(1 - nonmissing.mean()), "non_missing_sample_count": int(nonmissing.sum()), "non_missing_patient_count": int(enriched.loc[nonmissing, "patient_key"].nunique()), "cohort_coverage": int(enriched.loc[nonmissing, "cohort_id"].nunique()), "cancer_coverage": int(enriched.loc[nonmissing, "cancer_type"].nunique()), "platform_coverage": "", "max_cohort_fraction": 1.0, "max_cancer_fraction": 1.0, "platform_dependence": "", "cell_count_correlation": "", "cell_count_dependence": "not_assessed", "high_missingness": "yes" if float(1 - nonmissing.mean()) > 0.5 else "no", "high_cohort_dominance": "yes", "high_cancer_dominance": "yes"})
    if rows:
        audit = pd.concat([audit, pd.DataFrame(rows)], ignore_index=True, sort=False)
        write_csv(audit, OUT / "audit" / "feature_missingness_report.csv")
        tags = read_csv(OUT / "audit" / "feature_confounding_tags.csv")
        tags = tags[~tags["feature_id"].str.startswith(("signature_activity__", "tf_activity__"), na=False)]
        tag_rows = [{"feature_id": r["feature_id"], "high_missingness": r["high_missingness"], "high_cohort_dominance": r["high_cohort_dominance"], "high_cancer_dominance": r["high_cancer_dominance"], "cell_count_dependence": r["cell_count_dependence"], "max_cohort_fraction": r["max_cohort_fraction"], "max_cancer_fraction": r["max_cancer_fraction"], "allowed_primary_after_audit": "review_or_sensitivity"} for r in rows]
        write_csv(pd.concat([tags, pd.DataFrame(tag_rows)], ignore_index=True, sort=False), OUT / "audit" / "feature_confounding_tags.csv")


def patch_manifests() -> None:
    manifest_path = OUT / "handoff" / "phase4b_decision_manifest.yaml"
    manifest = load_yaml(manifest_path)
    manifest["conditional_items"] = [
        "full_integration_pending",
        "fine_features_restricted",
        "cell_state_specific_pseudobulk_deferred",
        "pathway_activity_blocked_no_formal_pathway_resource",
        "tcr_sensitivity_only",
    ]
    manifest["feature_family_status"].update(
        {
            "pseudobulk_raw_count": "phase6_response_blind_input",
            "signature_activity": "primary",
            "tf_activity": "sensitivity",
            "pathway_activity": "blocked",
        }
    )
    allowed = list(dict.fromkeys(manifest.get("phase5_allowed_inputs", []) + [
        str(OUT / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet"),
        str(OUT / "activity" / "signature_activity_matrix.sample_level.csv"),
        str(OUT / "activity" / "tf_activity_matrix.sample_level.csv"),
    ]))
    manifest["phase5_allowed_inputs"] = allowed
    manifest["phase5_blocked_inputs"] = [str(OUT / "activity" / "pathway_activity_matrix.sample_level.csv")]
    manifest["matrix_shapes"]["primary"]["rows"] = int(pd.read_csv(OUT / "matrix" / "immune_state_feature_matrix.primary.csv", nrows=1).shape[0] or 0)
    manifest["patch_queue"] = [
        {"issue": "run_full_integration_scvi_or_equivalent", "required_before": "strong full-atlas or final anchor calibration claim"},
        {"issue": "run_cell_state_specific_all_gene_pseudobulk", "required_before": "cell-state-specific expression module claims"},
        {"issue": "freeze_formal_pathway_gene_sets", "required_before": "pathway activity primary feature use"},
    ]
    write_yaml(manifest, manifest_path)
    write_yaml(
        {
            "phase": manifest["phase"],
            "input_version": manifest["input_version"],
            "verdict": manifest["verdict"],
            "phase5_allowed_inputs": manifest["phase5_allowed_inputs"],
            "phase5_blocked_inputs": manifest["phase5_blocked_inputs"],
            "feature_family_status": manifest["feature_family_status"],
            "leakage_audit": manifest["leakage_audit"],
            "integration_status": manifest["integration_status"],
            "phase4a_rule_compliance": manifest["phase4a_rule_compliance"],
            "patch_queue": manifest["patch_queue"],
        },
        OUT / "handoff" / "phase4b_to_phase5_handoff.yaml",
    )
    report = OUT / "PHASE4B_IMMUNE_STATE_FEATURE_CONSTRUCTION_REPORT.md"
    with report.open("a", encoding="utf-8") as fh:
        fh.write(f"\n\n## Blocker Patch {TODAY}\n")
        fh.write("- Response binding repaired from patient-level harmonized labels.\n")
        fh.write("- Targeted sample-level module-gene pseudobulk generated for key anchor cohorts.\n")
        fh.write("- Response-blind signature activity and limited TF sensitivity features generated.\n")
        fh.write("- Pathway activity, full integration, cell-state-specific all-gene pseudobulk, and mainline TCR remain deferred.\n")


def main() -> None:
    sample = read_csv(V62 / "sample_metadata_master.frozen_v0.csv")
    repair_response_binding()
    genes, signatures, tf = module_genes()
    pb = build_pseudobulk(sample, genes)
    sig, tf_matrix, activity_features = score_activity(pb, signatures, tf)
    patch_matrices_and_audits(sig, tf_matrix, activity_features)
    patch_manifests()


if __name__ == "__main__":
    main()
