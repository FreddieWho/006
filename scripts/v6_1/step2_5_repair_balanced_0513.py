from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "step2_v6_1_0505_0319"

IMMUNE_LINEAGES = {"T_NK", "Myeloid", "DC", "B", "Plasma", "Mast"}
MID_STATES = [
    "CD8_T",
    "CD4_T",
    "Treg",
    "NK",
    "NKT_like",
    "Mono_FCN1",
    "Macro_C1QC",
    "Macro_SPP1",
    "Macro_inflammatory",
    "Macro_suppressive_candidate",
    "cDC1",
    "cDC2",
    "pDC",
    "B_naive_memory",
    "Plasma",
    "Mast",
    "Immune_unspecified",
]

CATEGORY_TO_MID = {
    "CD8_T": "CD8_T",
    "CD4_T": "CD4_T",
    "Treg": "Treg",
    "NK": "NK",
    "NKT_like": "NKT_like",
    "monocyte": "Mono_FCN1",
    "macrophage": "Macro_C1QC",
    "TAM_like": "Macro_SPP1",
    "cDC1": "cDC1",
    "cDC2": "cDC2",
    "pDC": "pDC",
    "B": "B_naive_memory",
    "Plasma": "Plasma",
    "Mast": "Mast",
}

MARKER_RULE_TO_MID = {
    "marker_CD8_T": "CD8_T",
    "marker_cytotoxic": "CD8_T",
    "marker_CD4_T": "CD4_T",
    "marker_Treg": "Treg",
    "marker_NK": "NK",
    "marker_NKT_like": "NKT_like",
    "marker_monocyte": "Mono_FCN1",
    "marker_macrophage": "Macro_C1QC",
    "marker_TAM_like": "Macro_SPP1",
    "marker_cDC1": "cDC1",
    "marker_cDC2": "cDC2",
    "marker_pDC": "pDC",
    "marker_B": "B_naive_memory",
    "marker_Plasma": "Plasma",
    "marker_Mast": "Mast",
}

CATEGORY_REFS = {
    "pan_T": "REF02;REF03;REF07;REF08",
    "CD8_T": "REF02;REF07;REF08",
    "CD4_T": "REF02;REF03;REF08",
    "Treg": "REF02;REF08;REF11",
    "NK": "REF02;REF03;REF08",
    "NKT_like": "REF02;REF03;REF08",
    "B": "REF02;REF03;REF08",
    "Plasma": "REF02;REF03;REF09;REF13",
    "pan_myeloid": "REF02;REF09;REF10;REF12",
    "monocyte": "REF02;REF09;REF10;REF12",
    "macrophage": "REF09;REF10;REF12",
    "TAM_like": "REF09;REF10;REF12",
    "cDC1": "REF02;REF09;REF10",
    "cDC2": "REF02;REF09;REF10",
    "pDC": "REF02;REF09;REF10",
    "epithelial": "REF09;REF10;REF13",
    "fibroblast": "REF13;REF14",
    "endothelial": "REF02;REF09;REF13",
    "Mast": "REF02;REF03",
    "cycling": "REF08;REF10;REF13",
    "erythroid": "REF02;REF09;REF13",
    "NA": "manual_audit",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def manifest_ref(run_root: Path) -> str:
    patched = run_root / "00_manifest" / "step2_run_manifest.patched.yaml"
    original = run_root / "00_manifest" / "step2_run_manifest.yaml"
    return str(patched if patched.exists() else original)


def normalize_label(x: object) -> str:
    if pd.isna(x):
        return "NA"
    s = str(x).strip()
    return s if s else "NA"


def infer_mid_from_text(major: str, text: str) -> str:
    s = text.lower()
    if major == "T_NK":
        if any(k in s for k in ["treg", "regulatory", "foxp3"]):
            return "Treg"
        if any(k in s for k in ["nkt", "mait", "gdt", "gdt", "gd t", "ilc"]):
            return "NKT_like"
        if "nk" in s and "nkt" not in s and "cd8" not in s:
            return "NK"
        if any(k in s for k in ["cd8", "tex", "cytotoxic", "effector", "dysfunctional", "exhaust", "trm", "stem-like memory"]):
            return "CD8_T"
        if any(k in s for k in ["cd4", " th", "th1", "th2", "th17", "tfh", "helper"]):
            return "CD4_T"
        return "Immune_unspecified"
    if major == "Myeloid":
        if any(k in s for k in ["mono", "cd14", "cd16", "fcn1", "vcan"]):
            return "Mono_FCN1"
        if any(k in s for k in ["spp1", "trem2", "marco"]):
            return "Macro_SPP1"
        if any(k in s for k in ["c1qc", "c1qa", "apoe", "macro", "tam", "kupffer"]):
            return "Macro_C1QC"
        if any(k in s for k in ["inflammatory", "cxcl10", "cxcl9", "il1b"]):
            return "Macro_inflammatory"
        if any(k in s for k in ["suppress", "ccl18"]):
            return "Macro_suppressive_candidate"
        return "Immune_unspecified"
    if major == "DC":
        if any(k in s for k in ["cdc1", "clec9a", "xcr1"]):
            return "cDC1"
        if any(k in s for k in ["cdc2", "cd1c", "fcera", "fcer1a", "lamp3"]):
            return "cDC2"
        if any(k in s for k in ["pdc", "plasmacytoid"]):
            return "pDC"
        return "Immune_unspecified"
    if major == "B":
        return "B_naive_memory"
    if major == "Plasma":
        return "Plasma"
    if major == "Mast":
        return "Mast"
    return "NA"


def infer_functional(mid: str, text: str) -> str:
    s = text.lower()
    if mid == "CD8_T":
        if any(k in s for k in ["exhaust", "tex", "dysfunctional", "layn", "pdcd1", "havcr2"]):
            return "CD8_exhausted"
        if any(k in s for k in ["cytotoxic", "effector", "gzmb", "prf1", "nkg7"]):
            return "CD8_cytotoxic"
        if any(k in s for k in ["prolif", "mki67", "top2a"]):
            return "CD8_proliferating"
    if mid == "CD4_T" and any(k in s for k in ["helper", "th1", "th2", "th17", "tfh"]):
        return "CD4_helper"
    if mid == "Treg" and any(k in s for k in ["suppress", "ctla4", "foxp3"]):
        return "Treg_suppressive"
    if mid in {"Mono_FCN1", "Macro_C1QC", "Macro_SPP1", "Macro_inflammatory", "Macro_suppressive_candidate"}:
        if any(k in s for k in ["inflammatory", "cxcl10", "cxcl9", "il1b", "ifn"]):
            return "Myeloid_inflammatory"
        if any(k in s for k in ["suppress", "ccl18", "arg1", "spp1", "trem2", "marco"]):
            return "Myeloid_suppressive"
    if mid in {"cDC1", "cDC2", "pDC"} and any(k in s for k in ["apc", "hla", "lamp3"]):
        return "DC_APC_high"
    return "NA"


def expand_rules(rules: pd.DataFrame, created_at: str, mref: str) -> pd.DataFrame:
    rules = rules.copy()
    if "reference_ids" not in rules.columns:
        rules["reference_ids"] = rules["marker_category"].map(CATEGORY_REFS).fillna("manual_audit")
    rules["created_at"] = created_at
    rules["input_manifest_ref"] = mref

    rows = [rules]
    mid_rows = []
    func_rows = []
    next_id = len(rules) + 1
    for _, r in rules.iterrows():
        cat = str(r.get("marker_category", ""))
        mid = CATEGORY_TO_MID.get(cat)
        if mid is not None:
            rr = r.copy()
            rr["rule_id"] = f"RULE_MID_{next_id:05d}"
            rr["target_level"] = "immune_mid_state"
            rr["target_label"] = mid
            rr["reference_ids"] = CATEGORY_REFS.get(cat, str(r.get("reference_ids", "manual_audit")))
            rr["notes"] = f"Derived immune_mid_state from author label category: {cat}"
            mid_rows.append(rr)
            next_id += 1
        func = infer_functional(mid or "", str(r.get("source_label", "")) + " " + str(r.get("notes", "")))
        if func != "NA":
            rr = r.copy()
            rr["rule_id"] = f"RULE_FUNC_{next_id:05d}"
            rr["target_level"] = "functional_state_optional"
            rr["target_label"] = func
            rr["reference_ids"] = CATEGORY_REFS.get(cat, str(r.get("reference_ids", "manual_audit")))
            rr["notes"] = f"Derived functional_state_optional from author label text: {func}"
            func_rows.append(rr)
            next_id += 1
    if mid_rows:
        rows.append(pd.DataFrame(mid_rows))
    if func_rows:
        rows.append(pd.DataFrame(func_rows))
    return pd.concat(rows, ignore_index=True)


def add_table_provenance_csv(path: Path, run_id: str, created_at: str, mref: str) -> None:
    df = pd.read_csv(path)
    for col, val in [("run_id", run_id), ("created_at", created_at), ("input_manifest_ref", mref)]:
        if col in df.columns:
            df[col] = val
        else:
            df.insert(0, col, val)
    df.to_csv(path, index=False)


def main() -> None:
    run_root = ROOT / "results" / "v6_1" / "step2" / RUN_ID
    out = run_root / "04_annotation"
    report_dir = run_root / "10_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()
    mref = manifest_ref(run_root)

    annotation_path = out / "cell_state_annotation_v6_1.parquet"
    balanced_path = run_root / "03_qc" / "final_cell_inclusion_flags.qc0513_balanced.parquet"
    mapping_path = run_root / "01_mapping_audit" / "cell_to_sample_mapping_draft.parquet"
    rules_path = out / "cell_state_mapping_rules_v6_1.csv"
    marker_path = out / "marker_registry_v6_1.csv"

    annot = pd.read_parquet(annotation_path)
    flags = pd.read_parquet(balanced_path)
    flags = flags.drop_duplicates(subset=["source_h5ad", "cell_barcode", "sample_id"]).copy()
    assert len(annot) == len(flags), f"row count mismatch before join: {len(annot)} vs {len(flags)}"

    mapping = pd.read_parquet(
        mapping_path,
        columns=["source_h5ad", "cell_barcode", "mapped_sample_id", "mapping_status"],
    )
    mapping = mapping[mapping["mapping_status"] == "mapped"].copy()
    mapping["_occ"] = mapping.groupby(["source_h5ad", "cell_barcode"]).cumcount()
    annot["_occ"] = annot.groupby(["source_h5ad", "cell_barcode"]).cumcount()
    annot = annot.merge(
        mapping[["source_h5ad", "cell_barcode", "_occ", "mapped_sample_id"]],
        on=["source_h5ad", "cell_barcode", "_occ"],
        how="left",
        validate="many_to_one",
    )
    bad_sample = annot["sample_id"].isna() | annot["sample_id"].astype(str).isin(["", "unknown", "nan", "None"])
    annot.loc[bad_sample & annot["mapped_sample_id"].notna(), "sample_id"] = annot.loc[bad_sample & annot["mapped_sample_id"].notna(), "mapped_sample_id"].astype(str)
    annot = annot.drop(columns=["_occ", "mapped_sample_id"])
    pooled_sources = flags.groupby("source_h5ad")["sample_id"].nunique(dropna=False)
    pooled_sources = set(pooled_sources[pooled_sources == 1].index.astype(str))
    if pooled_sources:
        pooled_lookup = flags.drop_duplicates("source_h5ad").set_index("source_h5ad")["sample_id"].astype(str).to_dict()
        pooled_mask = annot["source_h5ad"].astype(str).isin(pooled_sources)
        annot.loc[pooled_mask, "sample_id"] = annot.loc[pooled_mask, "source_h5ad"].map(pooled_lookup).astype(str)

    key = ["source_h5ad", "cell_barcode", "sample_id"]
    repl_cols = key + ["final_inclusion", "inclusion_reason", "source_qc_flags", "doublet_call", "filter_plan"]
    merged = annot.drop(columns=[c for c in ["final_inclusion"] if c in annot.columns]).merge(
        flags[repl_cols],
        on=key,
        how="left",
        validate="many_to_one",
    )
    assert len(merged) == len(annot), "row count changed after balanced inclusion join"
    if merged["final_inclusion"].isna().any():
        missing = int(merged["final_inclusion"].isna().sum())
        missing_by_source = (
            merged.loc[merged["final_inclusion"].isna(), "source_h5ad"]
            .value_counts()
            .head(20)
            .to_string()
        )
        raise RuntimeError(f"balanced final_inclusion did not map for all annotation rows: missing={missing}\n{missing_by_source}")

    rules = pd.read_csv(rules_path)
    rules_expanded = expand_rules(rules, created_at, mref)
    rules_expanded.to_csv(rules_path, index=False)

    mid_rules = rules_expanded[rules_expanded["target_level"] == "immune_mid_state"].copy()
    mid_key = mid_rules[["source_h5ad", "source_annotation_column", "source_label", "target_label"]].drop_duplicates()
    mid_key = mid_key.rename(columns={
        "source_annotation_column": "original_annotation_column",
        "source_label": "original_annotation",
        "target_label": "mid_from_rule",
    })
    mid_key["original_annotation"] = mid_key["original_annotation"].map(normalize_label)

    merged["original_annotation"] = merged["original_annotation"].map(normalize_label)
    merged = merged.merge(
        mid_key,
        on=["source_h5ad", "original_annotation_column", "original_annotation"],
        how="left",
    )
    text = merged["original_annotation"].astype(str)
    inferred = [
        infer_mid_from_text(major, label)
        for major, label in zip(merged["major_lineage"].astype(str), text)
    ]
    merged["immune_mid_state_previous"] = merged["immune_mid_state"]
    merged["immune_mid_state"] = merged["mid_from_rule"].fillna(pd.Series(inferred, index=merged.index))
    marker_mid = merged["annotation_rule_id"].astype(str).map(MARKER_RULE_TO_MID)
    marker_mid_mask = (
        merged["major_lineage"].isin(IMMUNE_LINEAGES)
        & (merged["immune_mid_state"] == "Immune_unspecified")
        & marker_mid.notna()
    )
    merged.loc[marker_mid_mask, "immune_mid_state"] = marker_mid[marker_mid_mask]
    merged.loc[~merged["major_lineage"].isin(IMMUNE_LINEAGES), "immune_mid_state"] = "NA"
    merged["functional_state_optional"] = [
        infer_functional(mid, label)
        for mid, label in zip(merged["immune_mid_state"].astype(str), text)
    ]
    merged = merged.drop(columns=["mid_from_rule"])

    merged["run_id"] = RUN_ID
    merged["created_at"] = created_at
    merged["input_manifest_ref"] = mref
    merged["step2_5_repair_branch"] = "qc0513_balanced"
    merged["annotation_repair_reason"] = "balanced_qc_inclusion_and_rule_traceability_repair_0513"

    pq.write_table(pa.Table.from_pandas(merged, preserve_index=False), annotation_path)
    merged.to_csv(out / "cell_state_annotation_v6_1.csv.gz", index=False, compression="gzip")

    marker = pd.read_csv(marker_path)
    marker["run_id"] = RUN_ID
    marker["created_at"] = created_at
    marker["input_manifest_ref"] = mref
    marker.to_csv(marker_path, index=False)

    write_summaries(merged, out, created_at, mref)
    write_manifest_and_report(merged, out, report_dir, created_at, mref, len(rules), len(rules_expanded))
    update_current_state(run_root, created_at)


def write_summaries(annot: pd.DataFrame, out: Path, created_at: str, mref: str) -> None:
    rows = []
    for cohort, grp in annot.groupby("cohort_id", dropna=False):
        main = grp[grp["final_inclusion"] == "include_main"]
        immune = main[main["major_lineage"].isin(IMMUNE_LINEAGES)]
        row = {
            "run_id": RUN_ID,
            "created_at": created_at,
            "input_manifest_ref": mref,
            "cohort_id": cohort,
            "n_cells_total": len(grp),
            "n_include_main": len(main),
            "n_include_sensitivity_only": int((grp["final_inclusion"] == "include_sensitivity_only").sum()),
            "n_excluded_qc_or_doublet": int(grp["final_inclusion"].astype(str).str.startswith("exclude").sum()),
            "n_immune_include_main": len(immune),
            "unknown_fraction": round(float((main["major_lineage"] == "Unknown").mean()) if len(main) else 0.0, 4),
            "annotation_conflict_fraction": round(float((main["annotation_conflict_flag"].fillna("") != "").mean()) if len(main) else 0.0, 4),
            "high_confidence_fraction": round(float((main["annotation_confidence"] == "high").mean()) if len(main) else 0.0, 4),
        }
        lineage_counts = main["major_lineage"].value_counts()
        for lin in ["T_NK", "Myeloid", "DC", "B", "Plasma", "Tumor_Epithelial", "Stromal_Fibroblast", "Endothelial", "Mast", "Cycling", "Erythroid", "Unknown"]:
            row[f"frac_{lin}"] = round(float(lineage_counts.get(lin, 0) / len(main)) if len(main) else 0.0, 4)
        mid_counts = immune["immune_mid_state"].value_counts()
        for mid in MID_STATES:
            row[f"immune_frac_{mid}"] = round(float(mid_counts.get(mid, 0) / len(immune)) if len(immune) else 0.0, 4)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "annotation_summary_by_cohort.csv", index=False)

    rows = []
    for (cohort, sample), grp in annot.groupby(["cohort_id", "sample_id"], dropna=False):
        main = grp[grp["final_inclusion"] == "include_main"]
        rows.append({
            "run_id": RUN_ID,
            "created_at": created_at,
            "input_manifest_ref": mref,
            "cohort_id": cohort,
            "sample_id": sample,
            "n_cells_total": len(grp),
            "n_include_main": len(main),
            "n_include_sensitivity_only": int((grp["final_inclusion"] == "include_sensitivity_only").sum()),
            "n_excluded_qc_or_doublet": int(grp["final_inclusion"].astype(str).str.startswith("exclude").sum()),
            "unknown_fraction": round(float((main["major_lineage"] == "Unknown").mean()) if len(main) else 0.0, 4),
            "annotation_conflict_fraction": round(float((main["annotation_conflict_flag"].fillna("") != "").mean()) if len(main) else 0.0, 4),
        })
    pd.DataFrame(rows).to_csv(out / "annotation_summary_by_sample.csv", index=False)


def write_manifest_and_report(annot: pd.DataFrame, out: Path, report_dir: Path, created_at: str, mref: str, n_rules_before: int, n_rules_after: int) -> None:
    main = annot[annot["final_inclusion"] == "include_main"]
    immune = main[main["major_lineage"].isin(IMMUNE_LINEAGES)]
    mid_cov = float((immune["immune_mid_state"] != "Immune_unspecified").mean()) if len(immune) else 0.0
    major_cov = float((main["major_lineage"] != "Unknown").mean()) if len(main) else 0.0
    status = "main_ready" if major_cov >= 0.90 and mid_cov >= 0.80 else ("degraded_ready" if major_cov >= 0.90 and mid_cov >= 0.60 else "blocked")

    cohort = pd.read_csv(out / "annotation_summary_by_cohort.csv")
    high_unknown = cohort[cohort["unknown_fraction"] > 0.20]
    high_unspec = cohort[cohort["immune_frac_Immune_unspecified"] > 0.30]
    high_conflict = cohort[cohort["annotation_conflict_fraction"] > 0.05]
    lines = [
        "# Annotation Conflict Report",
        "",
        f"Run ID: {RUN_ID}",
        f"Created: {created_at}",
        f"Repair branch: qc0513_balanced",
        f"Acceptance status: {status}",
        "",
        "## Coverage",
        f"- major_lineage coverage among include_main: {major_cov:.4f}",
        f"- immune_mid_state specified fraction among immune include_main: {mid_cov:.4f}",
        f"- include_main cells: {int((annot['final_inclusion'] == 'include_main').sum())}",
        f"- include_sensitivity_only cells: {int((annot['final_inclusion'] == 'include_sensitivity_only').sum())}",
        "",
        "## Cohorts with Unknown fraction > 20%",
    ]
    lines.extend([f"- {r.cohort_id}: {r.unknown_fraction*100:.1f}% unknown" for r in high_unknown.itertuples()] or ["- None"])
    lines.append("")
    lines.append("## Cohorts with major author-vs-marker conflicts > 5%")
    lines.extend([f"- {r.cohort_id}: {r.annotation_conflict_fraction*100:.1f}% conflict flag" for r in high_conflict.itertuples()] or ["- None"])
    lines.append("")
    lines.append("## Cohorts where immune cells cannot reach mid-state annotation > 30%")
    lines.extend([f"- {r.cohort_id}: {r.immune_frac_Immune_unspecified*100:.1f}% immune unspecified" for r in high_unspec.itertuples()] or ["- None"])
    lines.extend([
        "",
        "## Automated Annotation Status",
        "- CellTypist: not_run",
        "- Reason: optional and no local model was run; current status is degraded if mid-state coverage remains below the main-ready threshold.",
        "",
        "## Compliance Statement",
        "- Response labels were not read or used in this repair.",
        "- `qc0513_balanced` final inclusion was used as the Step2.5 inclusion source.",
        "- Main labels remain derived from author annotations, marker evidence, or conservative fallback only.",
    ])
    (out / "annotation_conflict_report.md").write_text("\n".join(lines), encoding="utf-8")

    manifest = {
        "run_id": RUN_ID,
        "step_id": "Step2.5_cell_annotation_harmonization",
        "created_at": created_at,
        "input_manifest_ref": mref,
        "repair_branch": "qc0513_balanced",
        "step2_5_status": status,
        "cell_count_total": int(len(annot)),
        "cell_count_include_main": int((annot["final_inclusion"] == "include_main").sum()),
        "cell_count_include_sensitivity_only": int((annot["final_inclusion"] == "include_sensitivity_only").sum()),
        "major_lineage_coverage_include_main": major_cov,
        "immune_mid_state_specified_fraction_immune_include_main": mid_cov,
        "mapping_rules_before_repair": int(n_rules_before),
        "mapping_rules_after_repair": int(n_rules_after),
        "automated_annotation_status": "not_run",
        "not_run_reason": "CellTypist was not run; author labels plus marker audit and explicit rule expansion were used.",
        "outputs": {
            "cell_state_annotation_v6_1.parquet": str(out / "cell_state_annotation_v6_1.parquet"),
            "cell_state_annotation_v6_1.csv.gz": str(out / "cell_state_annotation_v6_1.csv.gz"),
            "cell_state_mapping_rules_v6_1.csv": str(out / "cell_state_mapping_rules_v6_1.csv"),
            "marker_registry_v6_1.csv": str(out / "marker_registry_v6_1.csv"),
            "annotation_conflict_report.md": str(out / "annotation_conflict_report.md"),
            "annotation_summary_by_cohort.csv": str(out / "annotation_summary_by_cohort.csv"),
            "annotation_summary_by_sample.csv": str(out / "annotation_summary_by_sample.csv"),
        },
    }
    with (out / "annotation_manifest.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False)

    audit = [
        "# Step2.5 Balanced Repair Audit 0513",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        f"- status: `{status}`",
        f"- total rows: `{len(annot)}`",
        f"- include_main: `{int((annot['final_inclusion'] == 'include_main').sum())}`",
        f"- include_sensitivity_only: `{int((annot['final_inclusion'] == 'include_sensitivity_only').sum())}`",
        f"- major_lineage coverage: `{major_cov:.4f}`",
        f"- immune_mid_state specified fraction: `{mid_cov:.4f}`",
        f"- mapping rules before/after: `{n_rules_before}` -> `{n_rules_after}`",
        "",
        "Canonical Step2.5 now uses `03_qc/final_cell_inclusion_flags.qc0513_balanced.parquet` as its inclusion source.",
        "",
        "Response labels were not read or used.",
    ]
    (report_dir / "step2_5_balanced_repair_audit_0513.md").write_text("\n".join(audit), encoding="utf-8")


def update_current_state(run_root: Path, created_at: str) -> None:
    path = run_root / "STEP2_CURRENT_STATE.md"
    text = path.read_text(encoding="utf-8")
    text = text.replace("- current completed scope: `Step2.0` to `Step2.4`", "- current completed scope: `Step2.0` to `Step2.5`")
    text = text.replace("- `Step2.5+`: not started in this run directory", "- `Step2.5`: repaired with `qc0513_balanced` inclusion; see `10_reports/step2_5_balanced_repair_audit_0513.md`\n- `Step2.6+`: not started in this run directory")
    if "04_annotation/cell_state_annotation_v6_1.parquet" not in text:
        text = text.replace(
            "- `03_qc/step2_4_merge_audit.rerun_0511final.md`",
            "- `03_qc/step2_4_merge_audit.rerun_0511final.md`\n- `03_qc/final_cell_inclusion_flags.qc0513_balanced.parquet`\n- `04_annotation/cell_state_annotation_v6_1.parquet`\n- `04_annotation/annotation_manifest.yaml`\n- `04_annotation/annotation_conflict_report.md`",
        )
    note = f"\n## Step2.5 Repair Note\n\n- repaired_at: `{created_at}`\n- Step2.5 inclusion source: `03_qc/final_cell_inclusion_flags.qc0513_balanced.parquet`\n- Do not use pre-repair Step2.5 summaries.\n"
    if "## Step2.5 Repair Note" not in text:
        text += note
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
