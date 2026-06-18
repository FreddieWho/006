#!/usr/bin/env python3
"""Phase6 response-blind module discovery.

Builds a unified expression foundation from raw-count h5ad objects, projects
support expression into the same gene space, and freezes one primary module
system without using response/outcome/split fields.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
import yaml

from phase4b_h5ad_expression_io import matrix_node, matrix_shape, read_axis_strings, read_matrix_chunk

ROOT = Path(__file__).resolve().parents[2]
PHASE35 = ROOT / "results/v6_2/phase3_5_single_cell_processing_qc_gate"
PHASE4A = ROOT / "results/v6_2/phase4a_cell_state_harmonization"
PHASE4B = ROOT / "results/v6_2/phase4b_immune_state_feature_construction"
PHASE5 = ROOT / "results/v6_2/phase5_strong_baseline_confounding_audit"
OUT = ROOT / "results/v6_2/phase6_response_blind_module_discovery"
RECHUNK_PATCH = ROOT / "results/v6_2/phase6_h5ad_rechunk_patch"

TODAY = str(date.today())
INPUT_VERSION = "frozen_v0"
PROCESSING_VERSION = "phase6_coarse_rescue_source_context_v1"
RANDOM_SEED = 1729
CHUNK_ROWS = 16384
CELLSTATE_MIN_CELLS = 20
SAMPLE_MIN_CELLS = 1
SMOKE_MAX_OBJECTS = 3
SMOKE_MAX_FIT_ROWS = 500
FULL_MAX_FIT_ROWS = 5000
NMF_RANKS = [8, 12, 16]
NMF_ITERS = 120
NMF_BOOT_ITERS = 60
SMOKE_BOOTSTRAPS = 5
FULL_BOOTSTRAPS = 20
TOP_GENES_PER_MODULE = 40
EPS = 1e-9

BANNED_TOKENS = [
    "response",
    "responder",
    "non_responder",
    "outcome",
    "survival",
    "progression",
    "recist",
    "mrecist",
    "trg",
    "death",
    "pfs",
    "os",
    "split",
    "train",
    "test",
    "label",
]
UNKNOWN_VALUES = {"", "unknown", "nan", "none", "na", "null", "low_quality_or_ambient"}
SOURCE_CONTEXT_OVERRIDES = {
    "GSE272734::object_1": {
        "harmonized_coarse_label": "T_NK",
        "harmonized_mid_label": "CD8_T",
        "reason": "published_sorted_cd8_t_cells_from_pbmc_anti_pd1_melanoma",
    }
}


def mkdirs() -> None:
    for sub in [
        "preflight",
        "gene_universe",
        "foundation/annotation_cache",
        "foundation/object_parts",
        "modules",
        "annotation",
        "audit",
        "handoff",
    ]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, **kwargs)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def write_md(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def load_yaml(path: Path) -> dict:
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def partial_sha256(path: Path, nbytes: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()


def safe_id(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")
    return token or "unknown"


def safe_file_id(value: object) -> str:
    return safe_id(value).replace("::", "__")


def norm_gene(value: object) -> str:
    return str(value).strip().upper()


ENSEMBL_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$", re.IGNORECASE)
GENE_SYMBOL_COLUMNS = ["gene_symbol", "symbol", "Gene", "gene_name", "feature_name", "gene", "features", "name"]
GENE_ALIAS_MAP = {
    "MALAT-1": "MALAT1",
    "PD-1": "PDCD1",
    "CD279": "PDCD1",
    "TIM3": "HAVCR2",
    "TIM-3": "HAVCR2",
    "CD366": "HAVCR2",
    "CD223": "LAG3",
    "CD152": "CTLA4",
    "CD28H": "TMIGD2",
    "B7H3": "CD276",
    "B7-H3": "CD276",
    "B7H4": "VTCN1",
    "B7-H4": "VTCN1",
    "CD11C": "ITGAX",
    "CD11B": "ITGAM",
    "CD14": "CD14",
    "CD16": "FCGR3A",
    "CD20": "MS4A1",
    "CD56": "NCAM1",
    "CD94": "KLRD1",
    "NKG2D": "KLRK1",
    "CD25": "IL2RA",
    "CD127": "IL7R",
    "KI67": "MKI67",
    "KRT8/18": "KRT8",
    "PANCK": "KRT8",
}
_ENSEMBL_MAP: dict[str, str] | None = None


def canonical_gene(value: object) -> str:
    gene = norm_gene(value).split(".")[0]
    return GENE_ALIAS_MAP.get(gene, gene)


def looks_like_ensembl(value: object) -> bool:
    return bool(ENSEMBL_RE.match(str(value).strip()))


def load_ensembl_to_symbol() -> dict[str, str]:
    candidates = [
        ROOT / "data/ref/gencode.v49.basic.annotation.gtf",
        ROOT / "data/ref/gencode.v49.annotation.gtf",
        ROOT / "data/db/tcga/gencode.v23.annotation.gtf",
    ]
    gtf = next((p for p in candidates if p.exists()), None)
    mapping: dict[str, str] = {}
    if gtf is None:
        return mapping
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gene_id = ""
            gene_name = ""
            for token in parts[8].split(";"):
                token = token.strip()
                if token.startswith("gene_id") and "\"" in token:
                    gene_id = token.split("\"", 2)[1].split(".")[0]
                elif token.startswith("gene_name") and "\"" in token:
                    gene_name = token.split("\"", 2)[1]
                if gene_id and gene_name:
                    mapping[gene_id] = canonical_gene(gene_name)
                    break
    return mapping


def ensembl_to_symbol_map() -> dict[str, str]:
    global _ENSEMBL_MAP
    if _ENSEMBL_MAP is None:
        _ENSEMBL_MAP = load_ensembl_to_symbol()
    return _ENSEMBL_MAP


def is_known(value: object) -> bool:
    return str(value).strip().lower() not in UNKNOWN_VALUES


def file_audit(paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for name, path in paths.items():
        exists = path.exists()
        rows.append(
            {
                "input_name": name,
                "path": str(path),
                "exists": "yes" if exists else "no",
                "readable": "yes" if exists and os.access(path, os.R_OK) else "no",
                "size_bytes": path.stat().st_size if exists and path.is_file() else "",
                "partial_sha256_4mb": partial_sha256(path) if exists and path.is_file() else "",
            }
        )
    return pd.DataFrame(rows)


def load_required_inputs() -> dict[str, Path]:
    paths = {
        "phase5_handoff": PHASE5 / "handoff/phase5_to_phase6_handoff.yaml",
        "matrix_layer_decision": PHASE35 / "matrix_layer_decision_table.frozen_v0.csv",
        "analysis_object_registry": PHASE35 / "analysis_object_registry.frozen_v0.csv",
        "cell_state_master": PHASE4A / "handoff/cell_state_annotation_master.csv.gz",
        "hvg_global": PHASE35 / "hvg_global_pan_cancer.tsv",
        "hvg_by_state": PHASE35 / "hvg_by_major_cell_state.tsv",
        "hvg_caution": PHASE35 / "hvg_caution_or_exclusion_list.tsv",
        "gene_overlap": PHASE35 / "gene_overlap_and_immune_core_coverage.csv",
        "program_signatures": ROOT / "mvp/outputs/program_signatures.json",
        "primary_matrix": PHASE4B / "matrix/immune_state_feature_matrix.primary.csv",
        "feature_dictionary": PHASE4B / "matrix/feature_dictionary_v6_2.csv",
        "normalized_support": PHASE4B / "pseudobulk/pseudobulk_matrix_normalized_expression.parquet",
        "resolution_eligibility": PHASE4B / "preflight/phase4b_resolution_aware_aggregation_eligibility.csv",
        "spatial_input_gate": ROOT / "results/v6_2/spatial_input_gate.csv",
        "phase5_leakage_audit": PHASE5 / "audit/phase5_label_leakage_audit.csv",
        "feature_family_eligibility": PHASE5 / "feature_family_downstream_eligibility.csv",
    }
    audit = file_audit(paths)
    write_csv(audit, OUT / "preflight/phase6_input_audit.csv")
    missing = audit[audit["exists"].ne("yes")]["input_name"].tolist()
    if missing:
        raise RuntimeError(f"Missing Phase6 required inputs: {missing}")
    handoff = load_yaml(paths["phase5_handoff"])
    if not handoff.get("verdict", "").startswith("CONDITIONAL_GO_TO_PHASE6"):
        raise RuntimeError(f"Phase5 handoff does not allow Phase6: {handoff.get('verdict')}")
    leakage = read_csv(paths["phase5_leakage_audit"])
    if "leakage_risk" in leakage and leakage["leakage_risk"].eq("yes").any():
        raise RuntimeError("Phase5 leakage audit has positive leakage risk.")
    return paths


def load_resolution_eligibility(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    cols = [
        "cohort_id",
        "sample_id",
        "harmonized_fine_label",
        "harmonized_mid_label",
        "harmonized_coarse_label",
        "eligible_coarse_pseudobulk",
        "eligible_mid_pseudobulk",
        "eligible_fine_pseudobulk",
        "resolution_repair_reason",
    ]
    df = pd.read_csv(path, dtype=str, keep_default_na=False, usecols=lambda c: c in cols)
    return df.drop_duplicates(
        ["cohort_id", "sample_id", "harmonized_fine_label", "harmonized_mid_label", "harmonized_coarse_label"],
        keep="first",
    )


def rebuild_annotation_object_ids() -> set[str]:
    ids = set(SOURCE_CONTEXT_OVERRIDES)
    qc_path = OUT / "foundation/object_join_qc.csv"
    if not qc_path.exists():
        return ids
    qc = read_csv(qc_path)
    if "status" not in qc or "object_id" not in qc:
        return ids
    version = qc["phase6_processing_version"] if "phase6_processing_version" in qc else pd.Series([""] * len(qc))
    needs = qc["status"].isin(["fallback_only_no_cellstate_units", "blocked_no_gene_overlap"]) & version.ne(PROCESSING_VERSION)
    ids.update(qc.loc[needs, "object_id"].astype(str).tolist())
    current = set()
    for obj in ids:
        path = annotation_cache_path(obj)
        if not path.exists():
            continue
        try:
            versions = pd.read_parquet(path, columns=["annotation_cache_version"])["annotation_cache_version"].astype(str)
            if not versions.empty and versions.eq(PROCESSING_VERSION).all():
                current.add(obj)
        except Exception:
            pass
    return ids - current


def apply_rechunk_manifest(raw: pd.DataFrame) -> pd.DataFrame:
    raw = raw.copy()
    raw["original_object_path"] = raw.get("object_path", "").astype(str)
    raw["original_selected_layer"] = raw.get("selected_layer", "").astype(str)
    raw["rechunk_applied"] = "no"
    raw["rechunk_manifest_path"] = ""
    raw["rechunk_validation_status"] = ""
    manifest_path = RECHUNK_PATCH / "rechunk_manifest.csv"
    if not manifest_path.exists():
        return raw
    manifest = read_csv(manifest_path)
    required = {
        "cohort_id",
        "object_id",
        "copy_status",
        "validation_status",
        "rechunked_h5ad",
        "phase6_selected_layer_after_rechunk",
    }
    if manifest.empty or not required.issubset(manifest.columns):
        return raw
    complete = manifest[
        manifest["copy_status"].eq("complete")
        & manifest["validation_status"].eq("pass")
        & manifest["rechunked_h5ad"].astype(str).map(lambda p: Path(p).exists())
    ].copy()
    if complete.empty:
        return raw
    complete = complete.drop_duplicates(["cohort_id", "object_id"], keep="last")
    keep = [
        "cohort_id",
        "object_id",
        "rechunked_h5ad",
        "phase6_selected_layer_after_rechunk",
        "validation_status",
    ]
    raw = raw.merge(complete[keep], on=["cohort_id", "object_id"], how="left")
    mask = raw["rechunked_h5ad"].fillna("").ne("")
    raw.loc[mask, "object_path"] = raw.loc[mask, "rechunked_h5ad"]
    raw.loc[mask, "selected_layer"] = raw.loc[mask, "phase6_selected_layer_after_rechunk"].replace("", "X")
    raw.loc[mask, "rechunk_applied"] = "yes"
    raw.loc[mask, "rechunk_manifest_path"] = str(manifest_path)
    raw.loc[mask, "rechunk_validation_status"] = raw.loc[mask, "validation_status"]
    return raw.drop(columns=["rechunked_h5ad", "phase6_selected_layer_after_rechunk", "validation_status"], errors="ignore")


def raw_object_plan(paths: dict[str, Path], mode: str) -> pd.DataFrame:
    layer = read_csv(paths["matrix_layer_decision"])
    registry = read_csv(paths["analysis_object_registry"])
    raw = layer[layer["raw_counts_allowed"].eq("yes")].copy()
    raw["registry_order"] = np.arange(len(raw))
    keep = [
        "cohort_id",
        "object_id",
        "object_path",
        "object_type",
        "object_read_status",
        "n_cells",
        "n_genes",
        "sample_column",
        "patient_column",
        "phase4b_eligibility",
        "downgrade_reason",
    ]
    raw = raw.merge(registry[[c for c in keep if c in registry.columns]], on=["cohort_id", "object_id"], how="left")
    raw["n_cells_num"] = pd.to_numeric(raw["n_cells"], errors="coerce").fillna(0).astype(int)
    raw["selected_for_phase6"] = np.where(
        raw["phase4b_eligibility"].eq("yes") & raw["object_read_status"].eq("readable"),
        "yes",
        "no",
    )
    raw["phase6_plan_status"] = np.where(raw["selected_for_phase6"].eq("yes"), "planned_raw_count_primary", "blocked")
    raw["phase6_plan_reason"] = np.where(raw["selected_for_phase6"].eq("yes"), "", "not_readable_or_not_phase4b_eligible")
    if mode == "smoke":
        selected_ids = set(
            raw[raw["selected_for_phase6"].eq("yes")]
            .sort_values("registry_order")
            .head(SMOKE_MAX_OBJECTS)["object_id"]
        )
        raw["run_in_current_mode"] = np.where(raw["object_id"].isin(selected_ids), "yes", "no")
    else:
        raw["run_in_current_mode"] = raw["selected_for_phase6"]
    raw = apply_rechunk_manifest(raw)
    raw = raw.sort_values(["run_in_current_mode", "registry_order"], ascending=[False, True]).reset_index(drop=True)
    write_csv(raw, OUT / "preflight/raw_count_object_plan.csv")
    return raw


def load_program_signature_genes(path: Path) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    data = json.loads(path.read_text())
    rows = []
    signatures: dict[str, set[str]] = {}
    for name, spec in data.items():
        if not isinstance(spec, dict):
            continue
        genes = {norm_gene(g) for key in ["up_genes", "down_genes"] for g in spec.get(key, []) if norm_gene(g)}
        signatures[name] = genes
        for gene in sorted(genes):
            rows.append({"gene": gene, "source": "program_signature", "source_detail": name})
    return pd.DataFrame(rows), signatures


def build_gene_universe(paths: dict[str, Path]) -> tuple[pd.DataFrame, list[str], dict[str, set[str]]]:
    rows = []
    global_hvg = pd.read_csv(paths["hvg_global"], sep="\t", dtype=str, keep_default_na=False)
    for _, row in global_hvg.iterrows():
        rows.append({"gene": norm_gene(row["gene"]), "source": "global_hvg", "source_detail": row.get("selection_rule", "")})
    by_state = pd.read_csv(paths["hvg_by_state"], sep="\t", dtype=str, keep_default_na=False)
    for _, row in by_state.iterrows():
        rows.append({"gene": norm_gene(row["gene"]), "source": "major_cell_state_hvg", "source_detail": row.get("major_cell_state", "")})
    sig_rows, signatures = load_program_signature_genes(paths["program_signatures"])
    rows.extend(sig_rows.to_dict("records"))
    raw = pd.DataFrame(rows)
    caution = pd.read_csv(paths["hvg_caution"], sep="\t", dtype=str, keep_default_na=False)
    caution_genes = {norm_gene(g) for g in caution["gene"]}
    raw = raw[raw["gene"].ne("")].copy()
    raw["excluded_by_caution"] = np.where(raw["gene"].isin(caution_genes), "yes", "no")
    write_csv(raw[raw["excluded_by_caution"].eq("yes")], OUT / "gene_universe/gene_universe_excluded_or_caution.csv")
    grouped = (
        raw[raw["excluded_by_caution"].eq("no")]
        .groupby("gene", as_index=False)
        .agg(
            sources=("source", lambda x: "|".join(sorted(set(x)))),
            source_details=("source_detail", lambda x: "|".join(sorted(set(str(v) for v in x if str(v))))),
        )
        .sort_values("gene")
    )
    grouped["input_version"] = INPUT_VERSION
    write_csv(grouped, OUT / "gene_universe/phase6_gene_universe.v0.csv")
    report = f"""# Phase6 Gene Universe QC Report

**Run date:** {TODAY}

- Included genes: {len(grouped)}
- Raw source rows: {len(raw)}
- Caution/excluded genes removed: {raw['excluded_by_caution'].eq('yes').sum()}
- Sources: global HVG, major-cell-state HVG, response-blind program signatures.
- Response, split, outcome fields were not used.
"""
    write_md(report, OUT / "gene_universe/gene_universe_qc_report.md")
    return grouped, grouped["gene"].tolist(), signatures


def annotation_cache_path(object_id: str) -> Path:
    return OUT / "foundation/annotation_cache" / f"{safe_file_id(object_id)}.parquet"


def prepare_annotation_cache(object_plan: pd.DataFrame, mode: str, resume: bool) -> pd.DataFrame:
    targets = object_plan[object_plan["run_in_current_mode"].eq("yes")].copy()
    all_target_ids = set(targets["object_id"])
    rebuild_ids = rebuild_annotation_object_ids()
    target_ids = all_target_ids if not resume else {obj for obj in all_target_ids if obj in rebuild_ids or not annotation_cache_path(obj).exists()}
    cache_rows = []
    for _, row in targets.iterrows():
        path = annotation_cache_path(row["object_id"])
        cache_rows.append(
            {
                "cohort_id": row["cohort_id"],
                "object_id": row["object_id"],
                "cache_path": str(path),
                "cache_exists_before": "yes" if path.exists() else "no",
                "cache_rebuild_requested": "yes" if row["object_id"] in target_ids else "no",
            }
        )
    if resume and not target_ids:
        out = pd.DataFrame(cache_rows)
        out["cache_status"] = "reused"
        write_csv(out, OUT / "foundation/annotation_cache_manifest.csv")
        return out

    cols = [
        "source_object_id",
        "cohort_id",
        "sample_id",
        "original_cell_barcode",
        "standardized_cell_barcode",
        "sample_key",
        "patient_key",
        "harmonized_coarse_label",
        "harmonized_mid_label",
        "harmonized_fine_label",
        "allowed_phase4b_pseudobulk",
        "allowed_downstream_use",
        "reliability_level",
        "low_quality_flag",
        "doublet_risk_flag",
        "ambient_rna_risk_flag",
        "stress_dissociation_flag",
    ]
    resolution = load_resolution_eligibility(PHASE4B / "preflight/phase4b_resolution_aware_aggregation_eligibility.csv")
    cache_dir = OUT / "foundation/annotation_cache"
    tmp_paths = {obj: cache_dir / f"{safe_file_id(obj)}.tmp.csv" for obj in target_ids}
    if not resume:
        for obj in target_ids:
            tmp_paths[obj].unlink(missing_ok=True)
            annotation_cache_path(obj).unlink(missing_ok=True)
    expected_counts = {
        row["object_id"]: int(row["n_cells_num"])
        for _, row in targets.iterrows()
    }
    seen_counts = {obj: 0 for obj in target_ids}
    for chunk in pd.read_csv(
        PHASE4A / "handoff/cell_state_annotation_master.csv.gz",
        dtype=str,
        keep_default_na=False,
        usecols=cols,
        chunksize=500_000,
    ):
        chunk = chunk[chunk["source_object_id"].isin(target_ids)].copy()
        if chunk.empty:
            continue
        if not resolution.empty:
            chunk = chunk.merge(
                resolution,
                on=["cohort_id", "sample_id", "harmonized_fine_label", "harmonized_mid_label", "harmonized_coarse_label"],
                how="left",
            )
        for obj, sub in chunk.groupby("source_object_id", sort=False):
            ann_part = transform_annotation_cache(sub.copy())
            tmp_path = tmp_paths[obj]
            ann_part.to_csv(
                tmp_path,
                mode="a",
                index=False,
                header=not tmp_path.exists(),
            )
            seen_counts[obj] += len(sub)
        if target_ids and all(seen_counts[obj] >= expected_counts.get(obj, 0) for obj in target_ids):
            break

    cache_out = []
    for obj in all_target_ids:
        path = annotation_cache_path(obj)
        tmp_path = tmp_paths.get(obj)
        if obj not in target_ids and path.exists():
            ann = pd.read_parquet(path, columns=["unit_use"])
            status = "reused"
            n_rows = len(ann)
            n_primary = int(ann["unit_use"].eq("primary_cellstate").sum())
        elif tmp_path is not None and tmp_path.exists():
            ann = pd.read_csv(tmp_path, dtype=str, keep_default_na=False)
            ann.to_parquet(path, index=False)
            tmp_path.unlink(missing_ok=True)
            status = "written"
            n_rows = len(ann)
            n_primary = int(ann["unit_use"].eq("primary_cellstate").sum())
        elif path.exists():
            status = "reused_empty_source"
            ann = pd.read_parquet(path)
            n_rows = len(ann)
            n_primary = int(ann.get("unit_use", pd.Series(dtype=str)).eq("primary_cellstate").sum())
        else:
            empty_cols = [
                "original_cell_barcode",
                "standardized_cell_barcode",
                "sample_key",
                "patient_key",
                "cell_state_level",
                "cell_state",
                "unit_use",
                "quality_flags",
            ]
            pd.DataFrame(columns=empty_cols).to_parquet(path, index=False)
            status = "empty_no_matching_cells"
            n_rows = 0
            n_primary = 0
        row = targets[targets["object_id"].eq(obj)].iloc[0]
        cache_out.append(
            {
                "cohort_id": row["cohort_id"],
                "object_id": obj,
                "cache_path": str(path),
                "cache_status": status,
                "n_annotation_rows": n_rows,
                "n_primary_cellstate_rows": n_primary,
            }
        )
    out = pd.DataFrame(cache_out)
    write_csv(out, OUT / "foundation/annotation_cache_manifest.csv")
    return out


def transform_annotation_cache(ann: pd.DataFrame) -> pd.DataFrame:
    ann["source_context_override"] = ""
    for object_id, spec in SOURCE_CONTEXT_OVERRIDES.items():
        mask = ann["source_object_id"].astype(str).eq(object_id)
        if mask.any():
            ann.loc[mask, "harmonized_coarse_label"] = spec["harmonized_coarse_label"]
            ann.loc[mask, "harmonized_mid_label"] = spec["harmonized_mid_label"]
            ann.loc[mask, "source_context_override"] = spec["reason"]
    mid_known = ann["harmonized_mid_label"].map(is_known)
    coarse_known = ann["harmonized_coarse_label"].map(is_known)
    for col in ["eligible_coarse_pseudobulk", "eligible_mid_pseudobulk", "eligible_fine_pseudobulk", "resolution_repair_reason"]:
        if col not in ann.columns:
            ann[col] = ""
    old_eligible = ann["allowed_phase4b_pseudobulk"].eq("yes") | ann["allowed_downstream_use"].isin(
        ["full", "mid_only", "coarse_or_support"]
    )
    resolution_eligible = ann[["eligible_coarse_pseudobulk", "eligible_mid_pseudobulk", "eligible_fine_pseudobulk"]].eq("yes").any(axis=1)
    source_context_eligible = ann["source_context_override"].ne("") & (mid_known | coarse_known)
    eligible = old_eligible | resolution_eligible | source_context_eligible
    quality_flags = []
    for _, row in ann.iterrows():
        flags = []
        for col in ["low_quality_flag", "doublet_risk_flag", "ambient_rna_risk_flag", "stress_dissociation_flag"]:
            if str(row.get(col, "")).lower() in {"1", "yes", "true"}:
                flags.append(col.replace("_flag", ""))
        quality_flags.append("|".join(flags))
    ann["cell_state_level"] = np.where(mid_known, "mid", np.where(coarse_known, "coarse", "unknown"))
    ann["cell_state"] = np.where(mid_known, ann["harmonized_mid_label"], np.where(coarse_known, ann["harmonized_coarse_label"], "Unknown"))
    ann["unit_use"] = np.where(eligible & ann["cell_state_level"].ne("unknown"), "primary_cellstate", "quality_flagged_pool")
    ann["quality_flags"] = quality_flags
    ann["annotation_cache_version"] = PROCESSING_VERSION
    keep = [
        "original_cell_barcode",
        "standardized_cell_barcode",
        "sample_key",
        "patient_key",
        "cell_state_level",
        "cell_state",
        "unit_use",
        "quality_flags",
        "allowed_phase4b_pseudobulk",
        "allowed_downstream_use",
        "reliability_level",
        "eligible_coarse_pseudobulk",
        "eligible_mid_pseudobulk",
        "eligible_fine_pseudobulk",
        "resolution_repair_reason",
        "source_context_override",
        "annotation_cache_version",
    ]
    return ann[keep].drop_duplicates(["original_cell_barcode"], keep="first")


def read_h5ad_var_column(handle: h5py.File, column: str) -> np.ndarray:
    return read_axis_strings(handle, "var", column).astype(str)


def resolve_var_gene_symbols(handle: h5py.File, target_genes: set[str]) -> tuple[np.ndarray, str, int]:
    var_names = read_h5ad_index(handle, "var").astype(str)
    candidates: list[tuple[str, np.ndarray]] = [("var_index", var_names)]
    for col in GENE_SYMBOL_COLUMNS:
        if col in handle["var"]:
            try:
                candidates.append((f"var/{col}", read_h5ad_var_column(handle, col)))
            except Exception:
                continue
    best_name = "var_index"
    best_values = var_names
    best_hits = -1
    for name, values in candidates:
        mapped = [canonical_gene(v) for v in values]
        hits = len(set(mapped) & target_genes)
        if hits > best_hits:
            best_name = name
            best_values = np.asarray(mapped, dtype=object)
            best_hits = hits
    if best_hits > 0:
        return best_values.astype(str), best_name, best_hits
    if any(looks_like_ensembl(v) for v in var_names[: min(1000, len(var_names))]):
        emap = ensembl_to_symbol_map()
        mapped = np.asarray([emap.get(str(v).split(".")[0], canonical_gene(v)) for v in var_names], dtype=object)
        hits = len(set(mapped) & target_genes)
        return mapped.astype(str), "ensembl_gtf_map", hits
    return np.asarray([canonical_gene(v) for v in var_names], dtype=object).astype(str), best_name, best_hits


def h5ad_gene_selection(handle: h5py.File, genes: list[str]) -> tuple[list[str], np.ndarray, dict[str, object]]:
    target_to_gene = {canonical_gene(g): g for g in genes}
    target_genes = set(target_to_gene)
    var_symbols, source, hits = resolve_var_gene_symbols(handle, target_genes)
    upper_to_idx: dict[str, int] = {}
    for idx, gene in enumerate(var_symbols):
        upper_to_idx.setdefault(canonical_gene(gene), idx)
    selected = sorted(
        [(target_to_gene[gene], upper_to_idx[gene]) for gene in target_genes if gene in upper_to_idx],
        key=lambda item: item[1],
    )
    selected_genes = [g for g, _ in selected]
    selected_cols = np.asarray([i for _, i in selected], dtype=int)
    return selected_genes, selected_cols, {
        "gene_match_source": source,
        "gene_match_hits": int(hits),
        "gene_alias_count": int(sum(1 for g in var_symbols if canonical_gene(g) != norm_gene(g))),
    }


def read_h5ad_index(handle: h5py.File, axis: str) -> np.ndarray:
    index_col = handle[axis].attrs.get("_index", "_index")
    if isinstance(index_col, bytes):
        index_col = index_col.decode()
    if index_col in handle[axis]:
        return read_axis_strings(handle, axis, str(index_col)).astype(str)
    if "_index" in handle[axis]:
        return read_axis_strings(handle, axis, "_index").astype(str)
    raise ValueError(f"Cannot resolve {axis} index column for h5ad object.")


def map_obs_to_annotation(obs_names: np.ndarray, ann: pd.DataFrame) -> pd.DataFrame:
    ann = ann.copy()
    ann["original_cell_barcode"] = ann["original_cell_barcode"].astype(str)
    ann["standardized_cell_barcode"] = ann["standardized_cell_barcode"].astype(str)
    original = ann.drop_duplicates("original_cell_barcode").set_index("original_cell_barcode")
    mapped = original.reindex(obs_names)
    missing = mapped["sample_key"].isna() if "sample_key" in mapped else pd.Series([True] * len(obs_names))
    if missing.any():
        standardized = ann.drop_duplicates("standardized_cell_barcode").set_index("standardized_cell_barcode")
        mapped_std = standardized.reindex(obs_names[missing.to_numpy()])
        mapped.loc[missing.to_numpy(), :] = mapped_std.to_numpy()
    mapped = mapped.reset_index(drop=True)
    return mapped.fillna("")


def hdf5_chunk_block_reason(node) -> str:
    if isinstance(node, h5py.Group) and {"data", "indices", "indptr"}.issubset(set(node.keys())):
        chunks = node["data"].chunks
        if chunks and int(chunks[0]) > 50_000_000:
            return f"blocked_hdf5_csr_data_chunk_too_large:{int(chunks[0])}"
    return ""


def registry_only_object(row: pd.Series, reason: str, n_genes: int, n_cells: int) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "expression_unit_id": f"{row['cohort_id']}::{row['object_id']}::registry_only",
                "cohort_id": row["cohort_id"],
                "object_id": row["object_id"],
                "sample_key": "",
                "patient_key": "",
                "cell_state_level": "object",
                "cell_state": "registry_only",
                "layer_family": "raw_count_direct_expression_blocked",
                "n_genes": n_genes,
                "library_size": "",
                "n_cells_used": n_cells,
                "inclusion_status": "registry_only",
                "downgrade_reason": reason,
            }
        ]
    )


def aggregate_object(row: pd.Series, genes: list[str], resume: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    object_id = row["object_id"]
    part_id = safe_file_id(object_id)
    cellstate_path = OUT / "foundation/object_parts" / f"{part_id}.cellstate.parquet"
    sample_path = OUT / "foundation/object_parts" / f"{part_id}.sample.parquet"
    registry_path = OUT / "foundation/object_parts" / f"{part_id}.registry.csv"
    qc_path = OUT / "foundation/object_parts" / f"{part_id}.qc.csv"
    if resume and cellstate_path.exists() and sample_path.exists() and registry_path.exists() and qc_path.exists():
        cached_qc = read_csv(qc_path)
        status = cached_qc.get("status", pd.Series(dtype=str)).astype(str)
        version = cached_qc.get("phase6_processing_version", pd.Series([""] * len(cached_qc))).astype(str)
        cached_path = cached_qc.get("object_path", pd.Series([""] * len(cached_qc))).astype(str)
        expected_path = str(Path(row["object_path"]))
        must_rebuild = (
            cached_qc.empty
            or status.str.startswith("blocked_hdf5", na=False).any()
            or status.isin(["blocked_no_gene_overlap", "fallback_only_no_cellstate_units"]).any()
            or version.ne(PROCESSING_VERSION).any()
            or cached_path.ne(expected_path).any()
        )
        if must_rebuild:
            pass
        else:
            return (
                pd.read_parquet(cellstate_path),
                pd.read_parquet(sample_path),
                read_csv(registry_path),
                cached_qc,
            )

    ann_path = annotation_cache_path(object_id)
    ann = pd.read_parquet(ann_path)
    object_path = Path(row["object_path"])
    cell_sums: dict[str, np.ndarray] = {}
    sample_sums: dict[str, np.ndarray] = {}
    cell_counts: dict[str, int] = defaultdict(int)
    sample_counts: dict[str, int] = defaultdict(int)
    cell_meta: dict[str, dict] = {}
    sample_meta: dict[str, dict] = {}
    quality_counter: Counter[str] = Counter()

    with h5py.File(object_path, "r") as handle:
        obs_names = read_h5ad_index(handle, "obs").astype(str)
        mapped = map_obs_to_annotation(obs_names, ann)
        selected_genes, selected_cols, gene_match = h5ad_gene_selection(handle, genes)
        if not selected_genes:
            qc = pd.DataFrame(
                [
                    {
                        "cohort_id": row["cohort_id"],
                        "object_id": object_id,
                        "object_path": str(object_path),
                        "n_obs": len(obs_names),
                        "n_joined": int(mapped["sample_key"].astype(str).ne("").sum()),
                        "join_rate": 0,
                        "n_selected_genes": 0,
                        "gene_match_source": gene_match["gene_match_source"],
                        "gene_match_hits": gene_match["gene_match_hits"],
                        "gene_alias_count": gene_match["gene_alias_count"],
                        "status": "blocked_no_gene_overlap",
                        "phase6_processing_version": PROCESSING_VERSION,
                    }
                ]
            )
            write_csv(qc, qc_path)
            empty = pd.DataFrame()
            empty.to_parquet(cellstate_path, index=False)
            empty.to_parquet(sample_path, index=False)
            write_csv(pd.DataFrame(), registry_path)
            return empty, empty, pd.DataFrame(), qc

        node = matrix_node(handle, row["selected_layer"])
        n_obs, n_vars = matrix_shape(node, len(obs_names), len(read_h5ad_index(handle, "var")))
        block_reason = hdf5_chunk_block_reason(node)
        if block_reason:
            qc = pd.DataFrame(
                [
                    {
                        "cohort_id": row["cohort_id"],
                        "object_id": object_id,
                        "object_path": str(object_path),
                        "n_obs": len(obs_names),
                        "n_joined": int(mapped["sample_key"].astype(str).ne("").sum()),
                        "join_rate": int(mapped["sample_key"].astype(str).ne("").sum()) / max(len(obs_names), 1),
                        "n_selected_genes": len(selected_genes),
                        "gene_match_source": gene_match["gene_match_source"],
                        "gene_match_hits": gene_match["gene_match_hits"],
                        "gene_alias_count": gene_match["gene_alias_count"],
                        "status": block_reason,
                        "phase6_processing_version": PROCESSING_VERSION,
                    }
                ]
            )
            write_csv(qc, qc_path)
            empty = pd.DataFrame()
            empty.to_parquet(cellstate_path, index=False)
            empty.to_parquet(sample_path, index=False)
            registry = registry_only_object(row, block_reason, len(selected_genes), len(obs_names))
            write_csv(registry, registry_path)
            return empty, empty, registry, qc
        sample_keys = mapped["sample_key"].astype(str).to_numpy()
        patient_keys = mapped["patient_key"].astype(str).to_numpy()
        cell_states = mapped["cell_state"].astype(str).to_numpy()
        cell_state_levels = mapped["cell_state_level"].astype(str).to_numpy()
        unit_use = mapped["unit_use"].astype(str).to_numpy()
        quality_flags = mapped["quality_flags"].astype(str).to_numpy()
        joined = sample_keys != ""
        for start in range(0, n_obs, CHUNK_ROWS):
            end = min(start + CHUNK_ROWS, n_obs)
            mat = read_matrix_chunk(node, start, end, selected_cols, n_vars)
            chunk_joined = joined[start:end]
            if not chunk_joined.any():
                continue
            # Sample fallback aggregates all joined cells, preserving data even when cell-state use is downgraded.
            sample_chunk_keys = sample_keys[start:end]
            for sample_key in np.unique(sample_chunk_keys[chunk_joined]):
                mask = sample_chunk_keys == sample_key
                summed = np.asarray(mat[mask].sum(axis=0)).ravel()
                sample_sums.setdefault(sample_key, np.zeros(len(selected_genes), dtype=np.float64))
                sample_sums[sample_key] += summed
                sample_counts[sample_key] += int(mask.sum())
                sample_meta.setdefault(
                    sample_key,
                    {
                        "cohort_id": row["cohort_id"],
                        "object_id": object_id,
                        "sample_key": sample_key,
                        "patient_key": str(patient_keys[start:end][mask][0]) if mask.any() else "",
                        "cell_state_level": "sample",
                        "cell_state": "sample_level",
                        "layer_family": "raw_count_sample_fallback",
                    },
                )
            # Primary cell-state aggregates only eligible, reliable cells.
            primary = chunk_joined & (unit_use[start:end] == "primary_cellstate")
            if primary.any():
                group_ids = np.asarray(
                    [
                        f"{row['cohort_id']}::{sample_keys[i]}::{cell_state_levels[i]}::{safe_id(cell_states[i])}::raw_count_cellstate"
                        for i in range(start, end)
                    ],
                    dtype=object,
                )
                for group_id in np.unique(group_ids[primary]):
                    rel = group_ids == group_id
                    summed = np.asarray(mat[rel].sum(axis=0)).ravel()
                    cell_sums.setdefault(group_id, np.zeros(len(selected_genes), dtype=np.float64))
                    cell_sums[group_id] += summed
                    cell_counts[group_id] += int(rel.sum())
                    first = np.where(rel)[0][0] + start
                    cell_meta.setdefault(
                        group_id,
                        {
                            "cohort_id": row["cohort_id"],
                            "object_id": object_id,
                            "sample_key": sample_keys[first],
                            "patient_key": patient_keys[first],
                            "cell_state_level": cell_state_levels[first],
                            "cell_state": cell_states[first],
                            "layer_family": "raw_count_cellstate",
                        },
                    )
            for flag in quality_flags[start:end][chunk_joined & ~primary]:
                quality_counter[flag or "unflagged_nonprimary"] += 1

    cell_df = sums_to_frame(cell_sums, cell_counts, cell_meta, selected_genes, min_cells=CELLSTATE_MIN_CELLS)
    sample_df = sums_to_frame(sample_sums, sample_counts, sample_meta, selected_genes, min_cells=SAMPLE_MIN_CELLS)
    registry = registry_from_parts(cell_df, sample_df, row, len(selected_genes), quality_counter)
    n_joined = int(joined.sum())
    qc = pd.DataFrame(
        [
            {
                "cohort_id": row["cohort_id"],
                "object_id": object_id,
                "object_path": str(object_path),
                "n_obs": int(len(obs_names)),
                "n_joined": n_joined,
                "join_rate": n_joined / len(obs_names) if len(obs_names) else 0,
                "n_selected_genes": int(len(selected_genes)),
                "gene_match_source": gene_match["gene_match_source"],
                "gene_match_hits": gene_match["gene_match_hits"],
                "gene_alias_count": gene_match["gene_alias_count"],
                "n_cellstate_units": int(len(cell_df)),
                "n_sample_units": int(len(sample_df)),
                "status": "complete" if len(cell_df) else "fallback_only_no_cellstate_units",
                "phase6_processing_version": PROCESSING_VERSION,
            }
        ]
    )
    cell_df.to_parquet(cellstate_path, index=False)
    sample_df.to_parquet(sample_path, index=False)
    write_csv(registry, registry_path)
    write_csv(qc, qc_path)
    return cell_df, sample_df, registry, qc


def sums_to_frame(
    sums: dict[str, np.ndarray],
    counts: dict[str, int],
    meta: dict[str, dict],
    genes: list[str],
    min_cells: int,
) -> pd.DataFrame:
    rows = []
    for unit_id, values in sums.items():
        m = dict(meta[unit_id])
        row = {"expression_unit_id": unit_id, **m}
        row["n_cells_used"] = int(counts.get(unit_id, 0))
        row.update({gene: float(value) for gene, value in zip(genes, values)})
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    gene_cols = [c for c in df.columns if c in genes]
    df["library_size"] = df[gene_cols].sum(axis=1)
    df = df[df["library_size"].gt(0) & df["n_cells_used"].ge(min_cells)].copy()
    return df


def registry_from_parts(
    cell_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    row: pd.Series,
    n_genes: int,
    quality_counter: Counter[str],
) -> pd.DataFrame:
    rows = []
    for frame, role in [(cell_df, "primary"), (sample_df, "fallback_projection")]:
        if frame.empty:
            continue
        for _, item in frame[["expression_unit_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family", "library_size", "n_cells_used"]].iterrows():
            rows.append(
                {
                    "expression_unit_id": item["expression_unit_id"],
                    "cohort_id": row["cohort_id"],
                    "object_id": row["object_id"],
                    "sample_key": item["sample_key"],
                    "patient_key": item["patient_key"],
                    "cell_state_level": item["cell_state_level"],
                    "cell_state": item["cell_state"],
                    "layer_family": item["layer_family"],
                    "n_genes": n_genes,
                    "n_cells_used": item["n_cells_used"],
                    "library_size": item["library_size"],
                    "inclusion_status": role,
                    "downgrade_reason": "" if role == "primary" else "sample_level_fallback_not_cellstate_primary",
                }
            )
    for flag, count in quality_counter.items():
        rows.append(
            {
                "expression_unit_id": f"{row['cohort_id']}::{row['object_id']}::quality_flagged::{safe_id(flag)}",
                "cohort_id": row["cohort_id"],
                "object_id": row["object_id"],
                "sample_key": "",
                "patient_key": "",
                "cell_state_level": "quality_flagged",
                "cell_state": flag,
                "layer_family": "quality_flagged_pool",
                "n_genes": n_genes,
                "library_size": "",
                "inclusion_status": "registry_only",
                "downgrade_reason": "not_used_for_primary_module_learning",
            }
        )
    return pd.DataFrame(rows)


def combine_foundation(parts: list[pd.DataFrame], out_path: Path) -> pd.DataFrame:
    if parts:
        df = pd.concat([p for p in parts if not p.empty], ignore_index=True, sort=False)
    else:
        df = pd.DataFrame()
    if not df.empty:
        df.to_parquet(out_path, index=False)
    else:
        pd.DataFrame().to_parquet(out_path, index=False)
    return df


def logcpm(count_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    if count_df.empty:
        return count_df.copy()
    meta_cols = [c for c in count_df.columns if c not in genes]
    present_genes = [g for g in genes if g in count_df.columns]
    x = count_df[present_genes].astype(float)
    lib = x.sum(axis=1).replace(0, np.nan)
    y = np.log1p(x.div(lib, axis=0).fillna(0) * 1_000_000)
    return pd.concat([count_df[meta_cols].reset_index(drop=True), y.reset_index(drop=True)], axis=1)


def build_normalized_support_projection(paths: dict[str, Path], genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    support = pd.read_parquet(paths["normalized_support"])
    id_cols = ["cohort_id", "sample_key"]
    gene_map = {norm_gene(c): c for c in support.columns if c not in id_cols}
    selected = [g for g in genes if g in gene_map]
    gene_df = pd.DataFrame(
        {gene: pd.to_numeric(support[gene_map[gene]], errors="coerce").fillna(0.0) for gene in selected},
        index=support.index,
    )
    out = pd.concat([support[id_cols].copy(), gene_df], axis=1)
    out["expression_unit_id"] = out["cohort_id"] + "::" + out["sample_key"] + "::sample::normalized_support_projection"
    out["patient_key"] = out["sample_key"]
    out["cell_state_level"] = "sample"
    out["cell_state"] = "sample_level"
    out["layer_family"] = "normalized_support_projection"
    meta_cols = ["expression_unit_id", "cohort_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]
    out = out[meta_cols + selected]
    registry = pd.DataFrame(
        {
            "expression_unit_id": out["expression_unit_id"],
            "cohort_id": out["cohort_id"],
            "object_id": "normalized_support_projection",
            "sample_key": out["sample_key"],
            "patient_key": out["patient_key"],
            "cell_state_level": out["cell_state_level"],
            "cell_state": out["cell_state"],
            "layer_family": out["layer_family"],
            "n_genes": len(selected),
            "library_size": out[selected].sum(axis=1) if selected else 0,
            "inclusion_status": "support_projection",
            "downgrade_reason": "normalized_support_not_raw_count_primary",
        }
    )
    out.to_parquet(OUT / "foundation/normalized_support_projection.v0.parquet", index=False)
    return out, registry


def build_foundation(paths: dict[str, Path], object_plan: pd.DataFrame, genes: list[str], mode: str, resume: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prepare_annotation_cache(object_plan, mode, resume)
    run_rows = object_plan[object_plan["run_in_current_mode"].eq("yes")].copy()
    cell_parts = []
    sample_parts = []
    registry_parts = []
    qc_parts = []
    for _, row in run_rows.iterrows():
        cell_df, sample_df, registry, qc = aggregate_object(row, genes, resume)
        cell_parts.append(cell_df)
        sample_parts.append(sample_df)
        registry_parts.append(registry)
        qc_parts.append(qc)
    cell_counts = combine_foundation(cell_parts, OUT / "foundation/raw_count_cellstate_pseudobulk.v0.parquet")
    sample_counts = combine_foundation(sample_parts, OUT / "foundation/raw_count_sample_pseudobulk.v0.parquet")
    cell_log = logcpm(cell_counts, genes)
    sample_log = logcpm(sample_counts, genes)
    cell_log.to_parquet(OUT / "foundation/logcpm_cellstate_pseudobulk.v0.parquet", index=False)
    sample_log.to_parquet(OUT / "foundation/logcpm_sample_pseudobulk.v0.parquet", index=False)
    support, support_registry = build_normalized_support_projection(paths, genes)
    registry = pd.concat(registry_parts + [support_registry], ignore_index=True, sort=False) if registry_parts else support_registry
    write_csv(registry, OUT / "foundation/unified_expression_unit_registry.csv")
    qc = pd.concat(qc_parts, ignore_index=True, sort=False) if qc_parts else pd.DataFrame()
    write_csv(qc, OUT / "foundation/object_join_qc.csv")
    blocked = object_plan[object_plan["run_in_current_mode"].ne("yes")].copy()
    blocked["phase6_block_reason"] = np.where(
        blocked["selected_for_phase6"].eq("yes"),
        f"not_run_in_{mode}_mode",
        blocked["phase6_plan_reason"],
    )
    write_csv(blocked, OUT / "foundation/pseudobulk_blocked_objects.csv")
    report = f"""# Phase6 Pseudobulk Foundation Report

**Run date:** {TODAY}

- Mode: `{mode}`
- Raw objects run: {len(run_rows)}
- Raw cell-state units: {len(cell_counts)}
- Raw sample fallback units: {len(sample_counts)}
- Normalized support projection units: {len(support)}
- Unified registry rows: {len(registry)}
- Objects with join QC rows: {len(qc)}

No response, outcome, or split fields are used in foundation construction.
"""
    write_md(report, OUT / "foundation/pseudobulk_foundation_report.md")
    return cell_log, support, registry


def numeric_matrix(df: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    present = [g for g in genes if g in df.columns]
    x = df[present].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float)
    nonzero = x.sum(axis=0) > 0
    present = [g for g, keep in zip(present, nonzero) if keep]
    x = x[:, nonzero]
    return df, x, present


def balanced_fit_rows(meta: pd.DataFrame, max_rows: int) -> np.ndarray:
    rng = np.random.default_rng(RANDOM_SEED)
    if len(meta) <= max_rows:
        return np.arange(len(meta))
    idx = []
    for _, group in meta.groupby(["cohort_id", "cell_state"], sort=False):
        take = max(1, int(max_rows * len(group) / len(meta)))
        values = group.index.to_numpy()
        idx.extend(rng.choice(values, size=min(take, len(values)), replace=False).tolist())
    idx = np.asarray(sorted(set(idx)), dtype=int)
    if len(idx) > max_rows:
        idx = rng.choice(idx, size=max_rows, replace=False)
    return np.sort(idx)


def nmf_mu(x: np.ndarray, rank: int, iters: int, seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    rng = np.random.default_rng(seed)
    x = np.maximum(x, 0.0)
    w = rng.random((x.shape[0], rank)) + 0.1
    h = rng.random((rank, x.shape[1])) + 0.1
    for _ in range(iters):
        h *= (w.T @ x) / (w.T @ w @ h + EPS)
        w *= (x @ h.T) / (w @ h @ h.T + EPS)
    err = float(np.linalg.norm(x - w @ h) / (np.linalg.norm(x) + EPS))
    return w, h, err


def top_gene_sets(h: np.ndarray, genes: list[str], top_n: int = TOP_GENES_PER_MODULE) -> list[set[str]]:
    sets = []
    for k in range(h.shape[0]):
        order = np.argsort(-h[k])[: min(top_n, h.shape[1])]
        sets.append({genes[i] for i in order})
    return sets


def max_jaccard(base: list[set[str]], other: list[set[str]]) -> list[float]:
    values = []
    for b in base:
        scores = []
        for o in other:
            union = len(b | o)
            scores.append(len(b & o) / union if union else 0.0)
        values.append(max(scores) if scores else 0.0)
    return values


def choose_rank(x_fit: np.ndarray, genes: list[str], mode: str) -> tuple[int, dict[int, dict], dict[int, np.ndarray]]:
    boot_n = SMOKE_BOOTSTRAPS if mode == "smoke" else FULL_BOOTSTRAPS
    metrics: dict[int, dict] = {}
    hs: dict[int, np.ndarray] = {}
    rng = np.random.default_rng(RANDOM_SEED)
    for rank in NMF_RANKS:
        _, h, err = nmf_mu(x_fit, rank, NMF_ITERS if mode == "full" else 80, RANDOM_SEED + rank)
        base_sets = top_gene_sets(h, genes)
        stabilities = []
        for b in range(boot_n):
            rows = rng.integers(0, x_fit.shape[0], size=x_fit.shape[0])
            _, hb, _ = nmf_mu(x_fit[rows], rank, NMF_BOOT_ITERS, RANDOM_SEED + rank * 100 + b)
            stabilities.extend(max_jaccard(base_sets, top_gene_sets(hb, genes)))
        metrics[rank] = {
            "rank": rank,
            "reconstruction_error": err,
            "median_bootstrap_jaccard": float(np.median(stabilities)) if stabilities else 0.0,
            "n_bootstrap_values": len(stabilities),
        }
        hs[rank] = h
    best = sorted(metrics.values(), key=lambda r: (-r["median_bootstrap_jaccard"], r["reconstruction_error"], r["rank"]))[0]["rank"]
    return int(best), metrics, hs


def score_modules(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    weights = h / (h.sum(axis=1, keepdims=True) + EPS)
    return x @ weights.T


def build_modules(cell_log: pd.DataFrame, support: pd.DataFrame, genes: list[str], signatures: dict[str, set[str]], mode: str) -> dict:
    if cell_log.empty:
        raise RuntimeError("No raw-count cell-state pseudobulk rows available for module discovery.")
    meta, x, present_genes = numeric_matrix(cell_log, genes)
    if x.shape[0] < 5 or x.shape[1] < 20:
        raise RuntimeError(f"Insufficient matrix for module discovery: rows={x.shape[0]}, genes={x.shape[1]}")
    fit_idx = balanced_fit_rows(meta, SMOKE_MAX_FIT_ROWS if mode == "smoke" else FULL_MAX_FIT_ROWS)
    x_fit = x[fit_idx]
    best_rank, rank_metrics, hs = choose_rank(x_fit, present_genes, mode)
    h = hs[best_rank]
    scores = score_modules(x, h)
    module_ids = [f"module_M{i+1:02d}" for i in range(h.shape[0])]

    score_df = meta[
        ["expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]
    ].copy()
    for i, module_id in enumerate(module_ids):
        score_df[module_id] = scores[:, i]
    write_csv(score_df, OUT / "modules/module_score_matrix.v0.csv")

    pt = score_df.copy()
    for module_id in module_ids:
        pt[module_id] = pd.to_numeric(pt[module_id], errors="coerce")
    pt = pt.groupby(["patient_key", "cohort_id"], as_index=False)[module_ids].mean()
    write_csv(pt, OUT / "modules/module_score_matrix.patient_timepoint.v0.csv")

    membership_rows = []
    for module_idx, module_id in enumerate(module_ids):
        weights = h[module_idx]
        max_weight = float(weights.max()) or 1.0
        for gene, weight in zip(present_genes, weights):
            if weight <= 0:
                continue
            membership_rows.append(
                {
                    "module_id": module_id,
                    "gene": gene,
                    "membership_weight": float(weight),
                    "relative_weight": float(weight / max_weight),
                    "is_top_gene": "yes" if gene in top_gene_sets(h[[module_idx]], present_genes, TOP_GENES_PER_MODULE)[0] else "no",
                    "module_system": f"unified_weighted_nmf_modules_v0_rank_{best_rank}",
                }
            )
    membership = pd.DataFrame(membership_rows)
    write_csv(membership, OUT / "modules/module_membership.v0.csv")

    stability = pd.DataFrame(rank_metrics.values()).sort_values("rank")
    stability["selected_for_freeze"] = np.where(stability["rank"].eq(best_rank), "yes", "no")
    write_csv(stability, OUT / "modules/module_stability_metrics.csv")

    support_projection = project_support_modules(support, present_genes, h, module_ids)
    write_csv(support_projection, OUT / "modules/module_projection_support_matrix.v0.csv")

    attributes, carrier, overlap = annotate_modules(score_df, membership, module_ids, signatures)
    write_csv(attributes, OUT / "modules/module_attribute_table.v0.csv")
    write_csv(carrier, OUT / "annotation/module_carrier_cellstate_table.csv")
    write_csv(overlap, OUT / "annotation/module_signature_overlap.csv")
    freeze_text = f"""# Phase6 Module Freeze Decision

**Run date:** {TODAY}

- Module system: `unified_weighted_nmf_modules_v0`
- Selected rank: {best_rank}
- Fit rows: {len(fit_idx)}
- Scored raw cell-state units: {len(score_df)}
- Support projection rows: {len(support_projection)}
- Selection rule: highest median bootstrap Jaccard, then lower reconstruction error, then lower rank.
- Response, outcome, and split fields were not used.
"""
    write_md(freeze_text, OUT / "modules/module_freeze_decision.md")
    stability_report = "# Phase6 Module Stability Report\n\n" + stability.to_string(index=False) + "\n"
    write_md(stability_report, OUT / "modules/module_stability_report.md")
    return {
        "best_rank": best_rank,
        "module_ids": module_ids,
        "score_df": score_df,
        "membership": membership,
        "stability": stability,
        "support_projection": support_projection,
        "attributes": attributes,
    }


def project_support_modules(support: pd.DataFrame, genes: list[str], h: np.ndarray, module_ids: list[str]) -> pd.DataFrame:
    if support.empty:
        return pd.DataFrame(columns=["expression_unit_id", "cohort_id", "sample_key"])
    meta_cols = ["expression_unit_id", "cohort_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]
    x = pd.DataFrame(
        {
            gene: (
                pd.to_numeric(support[gene], errors="coerce").fillna(0.0)
                if gene in support.columns
                else pd.Series(0.0, index=support.index)
            )
            for gene in genes
        },
        index=support.index,
    )
    # Cohort-wise z-score keeps normalized support separate from raw-count scale.
    z_parts = []
    for _, group in support.groupby("cohort_id", sort=False):
        vals = x.loc[group.index]
        sd = vals.std(axis=0, ddof=0).replace(0, 1.0)
        z_parts.append((vals - vals.mean(axis=0)) / sd)
    z = pd.concat(z_parts).loc[support.index].fillna(0.0).to_numpy(float)
    scores = score_modules(np.maximum(z - z.min(axis=0), 0.0), h)
    out = support[[c for c in meta_cols if c in support.columns]].copy()
    return pd.concat([out, pd.DataFrame(scores, columns=module_ids, index=support.index)], axis=1)


def annotate_modules(score_df: pd.DataFrame, membership: pd.DataFrame, module_ids: list[str], signatures: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    carrier_rows = []
    attr_rows = []
    overlap_rows = []
    for module_id in module_ids:
        tmp = score_df[["cell_state", module_id]].copy()
        tmp[module_id] = pd.to_numeric(tmp[module_id], errors="coerce")
        carrier = tmp.groupby("cell_state", as_index=False)[module_id].mean().sort_values(module_id, ascending=False)
        top_carrier = carrier.iloc[0]["cell_state"] if not carrier.empty else ""
        carrier["module_id"] = module_id
        carrier = carrier.rename(columns={module_id: "mean_module_score"})
        carrier_rows.extend(carrier[["module_id", "cell_state", "mean_module_score"]].to_dict("records"))

        genes = set(membership[(membership["module_id"].eq(module_id)) & (membership["is_top_gene"].eq("yes"))]["gene"])
        best_sig = ""
        best_j = 0.0
        for sig, sig_genes in signatures.items():
            union = genes | sig_genes
            j = len(genes & sig_genes) / len(union) if union else 0.0
            overlap_rows.append(
                {
                    "module_id": module_id,
                    "signature_name": sig,
                    "n_overlap_genes": len(genes & sig_genes),
                    "jaccard": j,
                    "overlap_genes": "|".join(sorted(genes & sig_genes)),
                }
            )
            if j > best_j:
                best_j = j
                best_sig = sig
        attr_rows.append(
            {
                "module_id": module_id,
                "module_system": "unified_weighted_nmf_modules_v0",
                "top_carrier_cell_state": top_carrier,
                "best_signature_overlap": best_sig,
                "best_signature_jaccard": best_j,
                "top_genes": "|".join(sorted(genes)),
                "allowed_downstream_use": "phase7_barrier_identifiability_candidate",
                "claim_boundary": "response_blind_module_only_no_response_direction_claim",
            }
        )
    return pd.DataFrame(attr_rows), pd.DataFrame(carrier_rows), pd.DataFrame(overlap_rows)


def mappability_report(attributes: pd.DataFrame, paths: dict[str, Path]) -> None:
    spatial = read_csv(paths["spatial_input_gate"])
    rows = []
    for _, row in attributes.iterrows():
        rows.append(
            {
                "module_id": row["module_id"],
                "top_carrier_cell_state": row["top_carrier_cell_state"],
                "spatial_localization_support_available": "yes" if spatial["usable_for_localization_support"].eq("yes").any() else "no",
                "spatial_response_claim_allowed": "no",
                "perturbation_coverage_status": "coverage_unknown_no_formal_phase6_target_prior_registry",
                "mappability_boundary": "support_annotation_only",
            }
        )
    report = pd.DataFrame(rows)
    write_csv(report, OUT / "annotation/module_mappability_table.csv")
    text = f"""# Phase6 Module Mappability Report

**Run date:** {TODAY}

- Spatial localization support datasets available: {int(spatial['usable_for_localization_support'].eq('yes').sum())}
- Spatial response/failure claim allowed in Phase6: no
- Perturbation coverage: unknown unless later formal target prior registry is frozen.
"""
    write_md(text, OUT / "annotation/module_mappability_report.md")


def leakage_audit(output_paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in output_paths:
        if not path.exists() or path.suffix not in {".csv"}:
            continue
        try:
            df = pd.read_csv(path, nrows=1)
        except Exception:
            continue
        for col in df.columns:
            lower = col.lower()
            tokens = [tok for tok in BANNED_TOKENS if tok in lower]
            rows.append(
                {
                    "file": str(path),
                    "column": col,
                    "leakage_tokens": "|".join(tokens),
                    "leakage_risk": "yes" if tokens and not col.startswith("claim_boundary") else "no",
                }
            )
    audit = pd.DataFrame(rows)
    write_csv(audit, OUT / "audit/phase6_leakage_audit.csv")
    return audit


def module_confounding(score_df: pd.DataFrame, module_ids: list[str]) -> pd.DataFrame:
    rows = []
    for module_id in module_ids:
        max_cohort = score_df.groupby("cohort_id")[module_id].count().max() / max(len(score_df), 1)
        rows.append(
            {
                "module_id": module_id,
                "n_expression_units": len(score_df),
                "n_cohorts": score_df["cohort_id"].nunique(),
                "n_cell_states": score_df["cell_state"].nunique(),
                "max_cohort_fraction": max_cohort,
                "cohort_dominance_flag": "yes" if max_cohort > 0.5 else "no",
                "allowed_primary_after_audit": "no" if max_cohort > 0.5 else "yes",
            }
        )
    out = pd.DataFrame(rows)
    write_csv(out, OUT / "audit/module_confounding_tags.csv")
    return out


def write_final_outputs(
    mode: str,
    foundation_registry: pd.DataFrame,
    module_result: dict,
    leakage: pd.DataFrame,
    confounding: pd.DataFrame,
) -> None:
    hard_blockers = []
    if leakage["leakage_risk"].eq("yes").any():
        hard_blockers.append("response_or_split_leakage_detected")
    stable = module_result["stability"]
    selected = stable[stable["selected_for_freeze"].eq("yes")]
    median_j = float(selected["median_bootstrap_jaccard"].iloc[0]) if not selected.empty else 0.0
    if len(module_result["module_ids"]) < 6:
        hard_blockers.append("too_few_modules")
    verdict = "BLOCKED" if hard_blockers else ("GO_TO_PHASE7" if median_j >= 0.45 else "CONDITIONAL_GO_TO_PHASE7")
    if mode == "smoke" and verdict == "GO_TO_PHASE7":
        verdict = "CONDITIONAL_GO_TO_PHASE7"
    object_qc_path = OUT / "foundation/object_join_qc.csv"
    object_qc = read_csv(object_qc_path) if object_qc_path.exists() else pd.DataFrame()
    status = object_qc["status"] if "status" in object_qc else pd.Series(dtype=str)
    manifest = {
        "phase": "phase6_response_blind_module_discovery",
        "created_at": TODAY,
        "mode": mode,
        "input_version": INPUT_VERSION,
        "verdict": verdict,
        "hard_blockers": hard_blockers,
        "module_system": "unified_weighted_nmf_modules_v0",
        "n_modules": len(module_result["module_ids"]),
        "selected_rank": module_result["best_rank"],
        "median_bootstrap_jaccard": median_j,
        "foundation": {
            "n_expression_units": int(len(foundation_registry)),
            "n_primary_units": int(foundation_registry["inclusion_status"].eq("primary").sum()) if not foundation_registry.empty else 0,
            "n_raw_sample_fallback_units": int(foundation_registry["inclusion_status"].eq("fallback_projection").sum()) if not foundation_registry.empty else 0,
            "n_support_projection_units": int(foundation_registry["inclusion_status"].eq("support_projection").sum()) if not foundation_registry.empty else 0,
            "n_registry_only_units": int(foundation_registry["inclusion_status"].eq("registry_only").sum()) if not foundation_registry.empty else 0,
        },
        "data_surface": {
            "raw_h5ad_objects_audited": int(len(object_qc)),
            "complete_raw_expression_objects": int(status.eq("complete").sum()),
            "fallback_only_objects": int(status.eq("fallback_only_no_cellstate_units").sum()),
            "hdf5_chunk_blocked_objects": int(status.str.startswith("blocked_hdf5", na=False).sum()),
            "no_gene_overlap_objects": int(status.eq("blocked_no_gene_overlap").sum()),
            "primary_cohorts": int(foundation_registry[foundation_registry["inclusion_status"].eq("primary")]["cohort_id"].nunique()) if not foundation_registry.empty else 0,
            "total_registry_cohorts": int(foundation_registry["cohort_id"].nunique()) if not foundation_registry.empty else 0,
        },
        "known_limitations": {
            "hdf5_csr_chunk_blocked_raw_expression_objects": int(status.str.startswith("blocked_hdf5", na=False).sum()),
            "module_stability_moderate": median_j < 0.45,
            "full_single_cell_integration_pending": True,
            "response_direction_not_allowed_in_phase6": True,
        },
        "primary_outputs": {
            "module_membership": str(OUT / "modules/module_membership.v0.csv"),
            "module_score_matrix": str(OUT / "modules/module_score_matrix.v0.csv"),
            "module_attribute_table": str(OUT / "modules/module_attribute_table.v0.csv"),
            "module_stability_metrics": str(OUT / "modules/module_stability_metrics.csv"),
            "unified_expression_unit_registry": str(OUT / "foundation/unified_expression_unit_registry.csv"),
            "object_join_qc": str(object_qc_path),
            "leakage_audit": str(OUT / "audit/phase6_leakage_audit.csv"),
        },
        "rules": {
            "response_blind_module_discovery": True,
            "single_primary_module_system": True,
            "normalized_support_as_projection_only": True,
            "spatial_response_claim_allowed": False,
            "perturbation_negative_evidence_allowed": False,
        },
    }
    write_yaml(manifest, OUT / "handoff/module_freeze_manifest.yaml")
    handoff = {
        "phase": "phase6_to_phase7_handoff",
        "verdict": verdict,
        "phase7_allowed_inputs": [
            str(OUT / "modules/module_membership.v0.csv"),
            str(OUT / "modules/module_score_matrix.v0.csv"),
            str(OUT / "modules/module_attribute_table.v0.csv"),
            str(OUT / "modules/module_stability_metrics.csv"),
            str(OUT / "audit/module_confounding_tags.csv"),
            str(OUT / "annotation/module_mappability_table.csv"),
        ],
        "phase7_blocked_scope": [
            "response_direction_claim",
            "spatial_response_failure_claim",
            "perturbation_negative_evidence_claim",
        ],
        "dominant_barrier_claim_allowed": False,
        "next_gate": "Phase7 barrier identifiability must decide which modules are separable.",
    }
    write_yaml(handoff, OUT / "handoff/phase6_to_phase7_handoff.yaml")
    rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            rows.append({"relative_path": str(path.relative_to(OUT)), "path": str(path), "size_bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(OUT / "handoff/phase6_output_index.tsv", sep="\t", index=False)
    report = f"""# Phase6 Response-blind Module Discovery Report

**Run date:** {TODAY}

## Verdict
`{verdict}`

## What Phase6 Built
- Unified expression unit registry rows: {len(foundation_registry)}
- Modules frozen: {len(module_result['module_ids'])}
- Selected rank: {module_result['best_rank']}
- Median bootstrap Jaccard: {median_j}
- Leakage hard blocker: {bool(hard_blockers)}

## Boundaries
- No response, outcome, or split field was used for module discovery.
- Normalized support expression was projected for support only.
- Spatial and perturbation are mappability annotations only.
- Dominant barrier claim is deferred to Phase7.
"""
    write_md(report, OUT / "PHASE6_RESPONSE_BLIND_MODULE_DISCOVERY_REPORT.md")


def run(args: argparse.Namespace) -> None:
    mkdirs()
    paths = load_required_inputs()
    object_plan = raw_object_plan(paths, args.mode)
    gene_universe, genes, signatures = build_gene_universe(paths)
    write_yaml(
        {
            "phase": "phase6_preflight",
            "created_at": TODAY,
            "mode": args.mode,
            "raw_objects_planned": int(object_plan["run_in_current_mode"].eq("yes").sum()),
            "gene_universe_size": int(len(genes)),
            "status": "PASS",
        },
        OUT / "preflight/phase6_preflight_status.yaml",
    )
    cell_log, support_projection, registry = build_foundation(paths, object_plan, genes, args.mode, args.resume)
    module_result = build_modules(cell_log, support_projection, genes, signatures, args.mode)
    mappability_report(module_result["attributes"], paths)
    confounding = module_confounding(module_result["score_df"], module_result["module_ids"])
    write_csv(registry, OUT / "audit/module_data_inclusion_audit.csv")
    write_csv(pd.DataFrame(columns=["module_id", "conflict_type", "conflict_description", "status"]), OUT / "audit/module_conflict_log.csv")
    leakage = leakage_audit(
        [
            OUT / "modules/module_membership.v0.csv",
            OUT / "modules/module_score_matrix.v0.csv",
            OUT / "modules/module_attribute_table.v0.csv",
            OUT / "foundation/unified_expression_unit_registry.csv",
        ]
    )
    write_final_outputs(args.mode, registry, module_result, leakage, confounding)
    print(f"Phase6 {args.mode} complete: modules={len(module_result['module_ids'])}, units={len(registry)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run v6.2 Phase6 response-blind module discovery.")
    parser.add_argument("--mode", choices=["smoke", "full"], default="smoke")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
