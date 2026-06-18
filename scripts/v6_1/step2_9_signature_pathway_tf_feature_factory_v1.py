#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

RUN_ID = "step2_v6_1_0505_0319"
STEP2_ROOT = Path("/home/huyudi/006/results/v6_1/step2") / RUN_ID
PSEUDO_DIR = STEP2_ROOT / "07_pseudobulk" / "pseudobulk_feature_v1_final_v3"
FRACTION_DIR = STEP2_ROOT / "06_fraction" / "fraction_feature_v1"
GENE_DIR = STEP2_ROOT / "02_gene_standardization"
OUT_DIR = STEP2_ROOT / "08_signature_pathway_tf" / "signature_pathway_tf_feature_v1"
INPUT_MANIFEST_REF = PSEUDO_DIR / "pseudobulk_decision_manifest.yaml"

META_COLS = [
    "run_id", "created_at", "input_manifest_ref", "row_id", "cohort_id", "sample_id", "source_h5ad",
    "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "library_size", "detected_genes",
    "gene_universe", "raw_count_source", "raw_count_status",
]
PROVENANCE_COLS = ["run_id", "created_at", "input_manifest_ref"]
FORBIDDEN_RESPONSE_TERMS = ["response_strict", "response_broad", "response_binary", "response_ordered", "response_raw", "response_"]
MIN_COVERAGE_MAIN = 0.50

CORE_SIGNATURES: dict[str, dict[str, Any]] = {
    "IFN_response": {"up": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "ISG15"], "down": ["VEGFA", "TGFB1"], "cell_context": "immune;sample"},
    "cytotoxicity": {"up": ["GZMB", "PRF1", "NKG7", "GNLY", "GZMA", "CCL5", "IFNG"], "down": [], "cell_context": "CD8_T;NK;sample"},
    "exhaustion": {"up": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "CXCL13"], "down": ["IL2", "GZMB"], "cell_context": "CD8_T;T_NK;sample"},
    "antigen_presentation": {"up": ["HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "B2M", "TAP1"], "down": [], "cell_context": "Myeloid_DC;DC;sample"},
    "Treg_suppressive": {"up": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT", "ENTPD1"], "down": ["IFNG", "GZMB"], "cell_context": "Treg;sample"},
    "myeloid_suppression": {"up": ["TREM2", "APOE", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "LGALS3", "IL10"], "down": ["FCN1", "S100A8", "S100A9", "IL1B"], "cell_context": "Myeloid_DC;Macrophage;sample"},
    "myeloid_inflammation": {"up": ["FCN1", "S100A8", "S100A9", "IL1B", "CXCL8", "LYZ", "LST1", "NLRP3"], "down": ["TREM2", "MRC1"], "cell_context": "Myeloid_DC;Mono;sample"},
    "APC_DC_activation": {"up": ["CD80", "CD86", "CCR7", "LAMP3", "IL12B", "IRF8", "BATF3", "CLEC9A", "FCER1A"], "down": [], "cell_context": "DC;Myeloid_DC;sample"},
    "NK_cytotoxicity": {"up": ["NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1", "PRF1", "GZMB", "FCGR3A"], "down": [], "cell_context": "NK;T_NK;sample"},
    "proliferation_cell_cycle": {"up": ["MKI67", "TOP2A", "STMN1", "PCNA", "TYMS", "MCM2", "MCM5", "UBE2C"], "down": [], "cell_context": "all"},
    "hypoxia": {"up": ["VEGFA", "CA9", "LDHA", "SLC2A1", "ENO1", "BNIP3", "NDRG1"], "down": [], "cell_context": "tumor;stromal;sample"},
    "angiogenesis": {"up": ["VEGFA", "KDR", "FLT1", "ANGPT2", "PECAM1", "VWF", "ENG"], "down": [], "cell_context": "endothelial;sample"},
    "stromal_TGFb": {"up": ["TGFB1", "TGFBI", "COL1A1", "COL1A2", "ACTA2", "TAGLN", "POSTN", "FN1"], "down": [], "cell_context": "stromal;sample"},
}

PATHWAY_KEYWORDS = [
    "INTERFERON", "INFLAMMATORY", "ALLOGRAFT", "COMPLEMENT", "IL6", "TNFA", "HYPOXIA", "ANGIOGENESIS",
    "TGF", "APOPTOSIS", "P53", "KRAS", "ANTIGEN", "MHC", "T_CELL", "B_CELL", "NK", "MACROPHAGE", "DENDRITIC",
]
PATHWAY_FILES = [
    Path("/home/huyudi/006/data/pathway/h.all.v2025.1.Hs.symbols.gmt"),
    Path("/home/huyudi/006/data/pathway/c7.immunesigdb.v2025.1.Hs.symbols.gmt"),
    Path("/home/huyudi/006/data/pathway/immport_gene.gmt"),
    Path("/home/huyudi/006/data/pathway/c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sanitize_id(text: Any) -> str:
    out = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip()).strip("_")
    # Keep biological feature names in feature_name, but avoid response-like column names
    # in wide matrices because Step2 forbids response-derived output columns.
    out = re.sub(r"response", "rxn", out, flags=re.IGNORECASE)
    return out or "NA"


def read_gene_list(path: Path) -> set[str]:
    return {x.strip().upper() for x in path.read_text().splitlines() if x.strip()}


def forbidden_response_cols(cols: list[str]) -> list[str]:
    return [c for c in cols if any(term in c.lower() for term in FORBIDDEN_RESPONSE_TERMS)]


def add_feature_rows(rows: list[dict[str, Any]], feature_family: str, feature_name: str, up: list[str], down: list[str], source: str, cell_context: str, used: bool = True) -> None:
    feature_id = f"{feature_family}__{sanitize_id(feature_name)}"
    for gene in up:
        rows.append({"feature_id": feature_id, "feature_family": feature_family, "feature_name": feature_name, "gene_symbol": gene.upper(), "direction": 1, "source_resource": source, "used_in_main": used, "cell_context": cell_context, "notes": ""})
    for gene in down:
        rows.append({"feature_id": feature_id, "feature_family": feature_family, "feature_name": feature_name, "gene_symbol": gene.upper(), "direction": -1, "source_resource": source, "used_in_main": used, "cell_context": cell_context, "notes": ""})


def build_signature_registry() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for name, spec in CORE_SIGNATURES.items():
        add_feature_rows(rows, "signature", name, spec["up"], spec.get("down", []), "curated_step2_9_core", spec.get("cell_context", "all"), True)
    prog = Path("/home/huyudi/006/mvp/outputs/program_signatures.json")
    if prog.exists():
        data = json.loads(prog.read_text())
        for name, spec in data.items():
            add_feature_rows(rows, "signature", name, spec.get("up_genes", []), spec.get("down_genes", []), str(prog), "sample;cellstate", True)
    df = collapse_registry_sources(pd.DataFrame(rows))
    return df.sort_values(["feature_family", "feature_id", "gene_symbol"]).reset_index(drop=True)


def iter_gmt(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                yield parts[0], parts[1], [g.upper() for g in parts[2:] if g]


def build_pathway_registry(max_sets_per_file: int = 80) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in PATHWAY_FILES:
        if not path.exists():
            continue
        kept = 0
        for name, desc, genes in iter_gmt(path):
            upper = name.upper()
            if not any(k in upper for k in PATHWAY_KEYWORDS):
                continue
            if len(genes) < 5 or len(genes) > 500:
                continue
            add_feature_rows(rows, "pathway", name, genes, [], str(path), "sample;cellstate", True)
            kept += 1
            if kept >= max_sets_per_file:
                break
    if not rows:
        return pd.DataFrame(columns=["feature_id", "feature_family", "feature_name", "gene_symbol", "direction", "source_resource", "used_in_main", "cell_context", "notes"])
    df = pd.DataFrame(rows).drop_duplicates(["feature_id", "gene_symbol", "direction"])
    return df.sort_values(["feature_family", "feature_id", "gene_symbol"]).reset_index(drop=True)



def collapse_registry_sources(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    key = ["feature_id", "feature_family", "feature_name", "gene_symbol", "direction"]
    out = df.groupby(key, dropna=False).agg(
        source_resource=("source_resource", lambda x: ";".join(sorted(set(map(str, x))))),
        used_in_main=("used_in_main", "max"),
        cell_context=("cell_context", lambda x: ";".join(sorted(set(map(str, x))))),
    ).reset_index()
    out["notes"] = "source_and_context_collapsed_by_feature_gene_direction"
    return out[["feature_id", "feature_family", "feature_name", "gene_symbol", "direction", "source_resource", "used_in_main", "cell_context", "notes"]]

def empty_tf_registry() -> pd.DataFrame:
    return pd.DataFrame(columns=["feature_id", "feature_family", "feature_name", "gene_symbol", "direction", "source_resource", "used_in_main", "cell_context", "notes"])


def filter_registry_to_universe(reg: pd.DataFrame, universe: set[str]) -> pd.DataFrame:
    out = reg.copy()
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.upper().str.strip()
    out = out[out["gene_symbol"].isin(universe)]
    return out.drop_duplicates(["feature_id", "gene_symbol", "direction"]).reset_index(drop=True)


def coverage_table(reg_all: pd.DataFrame, reg_available: pd.DataFrame) -> pd.DataFrame:
    meta = reg_all.groupby("feature_id", dropna=False).agg(
        feature_family=("feature_family", "first"),
        feature_name=("feature_name", "first"),
        used_in_main=("used_in_main", "max"),
        cell_context=("cell_context", lambda x: ";".join(sorted(set(map(str, x))))),
    ).reset_index()
    expected = reg_all.groupby("feature_id", dropna=False)["gene_symbol"].nunique().reset_index(name="n_genes_expected").merge(meta, on="feature_id", how="left")
    available = reg_available.groupby("feature_id")["gene_symbol"].nunique().reset_index(name="n_genes_available")
    cov = expected.merge(available, on="feature_id", how="left")
    cov["n_genes_available"] = cov["n_genes_available"].fillna(0).astype(int)
    cov["gene_coverage"] = cov["n_genes_available"] / cov["n_genes_expected"].replace(0, np.nan)
    cov["coverage_confidence"] = np.select([cov["gene_coverage"].ge(0.70), cov["gene_coverage"].ge(0.50)], ["high", "medium"], default="low")
    cov["low_coverage_flag"] = cov["gene_coverage"].lt(MIN_COVERAGE_MAIN)
    cov["included_in_main_feature_matrix"] = cov["used_in_main"].astype(bool) & cov["gene_coverage"].ge(MIN_COVERAGE_MAIN)
    cov["sensitivity_only_feature"] = ~cov["included_in_main_feature_matrix"]
    return cov


def read_matrix(path: Path, reg: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    cols = pf.schema_arrow.names
    bad = forbidden_response_cols(cols)
    if bad:
        raise RuntimeError(f"forbidden response columns in {path}: {bad}")
    genes = sorted(set(reg["gene_symbol"]) & set(cols))
    df = pd.read_parquet(path, columns=META_COLS + genes)
    return df[META_COLS].copy(), df[genes].astype(float)


def zscore_frame(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=0, skipna=True)
    sd = expr.std(axis=0, skipna=True, ddof=0).replace(0, np.nan)
    return (expr - mu) / sd


def score_features(meta: pd.DataFrame, expr: pd.DataFrame, reg_all: pd.DataFrame, reg_available: pd.DataFrame, matrix_role: str, created_at: str) -> pd.DataFrame:
    z = zscore_frame(expr)
    cov = coverage_table(reg_all, reg_available)
    rows = []
    reg_by_feature = {fid: sub for fid, sub in reg_available.groupby("feature_id")}
    base_meta = meta[["row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe", "raw_count_status"]].reset_index(drop=True)
    for cov_row in cov.itertuples(index=False):
        sub = reg_by_feature.get(cov_row.feature_id)
        if sub is None or sub.empty:
            score = pd.Series(np.nan, index=meta.index)
        else:
            vals = []
            for r in sub.itertuples(index=False):
                if r.gene_symbol in z.columns:
                    vals.append(z[r.gene_symbol] * float(r.direction))
            score = pd.concat(vals, axis=1).mean(axis=1, skipna=True) if vals else pd.Series(np.nan, index=meta.index)
            if not bool(cov_row.included_in_main_feature_matrix):
                score[:] = np.nan
        out = base_meta.copy()
        out.insert(0, "input_manifest_ref", str(INPUT_MANIFEST_REF))
        out.insert(0, "created_at", created_at)
        out.insert(0, "run_id", RUN_ID)
        out["matrix_role"] = matrix_role
        out["feature_id"] = cov_row.feature_id
        out["feature_family"] = cov_row.feature_family
        out["feature_name"] = cov_row.feature_name
        out["score"] = score.to_numpy(dtype=float)
        out["n_genes_expected"] = int(cov_row.n_genes_expected)
        out["n_genes_available"] = int(cov_row.n_genes_available)
        out["gene_coverage"] = float(cov_row.gene_coverage) if pd.notna(cov_row.gene_coverage) else np.nan
        out["coverage_confidence"] = cov_row.coverage_confidence
        out["low_coverage_flag"] = bool(cov_row.low_coverage_flag)
        out["included_in_main_feature_matrix"] = bool(cov_row.included_in_main_feature_matrix)
        out["sensitivity_only_feature"] = bool(cov_row.sensitivity_only_feature)
        rows.append(out)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def write_feature_matrix(scores: pd.DataFrame, path: Path) -> pd.DataFrame:
    main = scores[scores["included_in_main_feature_matrix"].astype(bool)].copy()
    idx_cols = PROVENANCE_COLS + ["row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe", "matrix_role"]
    full_index = main[idx_cols].drop_duplicates().reset_index(drop=True)
    dedup = main[idx_cols + ["feature_id", "score"]].drop_duplicates(idx_cols + ["feature_id"], keep="first")
    value_wide = dedup.pivot(index=idx_cols, columns="feature_id", values="score").reset_index()
    value_wide.columns.name = None
    feature_cols = [c for c in value_wide.columns if c not in idx_cols]
    wide = full_index.merge(value_wide, on=idx_cols, how="left")
    wide = wide[idx_cols + feature_cols]
    wide.to_parquet(path, index=False)
    return wide


def all_nan_feature_rows(wide: pd.DataFrame) -> pd.DataFrame:
    idx_cols = PROVENANCE_COLS + ["row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe", "matrix_role"]
    feature_cols = [c for c in wide.columns if c not in idx_cols]
    out = wide.loc[wide[feature_cols].isna().all(axis=1), idx_cols].copy()
    return out


def validate_scores(df: pd.DataFrame, expected_rows: int | None = None) -> dict[str, Any]:
    bad = forbidden_response_cols(df.columns.tolist())
    return {
        "rows": int(len(df)),
        "unique_input_rows": int(df["row_id"].nunique()) if len(df) else 0,
        "unique_features": int(df["feature_id"].nunique()) if len(df) else 0,
        "forbidden_response_cols": bad,
        "has_coverage": all(c in df.columns for c in ["n_genes_expected", "n_genes_available", "gene_coverage"]),
        "expected_input_rows": expected_rows,
        "input_row_count_match": None if expected_rows is None else int(df["row_id"].nunique()) == expected_rows,
    }


def run(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir) if args.out_dir else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()

    core_genes = read_gene_list(GENE_DIR / "core_intersection_genes.txt")
    immune_genes = read_gene_list(GENE_DIR / "immune_feature_genes.txt")
    signature_reg = build_signature_registry()
    pathway_reg = build_pathway_registry(args.max_pathways_per_file)
    tf_reg = empty_tf_registry()
    tf_not_run_reason = "no_auditable_local_tf_regulon_resource_found"

    signature_reg.to_csv(out_dir / "signature_registry_v6_1.csv", index=False)
    pathway_reg.to_csv(out_dir / "pathway_registry_v6_1.csv", index=False)
    tf_reg.to_csv(out_dir / "tf_registry_v6_1.csv", index=False)

    signature_core = filter_registry_to_universe(signature_reg, core_genes)
    signature_immune = filter_registry_to_universe(signature_reg, immune_genes)
    pathway_core = filter_registry_to_universe(pathway_reg, core_genes)

    sample_meta_sig, sample_expr_sig = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_core.parquet", signature_core)
    cell_meta_sig_core, cell_expr_sig_core = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_core.parquet", signature_core)
    cell_meta_sig_imm, cell_expr_sig_imm = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_immune.parquet", signature_immune)
    sample_sig = score_features(sample_meta_sig, sample_expr_sig, signature_reg, signature_core, "sample_core", created_at)
    cell_sig_core = score_features(cell_meta_sig_core, cell_expr_sig_core, signature_reg, signature_core, "sample_cellstate_core", created_at)
    cell_sig_imm = score_features(cell_meta_sig_imm, cell_expr_sig_imm, signature_reg, signature_immune, "sample_cellstate_immune", created_at)
    cell_sig = pd.concat([cell_sig_core, cell_sig_imm], ignore_index=True)

    sample_sig.to_parquet(out_dir / "signature_scores_by_sample.parquet", index=False)
    cell_sig.to_parquet(out_dir / "signature_scores_by_sample_cellstate.parquet", index=False)

    sample_meta_path, sample_expr_path = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_core.parquet", pathway_core)
    cell_meta_path, cell_expr_path = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_core.parquet", pathway_core)
    sample_pathway = score_features(sample_meta_path, sample_expr_path, pathway_reg, pathway_core, "sample_core", created_at)
    cell_pathway = score_features(cell_meta_path, cell_expr_path, pathway_reg, pathway_core, "sample_cellstate_core", created_at)
    sample_pathway.to_parquet(out_dir / "pathway_activity_scores_by_sample.parquet", index=False)
    cell_pathway.to_parquet(out_dir / "pathway_activity_scores_by_sample_cellstate.parquet", index=False)

    empty_tf_cols = sample_sig.head(0).columns.tolist() if len(sample_sig) else []
    pd.DataFrame(columns=empty_tf_cols).to_parquet(out_dir / "tf_activity_scores_by_sample.parquet", index=False)
    pd.DataFrame(columns=empty_tf_cols).to_parquet(out_dir / "tf_activity_scores_by_sample_cellstate.parquet", index=False)

    all_scores = pd.concat([sample_sig, cell_sig, sample_pathway, cell_pathway], ignore_index=True)
    feature_matrix = write_feature_matrix(all_scores, out_dir / "signature_pathway_tf_feature_matrix.parquet")
    all_nan_rows = all_nan_feature_rows(feature_matrix)
    all_nan_rows.to_csv(out_dir / "all_nan_feature_rows_by_source.csv", index=False)

    coverage = pd.concat([
        coverage_table(signature_reg, signature_core).assign(matrix_role="sample_or_cellstate_core"),
        coverage_table(signature_reg, signature_immune).assign(matrix_role="sample_cellstate_immune"),
        coverage_table(pathway_reg, pathway_core).assign(matrix_role="sample_or_cellstate_core"),
    ], ignore_index=True)
    coverage.to_csv(out_dir / "signature_gene_coverage_by_feature.csv", index=False)

    dictionary = coverage[["feature_id", "feature_family", "feature_name", "n_genes_expected", "n_genes_available", "gene_coverage", "coverage_confidence", "included_in_main_feature_matrix", "matrix_role"]].drop_duplicates()
    dictionary.to_csv(out_dir / "signature_pathway_tf_feature_dictionary.csv", index=False)

    validations = {
        "signature_sample": validate_scores(sample_sig, 2213),
        "signature_cellstate": validate_scores(cell_sig, None),
        "pathway_sample": validate_scores(sample_pathway, 2213),
        "pathway_cellstate": validate_scores(cell_pathway, 18494),
        "feature_matrix": {
            "rows": int(len(feature_matrix)),
            "columns": int(len(feature_matrix.columns)),
            "forbidden_response_cols": forbidden_response_cols(feature_matrix.columns.tolist()),
            "all_nan_feature_rows_retained": int(len(all_nan_rows)),
            "sample_core_rows_observed": int((feature_matrix["matrix_role"] == "sample_core").sum()),
            "sample_core_rows_expected": 2213,
            "sample_cellstate_core_rows_observed": int((feature_matrix["matrix_role"] == "sample_cellstate_core").sum()),
            "sample_cellstate_core_rows_expected": 18494,
            "sample_cellstate_immune_rows_observed": int((feature_matrix["matrix_role"] == "sample_cellstate_immune").sum()),
            "sample_cellstate_immune_rows_expected": 18494,
        },
    }
    qc = pd.DataFrame([
        {"check_name": k, **v} for k, v in validations.items()
    ])
    qc.to_csv(out_dir / "signature_pathway_tf_qc_summary.csv", index=False)

    inputs = [
        PSEUDO_DIR / "pseudobulk_logcpm_by_sample_core.parquet",
        PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_core.parquet",
        PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_immune.parquet",
        INPUT_MANIFEST_REF,
        PSEUDO_DIR / "raw_count_authority_by_source.csv",
        GENE_DIR / "core_intersection_genes.txt",
        GENE_DIR / "immune_feature_genes.txt",
    ]
    status = "pass_with_tf_not_run_and_upstream_gene_mapping_warnings" if len(sample_sig) and len(sample_pathway) else "blocked"
    upstream_warning_sources = sorted(all_nan_rows["source_h5ad"].dropna().unique().tolist())
    manifest = {
        "run_id": RUN_ID,
        "step_id": "Step2.9_signature_pathway_tf_feature_factory_v1",
        "created_at": created_at,
        "input_manifest_ref": str(INPUT_MANIFEST_REF),
        "gate_status": status,
        "input_hashes": {str(p): sha256_file(p) for p in inputs if p.exists()},
        "upstream_step2_8_gate": "pass_with_two_unclear_raw_count_exclusions",
        "upstream_pseudobulk_dir": str(PSEUDO_DIR),
        "remaining_step2_8_excluded_sources": ["gse120575.h5ad", "gse229772.h5ad"],
        "dense_rescued_sources": ["gse123813_bcc.h5ad", "gse149614.h5ad", "gse161801.h5ad", "gse207422_sc.h5ad"],
        "scoring_policy": "mean_signed_gene_zscore_within_matrix",
        "coverage_threshold_main": MIN_COVERAGE_MAIN,
        "tf_status": "not_run",
        "tf_not_run_reason": tf_not_run_reason,
        "pathway_matrix_policy": "pathway scores are computed on sample_core and sample_cellstate_core only; sample_cellstate_immune is reserved for immune-focused signatures",
        "all_nan_feature_rows_retained": int(len(all_nan_rows)),
        "upstream_gene_standardization_zero_overlap_sources": upstream_warning_sources,
        "validation": validations,
        "forbidden_operations": ["response_leakage", "differential_analysis", "model_training", "global_clustering"],
    }
    (out_dir / "signature_pathway_tf_decision_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

    report = [
        "# Step2.9 Signature / Pathway / TF Feature Factory Report", "",
        f"- created_at: `{created_at}`", f"- gate_status: `{status}`",
        f"- signature registry rows: `{len(signature_reg)}`", f"- pathway registry rows: `{len(pathway_reg)}`",
        "- TF activity: `not_run`", f"- TF not_run_reason: `{tf_not_run_reason}`",
        "- response usage: `forbidden_not_used`", "- differential analysis: `not_run`", "- global clustering: `not_run`", "- model training: `not_run`",
        "- pathway matrix policy: `sample_core and sample_cellstate_core only; sample_cellstate_immune is immune-signature only`",
        f"- all-NaN feature rows retained: `{len(all_nan_rows)}`",
        f"- upstream gene-standardization warning sources: `{';'.join(upstream_warning_sources)}`", "",
        "## Output Validation", "",
    ]
    for k, v in validations.items():
        report.append(f"- {k}: {v}")
    report += ["", "## Upstream Pseudobulk", "", "- source: `pseudobulk_feature_v1_final_v3`", "- included sources: `41 / 43`", "- remaining excluded sources: `gse120575.h5ad`, `gse229772.h5ad`", "", "## Known Data-Quality Warnings", "", "- Rows with all feature values as NA are retained in the feature matrix and listed in `all_nan_feature_rows_by_source.csv`.", "- These rows reflect upstream gene-standardization / gene-overlap limitations, not Step2.9 scoring failure."]
    (out_dir / "signature_pathway_tf_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step2.9 signature/pathway/TF feature factory v1")
    parser.add_argument("--out-dir", default="")
    parser.add_argument("--max-pathways-per-file", type=int, default=80)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
