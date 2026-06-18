#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from scipy import sparse

RUN_ID = "step2_v6_1_0505_0319"
STEP2_ROOT = Path("/home/huyudi/006/results/v6_1/step2") / RUN_ID
BASE_DIR = STEP2_ROOT / "08_signature_pathway_tf" / "signature_pathway_tf_feature_v1"
PSEUDO_DIR = STEP2_ROOT / "07_pseudobulk" / "pseudobulk_feature_v1_final_v3"
GENE_DIR = STEP2_ROOT / "02_gene_standardization"
REGISTRY_PATH = Path("/home/huyudi/006/data/tf_regulon/consensus_tf_regulon_v1/saezlab_tf_regulon_consensus_v1.csv.gz")
RESOURCE_MANIFEST = Path("/home/huyudi/006/data/tf_regulon/consensus_tf_regulon_v1/tf_regulon_resource_manifest.yaml")
OUT_DIR = STEP2_ROOT / "08_signature_pathway_tf" / "signature_pathway_tf_feature_v1_tf_rescue_v1"
INPUT_MANIFEST_REF = PSEUDO_DIR / "pseudobulk_decision_manifest.yaml"
META_COLS = [
    "run_id", "created_at", "input_manifest_ref", "row_id", "cohort_id", "sample_id", "source_h5ad",
    "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "library_size", "detected_genes",
    "gene_universe", "raw_count_source", "raw_count_status",
]
PROVENANCE_COLS = ["run_id", "created_at", "input_manifest_ref"]
IDX_COLS = PROVENANCE_COLS + ["row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe", "matrix_role"]
FORBIDDEN_RESPONSE_TERMS = ["response_strict", "response_broad", "response_binary", "response_ordered", "response_raw", "response_"]
MIN_COVERAGE = 0.50
MIN_TARGETS = 5


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def forbidden_response_cols(cols: list[str]) -> list[str]:
    return [c for c in cols if any(term in c.lower() for term in FORBIDDEN_RESPONSE_TERMS)]


def copy_base(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for p in BASE_DIR.iterdir():
        if p.is_file():
            shutil.copy2(p, out_dir / p.name)


def read_gene_list(path: Path) -> set[str]:
    return {x.strip().upper() for x in path.read_text().splitlines() if x.strip()}


def load_registry(path: Path, core_genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reg = pd.read_csv(path)
    reg["gene_symbol"] = reg["gene_symbol"].astype(str).str.upper().str.strip()
    reg["feature_id"] = reg["feature_id"].astype(str)
    reg["feature_family"] = "tf_activity"
    reg["feature_name"] = reg["feature_name"].astype(str)
    reg["direction"] = pd.to_numeric(reg["direction"], errors="coerce").fillna(1).astype(float)
    if "signed_weight" not in reg.columns:
        reg["signed_weight"] = reg["direction"]
    reg["signed_weight"] = pd.to_numeric(reg["signed_weight"], errors="coerce").fillna(reg["direction"])
    reg = reg.drop_duplicates(["feature_id", "gene_symbol", "direction", "signed_weight"])
    avail = reg[reg["gene_symbol"].isin(core_genes)].copy()
    meta = reg.groupby("feature_id", dropna=False).agg(
        feature_family=("feature_family", "first"),
        feature_name=("feature_name", "first"),
        n_genes_expected=("gene_symbol", "nunique"),
        used_in_main=("used_in_main", "max"),
        cell_context=("cell_context", lambda x: ";".join(sorted(set(map(str, x))))),
        resource=("resource", lambda x: ";".join(sorted(set(map(str, x))))),
    ).reset_index()
    got = avail.groupby("feature_id")["gene_symbol"].nunique().reset_index(name="n_genes_available")
    cov = meta.merge(got, on="feature_id", how="left")
    cov["n_genes_available"] = cov["n_genes_available"].fillna(0).astype(int)
    cov["gene_coverage"] = cov["n_genes_available"] / cov["n_genes_expected"].replace(0, np.nan)
    cov["coverage_confidence"] = np.select([cov["gene_coverage"].ge(0.70), cov["gene_coverage"].ge(0.50)], ["high", "medium"], default="low")
    cov["low_coverage_flag"] = cov["gene_coverage"].lt(MIN_COVERAGE) | cov["n_genes_available"].lt(MIN_TARGETS)
    cov["included_in_main_feature_matrix"] = cov["used_in_main"].astype(bool) & ~cov["low_coverage_flag"]
    cov["sensitivity_only_feature"] = ~cov["included_in_main_feature_matrix"]
    return reg, avail, cov


def read_matrix(path: Path, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    pf = pq.ParquetFile(path)
    cols = pf.schema_arrow.names
    bad = forbidden_response_cols(cols)
    if bad:
        raise RuntimeError(f"forbidden response columns in {path}: {bad}")
    genes = [g for g in genes if g in cols]
    df = pd.read_parquet(path, columns=META_COLS + genes)
    return df[META_COLS].copy(), df[genes].astype("float32")


def zscore(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=0, skipna=True)
    sd = expr.std(axis=0, skipna=True, ddof=0).replace(0, np.nan)
    return ((expr - mu) / sd).astype("float32")


def score_wide(meta: pd.DataFrame, expr: pd.DataFrame, avail: pd.DataFrame, cov: pd.DataFrame, matrix_role: str, created_at: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(expr.columns)
    feat = cov.sort_values("feature_id").reset_index(drop=True)
    gene_idx = {g: i for i, g in enumerate(genes)}
    feat_idx = {f: i for i, f in enumerate(feat["feature_id"])}
    sub = avail[avail["feature_id"].isin(feat_idx) & avail["gene_symbol"].isin(gene_idx)]
    rows = sub["gene_symbol"].map(gene_idx).to_numpy()
    cols = sub["feature_id"].map(feat_idx).to_numpy()
    vals = sub["signed_weight"].astype(float).to_numpy()
    mat = sparse.csr_matrix((vals, (rows, cols)), shape=(len(genes), len(feat)), dtype=np.float32)
    denom = np.asarray(np.abs(mat).sum(axis=0)).ravel().astype("float32")
    denom[denom == 0] = np.nan
    scores = zscore(expr).to_numpy(dtype=np.float32) @ mat
    scores = scores / denom
    main_mask = feat["included_in_main_feature_matrix"].astype(bool).to_numpy()
    scores[:, ~main_mask] = np.nan
    base = meta[["row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe"]].reset_index(drop=True)
    base.insert(0, "input_manifest_ref", str(INPUT_MANIFEST_REF))
    base.insert(0, "created_at", created_at)
    base.insert(0, "run_id", RUN_ID)
    base["matrix_role"] = matrix_role
    wide = pd.concat([base, pd.DataFrame(scores, columns=feat["feature_id"].tolist())], axis=1)
    return wide, feat


def write_long(wide: pd.DataFrame, feat: pd.DataFrame, path: Path, chunk_size: int = 50) -> int:
    writer = None
    total = 0
    meta = wide[IDX_COLS]
    features = feat["feature_id"].tolist()
    feat_meta = feat.set_index("feature_id")
    for start in range(0, len(features), chunk_size):
        cols = features[start:start + chunk_size]
        chunk = wide[IDX_COLS + cols].melt(id_vars=IDX_COLS, var_name="feature_id", value_name="score")
        add = feat_meta.loc[chunk["feature_id"], [
            "feature_family", "feature_name", "n_genes_expected", "n_genes_available", "gene_coverage",
            "coverage_confidence", "low_coverage_flag", "included_in_main_feature_matrix", "sensitivity_only_feature"
        ]].reset_index(drop=True)
        chunk = pd.concat([chunk.reset_index(drop=True), add], axis=1)
        table = pa.Table.from_pandas(chunk, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema, compression="snappy")
        writer.write_table(table)
        total += len(chunk)
    if writer is not None:
        writer.close()
    return total


def merge_feature_matrix(base_matrix_path: Path, sample_wide: pd.DataFrame, cell_wide: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    base = pd.read_parquet(base_matrix_path)
    tf_wide = pd.concat([sample_wide, cell_wide], ignore_index=True)
    tf_cols = [c for c in tf_wide.columns if c not in IDX_COLS]
    merged = base.merge(tf_wide[IDX_COLS + tf_cols], on=IDX_COLS, how="left")
    merged.to_parquet(out_path, index=False)
    return merged


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--tf-registry", default=str(REGISTRY_PATH))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    created_at = now_iso()
    copy_base(out_dir)
    core_genes = read_gene_list(GENE_DIR / "core_intersection_genes.txt")
    reg, avail, cov = load_registry(Path(args.tf_registry), core_genes)
    reg.to_csv(out_dir / "tf_registry_v6_1.csv", index=False)
    cov.assign(matrix_role="sample_or_cellstate_core").to_csv(out_dir / "tf_gene_coverage_by_feature.csv", index=False)
    target_genes = sorted(avail["gene_symbol"].unique())
    sample_meta, sample_expr = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_core.parquet", target_genes)
    cell_meta, cell_expr = read_matrix(PSEUDO_DIR / "pseudobulk_logcpm_by_sample_cellstate_core.parquet", target_genes)
    sample_wide, feat = score_wide(sample_meta, sample_expr, avail, cov, "sample_core", created_at)
    cell_wide, _ = score_wide(cell_meta, cell_expr, avail, cov, "sample_cellstate_core", created_at)
    sample_long_rows = write_long(sample_wide, feat, out_dir / "tf_activity_scores_by_sample.parquet")
    cell_long_rows = write_long(cell_wide, feat, out_dir / "tf_activity_scores_by_sample_cellstate.parquet")
    merged = merge_feature_matrix(out_dir / "signature_pathway_tf_feature_matrix.parquet", sample_wide, cell_wide, out_dir / "signature_pathway_tf_feature_matrix.parquet")
    all_nan = merged.loc[merged[[c for c in merged.columns if c not in IDX_COLS]].isna().all(axis=1), IDX_COLS].copy()
    all_nan.to_csv(out_dir / "all_nan_feature_rows_by_source.csv", index=False)
    base_cov = pd.read_csv(out_dir / "signature_gene_coverage_by_feature.csv")
    tf_cov = cov.assign(matrix_role="sample_or_cellstate_core")
    combined_cov = pd.concat([base_cov, tf_cov[base_cov.columns.intersection(tf_cov.columns).tolist()]], ignore_index=True, sort=False)
    combined_cov.to_csv(out_dir / "signature_gene_coverage_by_feature.csv", index=False)
    dictionary = combined_cov[["feature_id", "feature_family", "feature_name", "n_genes_expected", "n_genes_available", "gene_coverage", "coverage_confidence", "included_in_main_feature_matrix", "matrix_role"]].drop_duplicates()
    dictionary.to_csv(out_dir / "signature_pathway_tf_feature_dictionary.csv", index=False)
    qcs = []
    old_qc = pd.read_csv(out_dir / "signature_pathway_tf_qc_summary.csv")
    qcs.append(old_qc)
    qcs.append(pd.DataFrame([
        {"check_name": "tf_sample", "rows": sample_long_rows, "unique_input_rows": int(sample_wide["row_id"].nunique()), "unique_features": int(feat["feature_id"].nunique()), "forbidden_response_cols": forbidden_response_cols(sample_wide.columns.tolist()), "has_coverage": True, "expected_input_rows": 2213, "input_row_count_match": int(sample_wide["row_id"].nunique()) == 2213},
        {"check_name": "tf_cellstate", "rows": cell_long_rows, "unique_input_rows": int(cell_wide["row_id"].nunique()), "unique_features": int(feat["feature_id"].nunique()), "forbidden_response_cols": forbidden_response_cols(cell_wide.columns.tolist()), "has_coverage": True, "expected_input_rows": 18494, "input_row_count_match": int(cell_wide["row_id"].nunique()) == 18494},
        {"check_name": "feature_matrix_with_tf", "rows": int(len(merged)), "unique_input_rows": int(merged["row_id"].nunique()), "unique_features": int(len([c for c in merged.columns if c not in IDX_COLS])), "forbidden_response_cols": forbidden_response_cols(merged.columns.tolist()), "has_coverage": True, "expected_input_rows": None, "input_row_count_match": None},
    ]))
    pd.concat(qcs, ignore_index=True, sort=False).to_csv(out_dir / "signature_pathway_tf_qc_summary.csv", index=False)
    manifest_path = out_dir / "signature_pathway_tf_decision_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.update({
        "created_at": created_at,
        "gate_status": "pass_with_saezlab_tf_activity_and_upstream_gene_mapping_warnings",
        "tf_status": "run",
        "tf_method": "saezlab_omnipath_regulon_weighted_signed_mean_gene_zscore",
        "tf_resource_registry": str(Path(args.tf_registry)),
        "tf_resource_manifest": str(RESOURCE_MANIFEST),
        "tf_resource_registry_sha256": sha256_file(Path(args.tf_registry)),
        "tf_feature_count_total": int(feat["feature_id"].nunique()),
        "tf_feature_count_in_main": int(feat["included_in_main_feature_matrix"].sum()),
        "tf_target_genes_available": int(len(target_genes)),
        "tf_min_target_coverage_main": MIN_COVERAGE,
        "tf_min_targets_main": MIN_TARGETS,
        "trrust_included": False,
        "scenic_included": False,
        "feature_matrix_with_tf_shape": {"rows": int(len(merged)), "columns": int(len(merged.columns))},
    })
    manifest.pop("tf_not_run_reason", None)
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# Step2.9 Signature / Pathway / TF Feature Factory Report", "",
        f"- created_at: `{created_at}`",
        "- gate_status: `pass_with_saezlab_tf_activity_and_upstream_gene_mapping_warnings`",
        "- TF activity: `run`",
        "- TF resource: `saezlab/OmniPath CollecTRI + DoRothEA A/B/C`",
        "- TRRUST included: `false`",
        "- SCENIC included: `false`",
        f"- TF features total: `{feat['feature_id'].nunique()}`",
        f"- TF features in main matrix: `{int(feat['included_in_main_feature_matrix'].sum())}`",
        f"- TF target genes available in core matrix: `{len(target_genes)}`",
        f"- sample TF long rows: `{sample_long_rows}`",
        f"- cellstate TF long rows: `{cell_long_rows}`",
        f"- merged feature matrix rows: `{len(merged)}`",
        f"- merged feature matrix columns: `{len(merged.columns)}`",
        "- response usage: `forbidden_not_used`",
        "- differential analysis: `not_run`",
        "- global clustering: `not_run`",
        "- model training: `not_run`",
        "",
        "## Rescue Notes",
        "",
        "- This rescue branch starts from `signature_pathway_tf_feature_v1` and adds TF activity only.",
        "- Original signature and pathway score files are copied unchanged.",
        "- TF activity is computed on `sample_core` and `sample_cellstate_core`; `sample_cellstate_immune` receives no TF-specific columns.",
    ]
    report_text = "\n".join(report) + "\n"
    (out_dir / "tf_activity_rescue_report.md").write_text(report_text, encoding="utf-8")
    (out_dir / "signature_pathway_tf_report.md").write_text(report_text, encoding="utf-8")
    print("TF rescue complete")
    print(f"out_dir={out_dir}")
    print(f"tf_features={feat['feature_id'].nunique()} main={int(feat['included_in_main_feature_matrix'].sum())}")
    print(f"matrix_shape={merged.shape}")


if __name__ == "__main__":
    main()
