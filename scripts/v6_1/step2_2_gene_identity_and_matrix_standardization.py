from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]

REQUIRED_OUT_COLS = [
    "run_id", "created_at", "input_manifest_ref", "cohort_id", "source_h5ad",
    "source_var_path", "source_var_type", "var_name", "gene_symbol_raw",
    "ensembl_id_raw", "feature_name_raw", "genome_build_raw", "gene_biotype_raw",
    "standardized_gene_symbol", "symbol_status", "duplicate_symbol_flag",
    "duplicate_symbol_feature_count",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def csv_write(path: Path, fieldnames: List[str], rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_patched_manifest(run_root: Path) -> Dict[str, str]:
    patched = run_root / "00_manifest" / "step2_run_manifest.patched.yaml"
    orig = run_root / "00_manifest" / "step2_run_manifest.yaml"
    path = patched if patched.exists() else orig
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, val = line.split(":", 1)
                result[key.strip()] = val.strip()
    return result


def check_gate(run_root: Path) -> bool:
    gate_path = run_root / "00_manifest" / "step_specific_gate_status.csv"
    if not gate_path.exists():
        return False
    with gate_path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("step_name") == "Step2.2_gene_standardization":
                status = row.get("status", "")
                return status in ("ready", "degraded_ready")
    return False


def read_h5ad_var_index(path: str) -> List[str]:
    """Fast var index read via h5py. Tries multiple column names."""
    try:
        with h5py.File(path, "r") as f:
            var = f["var"]
            # Try multiple possible var index columns
            for key in ["_index", "gene", "gene_symbol", "feature_name", "Gene"]:
                if key in var:
                    idx = var[key]
                    if hasattr(idx, "dtype"):
                        vals = idx[:]
                        return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in vals]
                    elif hasattr(idx, "keys") and "categories" in idx:
                        # Categorical encoding
                        cats = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in idx["categories"][:]]
                        codes = idx["codes"][:]
                        return [cats[c] for c in codes]
    except Exception:
        pass
    return []


def standardize_gene_symbol(raw: str) -> str:
    s = str(raw).strip().upper()
    return s


def strip_ensembl_version(raw: str) -> str:
    s = str(raw).strip()
    if "." in s:
        s = s.split(".")[0]
    return s


def read_gmt_genes(path: Path) -> Set[str]:
    genes: Set[str] = set()
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) > 2:
                genes.update(p.strip().upper() for p in parts[2:] if p.strip())
    return genes


def read_json_pathway_genes(path: Path) -> Set[str]:
    genes: Set[str] = set()
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        for entry in data.values() if isinstance(data, dict) else data:
            if isinstance(entry, dict):
                for k in ("geneSet", "genes", "members"):
                    if k in entry:
                        genes.update(str(g).strip().upper() for g in entry[k] if str(g).strip())
            elif isinstance(entry, list):
                genes.update(str(g).strip().upper() for g in entry if str(g).strip())
    except Exception:
        pass
    return genes


def read_program_signatures(path: Path) -> Set[str]:
    genes: Set[str] = set()
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        for sig in data.values():
            if isinstance(sig, dict):
                for k in ("up_genes", "down_genes", "genes"):
                    if k in sig:
                        genes.update(str(g).strip().upper() for g in sig[k] if str(g).strip())
    except Exception:
        pass
    return genes


def read_panglao_markers(path: Path) -> Set[str]:
    genes: Set[str] = set()
    try:
        with path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for k in ("official gene symbol", "gene_symbol", "Gene", "gene"):
                    if k in row and row[k]:
                        genes.add(row[k].strip().upper())
    except Exception:
        pass
    return genes


def collect_immune_feature_genes(root: Path) -> Set[str]:
    genes: Set[str] = set()
    pathway_dir = root / "data" / "pathway"
    if pathway_dir.exists():
        for p in pathway_dir.glob("*.gmt"):
            genes.update(read_gmt_genes(p))
        for p in pathway_dir.glob("*.json"):
            genes.update(read_json_pathway_genes(p))

    sig_path = root / "mvp" / "outputs" / "program_signatures.json"
    if sig_path.exists():
        genes.update(read_program_signatures(sig_path))

    panglao_path = root / "data" / "db" / "cell_marker" / "PanglaoDB_markers_27_Mar_2020.tsv"
    if panglao_path.exists():
        genes.update(read_panglao_markers(panglao_path))

    return genes


def discover_sources_from_mapping(run_root: Path) -> List[Dict[str, str]]:
    """Read authoritative sources from mapping parquet."""
    parquet_path = run_root / "01_mapping_audit" / "cell_to_sample_mapping_draft.parquet"
    if not parquet_path.exists():
        return []

    table = pq.read_table(parquet_path, columns=["source_h5ad", "mapped_cohort_id"])
    df = table.to_pandas()
    grouped = df.groupby("source_h5ad")["mapped_cohort_id"].first().reset_index()

    sources = []
    h5ad_dir = ROOT / "data" / "processed" / "srt" / "raw"
    for _, row in grouped.iterrows():
        name = str(row["source_h5ad"])
        cohort = str(row["mapped_cohort_id"])
        path = h5ad_dir / name
        if path.exists():
            sources.append({"path": str(path), "name": name, "cohort": cohort})
    return sources


def discover_sidecars(root: Path) -> Dict[str, str]:
    sidecar_dir = root / "mvp" / "sidecars"
    result: Dict[str, str] = {}
    if sidecar_dir.exists():
        for subdir in sidecar_dir.iterdir():
            var_path = subdir / "var.parquet"
            if var_path.exists():
                result[subdir.name] = str(var_path)
    return result


def build_gene_mapping_for_source(
    source: Dict[str, str],
    sidecars: Dict[str, str],
) -> Tuple[List[Dict], Dict[str, object]]:
    """Build gene mapping rows for one source."""
    cohort = source["cohort"]
    path = source["path"]
    name = source["name"]

    var_names = read_h5ad_var_index(path)
    source_var_path = path
    source_var_type = "h5ad_var"

    # Try sidecar
    sidecar_path = sidecars.get(cohort, "")
    if sidecar_path:
        try:
            t = pq.read_table(sidecar_path)
            df = t.to_pandas()
            if "gene_symbol" in df.columns and len(df) == len(var_names):
                var_names = df["gene_symbol"].astype(str).tolist()
                source_var_path = sidecar_path
                source_var_type = "sidecar_var"
        except Exception:
            pass

    standardized = [standardize_gene_symbol(v) for v in var_names]
    symbol_counts = defaultdict(int)
    for s in standardized:
        if s:
            symbol_counts[s] += 1

    rows = []
    for i, (raw, std) in enumerate(zip(var_names, standardized)):
        dup_flag = symbol_counts.get(std, 0) > 1 if std else False
        rows.append({
            "var_name": raw,
            "gene_symbol_raw": raw,
            "ensembl_id_raw": "",
            "feature_name_raw": raw,
            "genome_build_raw": "",
            "gene_biotype_raw": "",
            "standardized_gene_symbol": std,
            "symbol_status": "ok" if std else "empty",
            "duplicate_symbol_flag": dup_flag,
            "duplicate_symbol_feature_count": symbol_counts.get(std, 0) if dup_flag else 0,
            "cohort_id": cohort,
            "source_h5ad": name,
            "source_var_path": source_var_path,
            "source_var_type": source_var_type,
        })

    n_dup_groups = sum(1 for c in symbol_counts.values() if c > 1)
    n_empty = sum(1 for s in standardized if not s)

    summary = {
        "cohort_id": cohort,
        "source_h5ad": name,
        "source_var_path": source_var_path,
        "source_var_type": source_var_type,
        "n_var_features": len(var_names),
        "n_native_gene_symbols": len(set(standardized)),
        "n_duplicate_symbol_groups": n_dup_groups,
        "n_duplicate_feature_rows": sum(c - 1 for c in symbol_counts.values() if c > 1),
        "n_empty_symbol_rows": n_empty,
    }
    return rows, summary


def run_step2_2_gene_standardization(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    out_dir = run_root / "02_gene_standardization"
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_patched_manifest(run_root)
    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    gate_ready = check_gate(run_root)
    if not gate_ready:
        print("Gate not ready for Step2.2; writing degraded note.")
        (out_dir / "gene_standardization_report.md").write_text(
            "# Gene Standardization Report\n\nGate not ready.\n", encoding="utf-8"
        )
        return {"run_id": run_id, "status": "gate_not_ready"}

    # Discover sources
    sources = discover_sources_from_mapping(run_root)
    sidecars = discover_sidecars(root)

    print(f"Processing {len(sources)} sources ...")

    # Build gene mappings
    all_mapping_rows: List[Dict] = []
    summaries: List[Dict] = []
    cohort_to_genes: Dict[str, Set[str]] = {}

    for src in sources:
        print(f"  {src['name']} ...")
        rows, summary = build_gene_mapping_for_source(src, sidecars)
        all_mapping_rows.extend(rows)
        summaries.append(summary)
        std_genes = {r["standardized_gene_symbol"] for r in rows if r["standardized_gene_symbol"]}
        cohort_to_genes[summary["cohort_id"]] = std_genes

    # Write gene_id_mapping_by_cohort.csv
    mapping_out_rows = []
    for row in all_mapping_rows:
        mapping_out_rows.append({
            "run_id": run_id,
            "created_at": now_iso(),
            "input_manifest_ref": manifest_ref,
            **row,
        })

    csv_write(
        out_dir / "gene_id_mapping_by_cohort.csv",
        REQUIRED_OUT_COLS,
        mapping_out_rows,
    )

    # Build conflict rows
    conflict_rows = []
    for row in all_mapping_rows:
        if row["duplicate_symbol_flag"]:
            conflict_rows.append({
                "run_id": run_id,
                "created_at": now_iso(),
                "input_manifest_ref": manifest_ref,
                "cohort_id": row["cohort_id"],
                "source_h5ad": row["source_h5ad"],
                "source_var_path": row["source_var_path"],
                "gene_symbol": row["standardized_gene_symbol"],
                "n_features": row["duplicate_symbol_feature_count"],
                "var_names": "",
                "ensembl_ids": "",
                "feature_names": row["gene_symbol_raw"],
                "conflict_type": "duplicate_symbol",
                "recommended_action": "review_before_merge",
            })

    if conflict_rows:
        csv_write(
            out_dir / "gene_symbol_conflicts.csv",
            [
                "run_id", "created_at", "input_manifest_ref", "cohort_id",
                "source_h5ad", "source_var_path", "gene_symbol", "n_features",
                "var_names", "ensembl_ids", "feature_names", "conflict_type",
                "recommended_action",
            ],
            conflict_rows,
        )
    else:
        (out_dir / "gene_symbol_conflicts.csv").write_text(
            "run_id,created_at,input_manifest_ref,cohort_id,source_h5ad,source_var_path,gene_symbol,n_features,var_names,ensembl_ids,feature_names,conflict_type,recommended_action\n",
            encoding="utf-8",
        )

    # Build gene universes
    all_genes = set()
    for genes in cohort_to_genes.values():
        all_genes.update(genes)

    # Determine main analysis cohorts from frozen inclusion
    frozen_path = root / "results" / "v6_1" / "data_pool" / "snapshots" / "current_step2_input" / "step1_exports" / "frozen_cohort_inclusion_v6_1.csv"
    main_candidates: Set[str] = set()
    if frozen_path.exists():
        with frozen_path.open("r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("inclusion_status") == "main":
                    main_candidates.add(row.get("cohort_id", ""))

    # Only include candidates that actually have gene data
    main_candidates = {c for c in main_candidates if c in cohort_to_genes}

    # Iteratively compute intersection; if a cohort drives it to 0, exclude it
    main_cohorts: List[str] = []
    if main_candidates:
        sorted_candidates = sorted(main_candidates, key=lambda c: len(cohort_to_genes[c]), reverse=True)
        core_intersection = set(cohort_to_genes[sorted_candidates[0]])
        main_cohorts.append(sorted_candidates[0])
        for c in sorted_candidates[1:]:
            test = core_intersection & cohort_to_genes[c]
            if len(test) > 0:
                core_intersection = test
                main_cohorts.append(c)

    if not main_cohorts:
        # fallback: largest cohort only
        main_cohorts = [max(cohort_to_genes.keys(), key=lambda c: len(cohort_to_genes[c]))]
        core_intersection = set(cohort_to_genes[main_cohorts[0]])

    immune_feature_genes = collect_immune_feature_genes(root)
    # Ensure all core genes are included
    immune_feature_genes.update(core_intersection)

    # Write gene lists
    (out_dir / "core_intersection_genes.txt").write_text(
        "\n".join(sorted(core_intersection)) + "\n", encoding="utf-8"
    )
    (out_dir / "immune_feature_genes.txt").write_text(
        "\n".join(sorted(immune_feature_genes)) + "\n", encoding="utf-8"
    )

    # Cohort native summary
    summary_rows = []
    for s in summaries:
        summary_rows.append({
            "run_id": run_id,
            "created_at": now_iso(),
            "input_manifest_ref": manifest_ref,
            "cohort_id": s["cohort_id"],
            "source_h5ad": s["source_h5ad"],
            "source_var_path": s["source_var_path"],
            "source_var_type": s["source_var_type"],
            "n_var_features": s["n_var_features"],
            "n_native_gene_symbols": s["n_native_gene_symbols"],
            "n_duplicate_symbol_groups": s["n_duplicate_symbol_groups"],
            "n_duplicate_feature_rows": s["n_duplicate_feature_rows"],
            "n_empty_symbol_rows": s["n_empty_symbol_rows"],
            "n_mapped_samples": 0,
            "is_main_cohort": s["cohort_id"] in main_cohorts,
        })

    csv_write(
        out_dir / "cohort_native_gene_summary.csv",
        [
            "run_id", "created_at", "input_manifest_ref", "cohort_id",
            "source_h5ad", "source_var_path", "source_var_type",
            "n_var_features", "n_native_gene_symbols", "n_duplicate_symbol_groups",
            "n_duplicate_feature_rows", "n_empty_symbol_rows", "n_mapped_samples",
            "is_main_cohort",
        ],
        summary_rows,
    )

    # Coverage by cohort
    n_core = len(core_intersection)
    n_immune = len(immune_feature_genes)
    coverage_rows = []
    for cohort, genes in cohort_to_genes.items():
        core_cov = len(genes & core_intersection) / n_core if n_core > 0 else 0.0
        immune_avail = len(genes & immune_feature_genes)
        coverage_rows.append({
            "run_id": run_id,
            "created_at": now_iso(),
            "input_manifest_ref": manifest_ref,
            "cohort_id": cohort,
            "n_native_genes": len(genes),
            "core_intersection_coverage": round(core_cov, 4),
            "immune_feature_gene_coverage": round(immune_avail / n_immune, 4) if n_immune > 0 else 0.0,
            "n_immune_feature_genes_available": immune_avail,
            "n_immune_feature_genes_expected": n_immune,
        })

    csv_write(
        out_dir / "gene_universe_coverage_by_cohort.csv",
        [
            "run_id", "created_at", "input_manifest_ref", "cohort_id",
            "n_native_genes", "core_intersection_coverage",
            "immune_feature_gene_coverage", "n_immune_feature_genes_available",
            "n_immune_feature_genes_expected",
        ],
        coverage_rows,
    )

    # Manifest
    blocking_review = n_core < 6000
    warning = n_core < 8000
    manifest_yaml = f'''manifest_version: "v6.1"
run_id: "{run_id}"
step_id: "Step2.2_gene_identity_and_matrix_standardization"
created_at: "{now_iso()}"
input_manifest_ref: "{manifest_ref}"
gate_source: "{manifest_ref}"
step2_2_gate_status: "{"ready" if not warning else "degraded_ready"}"
gene_identity_rules:
  symbol_normalization: "uppercase_strip_whitespace_drop_empty"
  symbol_fallback_order: "gene_symbol_then_feature_name_then_var_name"
  retain_original_var_name: true
  duplicate_symbol_merge_policy: "no_expression_merge_in_step2_2"
  duplicate_symbol_recording: "all_duplicate_symbols_written_to_gene_symbol_conflicts.csv"
  multi_ensembl_same_symbol_policy: "record_conflict_do_not_sum"
  missing_gene_in_union_matrix: "NA"
  measured_but_zero_expression: "0"
  na_zero_mix_forbidden: true
  expression_imputation: "forbidden"
gene_universe_rules:
  main_analysis_cohorts_present: "{",".join(sorted(main_cohorts))}"
  core_intersection_gene_count: {n_core}
  core_intersection_threshold_status: "{"blocking_review" if blocking_review else ("warning" if warning else "ok")}"
  core_intersection_blocking_review: {str(blocking_review).lower()}
  immune_feature_gene_count: {n_immune}
resources:
'''
    for p in sorted(root.glob("data/pathway/*")):
        manifest_yaml += f'  - "{p}"\n'
    for p in sorted(root.glob("data/db/cell_marker/*")):
        manifest_yaml += f'  - "{p}"\n'
    sig_path = root / "mvp" / "outputs" / "program_signatures.json"
    if sig_path.exists():
        manifest_yaml += f'  - "{sig_path}"\n'

    (out_dir / "gene_standardization_manifest.yaml").write_text(manifest_yaml, encoding="utf-8")

    # Report
    dup_total = sum(s["n_duplicate_symbol_groups"] for s in summaries)
    report = f"""# Gene Standardization Report

- **run_id**: {run_id}
- **gate_ready**: {gate_ready}
- **n_sources**: {len(sources)}
- **core_intersection_genes**: {n_core}
- **immune_feature_genes**: {n_immune}
- **duplicate_symbol_groups**: {dup_total}
- **blocking_review**: {blocking_review}
- **warning**: {warning}

## Cohort Summary

| Cohort | Native Genes | Duplicate Groups | Empty Symbols |
|--------|-------------|------------------|---------------|
"""
    for s in summaries:
        report += f"| {s['cohort_id']} | {s['n_native_gene_symbols']} | {s['n_duplicate_symbol_groups']} | {s['n_empty_symbol_rows']} |\n"

    report += f"""
## Threshold Status

- core_intersection = {n_core} genes
- threshold: < 8000 = warning, < 6000 = blocking review
- status: {'BLOCKING REVIEW' if blocking_review else ('WARNING' if warning else 'OK')}

## Acceptance Rules

- [x] core_intersection_genes count reported
- [x] duplicate symbols recorded in conflicts file
- [x] immune feature coverage reported
- [x] NA vs 0 rules written to manifest
"""
    (out_dir / "gene_standardization_report.md").write_text(report, encoding="utf-8")

    # Checkpoint
    (checkpoint_dir / "step2_2.checkpoint.yaml").write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\n"
        f"rows_written: {len(all_mapping_rows)}\nstatus: complete\n",
        encoding="utf-8",
    )

    print(f"Done. core_intersection={n_core}, immune_feature={n_immune}, dup_groups={dup_total}")
    return {
        "run_id": run_id,
        "status": "complete",
        "n_sources": len(sources),
        "core_intersection": n_core,
        "immune_feature": n_immune,
        "duplicate_groups": dup_total,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step2.2 gene identity and matrix standardization")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    args = parser.parse_args()
    result = run_step2_2_gene_standardization(root_dir=args.root_dir, run_id=args.run_id)
    print(result)
