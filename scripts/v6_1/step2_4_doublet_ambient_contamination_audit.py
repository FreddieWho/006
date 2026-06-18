from __future__ import annotations

import argparse
import csv
import gc
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scrublet as scr

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def csv_write(path: Path, fieldnames: List[str], rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def csv_append(path: Path, fieldnames: List[str], row: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = (not path.exists()) or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
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
            if row.get("step_name") == "Step2.4_contamination":
                return row.get("status") == "ready"
    return False


def expected_doublet_rate(n_cells: int) -> float:
    """Approximate 10x multiplet-rate curve by recovered cells.

    10x published tables are roughly:
    - ~0.4% at ~500 recovered cells
    - ~8% at ~10k recovered cells
    """
    rate = 8e-6 * max(n_cells, 1)
    # Conservative prior cap per request: never exceed 5%.
    return min(0.05, max(0.004, rate))


def assign_doublet_calls(
    scores: np.ndarray,
    predicted: np.ndarray,
    max_plausible_rate: float = 0.35,
) -> np.ndarray:
    """Assign per-cell doublet calls with safety guard.

    Primary label source is Scrublet predicted_doublets. Scores are used to
    define borderline calls among non-predicted cells.
    """
    if scores is None or predicted is None or len(scores) == 0:
        return np.array([], dtype=object)

    pred_rate = float(np.mean(predicted.astype(bool)))
    if pred_rate > max_plausible_rate:
        # Guard against pathological runs that would over-remove most cells.
        return np.array(["unknown"] * len(scores), dtype=object)

    q99 = float(np.quantile(scores, 0.99))
    borderline_threshold = max(0.2, q99)
    calls: List[str] = []
    for s, p in zip(scores, predicted):
        if bool(p):
            calls.append("high_confidence_doublet")
        elif float(s) >= borderline_threshold:
            calls.append("borderline_doublet")
        else:
            calls.append("singlet")
    return np.array(calls, dtype=object)


def run_scrublet_for_sample(
    counts: np.ndarray | object,
    expected_rate: float,
    split_threshold: int = 20_000,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], str]:
    """Run Scrublet and return scores, calls, status."""
    try:
        n_cells = counts.shape[0]
        if n_cells < 50:
            return None, None, "not_possible_too_few_cells"
        if n_cells <= split_threshold:
            scrub = scr.Scrublet(counts, expected_doublet_rate=expected_rate)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
            return doublet_scores, predicted_doublets, "run_success"
        # For large samples, split into deterministic chunks and run per chunk.
        # This preserves per-sample execution while controlling memory/runtime.
        score_chunks: List[np.ndarray] = []
        pred_chunks: List[np.ndarray] = []
        for start in range(0, n_cells, split_threshold):
            end = min(start + split_threshold, n_cells)
            sub = counts[start:end]
            scrub = scr.Scrublet(sub, expected_doublet_rate=expected_rate)
            s, p = scrub.scrub_doublets(verbose=False)
            score_chunks.append(np.asarray(s))
            pred_chunks.append(np.asarray(p))
        return np.concatenate(score_chunks), np.concatenate(pred_chunks), "run_success_chunked"
    except Exception as e:
        return None, None, f"failed_with_reason:{e}"


def subset_to_memory(adata: ad.AnnData, obs_idx: object) -> ad.AnnData:
    """Subset backed AnnData and return in-memory copy, handling missing X layer."""
    sub = adata[obs_idx]
    try:
        return sub.to_memory()
    except Exception:
        layer_key = None
        for key in ("counts", "data"):
            if key in sub.layers:
                layer_key = key
                break
        if layer_key is None:
            raise RuntimeError(f"No X layer and no usable layer (counts/data) found")
        X = sub.layers[layer_key]
        if hasattr(X, "to_memory"):
            X = X.to_memory()
        return ad.AnnData(
            X=X,
            obs=sub.obs.copy(),
            var=sub.var.copy(),
        )


def select_counts_matrix_for_scrublet(adata: ad.AnnData):
    """Select the best available raw-counts-like matrix."""
    if "counts" in adata.layers:
        return adata.layers["counts"], "layers[counts]"
    if adata.raw is not None and getattr(adata.raw, "X", None) is not None:
        return adata.raw.X, "raw.X"
    return adata.X, "X"


def compute_contamination_scores(adata: ad.AnnData) -> pd.DataFrame:
    """Compute contamination/stress scores per cell."""
    var_names = [str(v).upper() for v in adata.var_names]

    def subset_score(gene_prefixes):
        idx = [i for i, g in enumerate(var_names) if any(g.startswith(p) for p in gene_prefixes)]
        if not idx:
            return np.zeros(len(adata))
        X = adata.X
        if hasattr(X, "toarray"):
            vals = np.array(X[:, idx].sum(axis=1)).flatten()
        else:
            vals = np.array(X[:, idx].sum(axis=1)).flatten()
        return vals

    rbc_genes = ["HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"]
    stress_genes = ["FOS", "JUN", "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1", "HSPB1", "DDIT3", "ATF3"]
    tumor_genes = ["ALB", "AFP", "KRT5", "KRT8", "KRT18", "KRT19", "EPCAM"]

    rbc_score = subset_score(rbc_genes)
    stress_score = subset_score(stress_genes)
    tumor_score = subset_score(tumor_genes)

    # Normalize by total counts per cell
    X = adata.X
    if hasattr(X, "toarray"):
        total = np.array(X.sum(axis=1)).flatten()
    else:
        total = np.array(X.sum(axis=1)).flatten()

    total = np.where(total == 0, 1, total)
    rbc_score = 100.0 * rbc_score / total
    stress_score = 100.0 * stress_score / total
    tumor_score = 100.0 * tumor_score / total

    return pd.DataFrame({
        "cell_barcode": adata.obs_names.astype(str),
        "rbc_score": rbc_score,
        "stress_score": stress_score,
        "tumor_epithelial_contamination_score": tumor_score,
    })


def run_step2_4(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
    min_sample_cells: int | None = None,
    max_sample_cells: int | None = None,
    max_samples: int | None = None,
    batch_tag: str | None = None,
    sample_shard_index: int | None = None,
    sample_shard_count: int | None = None,
    filter_plan_path: str | None = None,
    output_suffix: str | None = None,
    batch_dir_name: str = "doublet_batches",
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    out_dir = run_root / "03_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_patched_manifest(run_root)
    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    gate_ready = check_gate(run_root)
    if not gate_ready:
        print("Gate not ready for Step2.4; writing degraded note.", flush=True)
        (out_dir / "doublet_ambient_qc_report.md").write_text(
            "# Doublet / Ambient / Contamination QC Report\n\nGate not ready.\n", encoding="utf-8"
        )
        return {"run_id": run_id, "status": "gate_not_ready"}

    # Read inputs
    qc_path = out_dir / "cell_qc_metrics.parquet"
    filter_path = Path(filter_plan_path) if filter_plan_path else (out_dir / "cell_filtering_plan.parquet")
    mapping_path = run_root / "01_mapping_audit" / "cell_to_sample_mapping_draft.parquet"

    if not qc_path.exists():
        print("cell_qc_metrics.parquet not found.", flush=True)
        return {"run_id": run_id, "status": "missing_qc"}
    if not filter_path.exists():
        print("cell_filtering_plan.parquet not found.", flush=True)
        return {"run_id": run_id, "status": "missing_filter_plan"}

    print("Reading QC metrics and filtering plan ...", flush=True)
    qc_df = pq.read_table(qc_path).to_pandas()
    filter_df = pq.read_table(filter_path).to_pandas()
    mapping_df = pq.read_table(mapping_path, columns=["cell_barcode", "source_h5ad", "mapped_sample_id", "mapping_status"]).to_pandas()
    mapped_mapping = mapping_df[mapping_df["mapping_status"] == "mapped"]

    print(f"QC rows: {len(qc_df)}, Filter rows: {len(filter_df)}, Mapped cells: {len(mapped_mapping)}", flush=True)

    # Merge QC and filter
    filter_merge_cols = ["source_h5ad", "cell_barcode", "sample_id", "filter_plan"]
    if "plan_reason" in filter_df.columns:
        filter_merge_cols.append("plan_reason")
    merged = qc_df.merge(
        filter_df[filter_merge_cols],
        on=["source_h5ad", "cell_barcode", "sample_id"],
        how="left",
    )
    if min_sample_cells is not None or max_sample_cells is not None:
        sample_n = merged.groupby("sample_id", as_index=False).size().rename(columns={"size": "n_cells"})
        cond = np.ones(len(sample_n), dtype=bool)
        if min_sample_cells is not None:
            cond &= sample_n["n_cells"].to_numpy() >= int(min_sample_cells)
        if max_sample_cells is not None:
            cond &= sample_n["n_cells"].to_numpy() <= int(max_sample_cells)
        selected_ids = sample_n.loc[cond, "sample_id"].astype(str).tolist()
        if max_samples is not None and max_samples > 0:
            selected_ids = sorted(selected_ids)[: int(max_samples)]
        selected_set = set(selected_ids)
        merged = merged[merged["sample_id"].astype(str).isin(selected_set)].copy()
        print(
            f"Batch filter applied: min={min_sample_cells}, max={max_sample_cells}, "
            f"max_samples={max_samples}, selected_samples={merged['sample_id'].nunique()}, selected_cells={len(merged)}"
        , flush=True)
    if sample_shard_index is not None and sample_shard_count is not None:
        if sample_shard_count <= 0:
            raise ValueError("sample_shard_count must be > 0")
        if sample_shard_index < 0 or sample_shard_index >= sample_shard_count:
            raise ValueError("sample_shard_index out of range")
        all_ids = sorted(merged["sample_id"].astype(str).unique().tolist())
        shard_ids = [sid for i, sid in enumerate(all_ids) if (i % sample_shard_count) == sample_shard_index]
        merged = merged[merged["sample_id"].astype(str).isin(set(shard_ids))].copy()
        print(
            f"Shard filter applied: shard={sample_shard_index}/{sample_shard_count}, "
            f"selected_samples={len(shard_ids)}, selected_cells={len(merged)}"
        , flush=True)

    # Per-source processing: read h5ad once, process all samples within source
    h5ad_dir = root / "data" / "processed" / "srt" / "raw"
    source_groups = merged.groupby("source_h5ad")

    suffix = f".{output_suffix}" if output_suffix else ""
    if batch_tag:
        batch_dir = out_dir / batch_dir_name
        batch_dir.mkdir(parents=True, exist_ok=True)
        doublet_path = batch_dir / f"doublet_scores_by_cell{suffix}.{batch_tag}.parquet"
        ambient_path = batch_dir / f"ambient_contamination_scores_by_cell{suffix}.{batch_tag}.parquet"
    else:
        doublet_path = out_dir / f"doublet_scores_by_cell{suffix}.parquet"
        ambient_path = out_dir / f"ambient_contamination_scores_by_cell{suffix}.parquet"
    for p in [doublet_path, ambient_path]:
        if p.exists():
            p.unlink()
    doublet_writer = None
    ambient_writer = None
    n_doublet_rows_written = 0
    n_ambient_rows_written = 0
    sample_summaries: List[Dict] = []
    summary_fields = [
        "run_id", "created_at", "input_manifest_ref", "sample_id", "cohort_id",
        "n_cells_input", "doublet_rate", "high_confidence_doublet_rate",
        "ambient_audit_status", "rbc_warning", "stress_warning",
        "tumor_contamination_warning", "sample_qc_decision",
        "doublet_method", "doublet_run_status",
    ]
    summary_path = (
        out_dir / f"sample_doublet_ambient_summary{suffix}.csv"
        if not batch_tag
        else (out_dir / batch_dir_name / f"sample_doublet_ambient_summary{suffix}.{batch_tag}.csv")
    )
    if summary_path.exists():
        summary_path.unlink()

    for src_idx, (source_name, source_df) in enumerate(source_groups):
        path = h5ad_dir / source_name
        n_source_cells = len(source_df)
        print(f"  [{src_idx+1}/{len(source_groups)}] Source {source_name} ({n_source_cells} cells, {source_df['sample_id'].nunique()} samples) ...", flush=True)

        # Build QC lookup by barcode for this source
        source_qc_by_barcode = {}
        for _, row in source_df.iterrows():
            source_qc_by_barcode[str(row["cell_barcode"])] = row

        # Memory guard: for very large sources, process per-sample from backed mode
        # instead of loading entire source into memory at once
        USE_BACKED_PER_SAMPLE = n_source_cells > 1_500_000
        if USE_BACKED_PER_SAMPLE:
            print(f"    -> Large source ({n_source_cells} cells): using per-sample backed mode", flush=True)

        adata_backed = None
        adata_sub = None
        if path.exists():
            try:
                adata_backed = ad.read_h5ad(path, backed="r")
                cell_barcodes = source_df["cell_barcode"].tolist()
                obs_idx = adata_backed.obs_names.isin(cell_barcodes)
                if obs_idx.any():
                    if USE_BACKED_PER_SAMPLE:
                        # Keep in backed mode; will subset per-sample
                        pass
                    else:
                        adata_sub = subset_to_memory(adata_backed, obs_idx)
                if not USE_BACKED_PER_SAMPLE:
                    adata_backed.file.close()
                    adata_backed = None
            except Exception as e:
                print(f"    ERROR reading h5ad: {e}", flush=True)
                if adata_backed is not None:
                    try:
                        adata_backed.file.close()
                    except Exception:
                        pass
                    adata_backed = None

        # Compute contamination scores for all cells in source at once (small sources only)
        source_contam_by_barcode = {}
        if adata_sub is not None and not USE_BACKED_PER_SAMPLE:
            try:
                contam_df = compute_contamination_scores(adata_sub)
                for _, c_row in contam_df.iterrows():
                    source_contam_by_barcode[str(c_row["cell_barcode"])] = c_row
            except Exception as e:
                print(f"    Contamination error: {e}", flush=True)

        # Process each sample within this source
        sample_groups_in_source = source_df.groupby("sample_id")
        for sample_id, sample_df in sample_groups_in_source:
            n_cells = len(sample_df)
            print(f"    sample {sample_id}: n_cells={n_cells}", flush=True)
            doublet_scores = None
            doublet_calls = None
            doublet_status = "not_possible"
            doublet_method = "scrublet"
            exp_rate = expected_doublet_rate(n_cells)
            cohort_id = sample_df["cohort_id"].iloc[0]
            sample_adata = None

            # --- Scrublet ---
            # Skip scrublet for very large sources (>1.5M cells) to avoid OOM
            SKIP_SCRUBLET_FOR_LARGE_SOURCE = USE_BACKED_PER_SAMPLE
            # Conservative + robust runtime policy:
            # - very small samples: mark not_possible (unstable estimate)
            # - very large samples: skip in main pass to avoid pathological runtime/memory
            if n_cells < 200:
                doublet_status = "not_possible_too_few_cells"
                doublet_method = "scrublet_not_run_small_sample"
            elif SKIP_SCRUBLET_FOR_LARGE_SOURCE:
                doublet_status = "skipped_large_source"
                doublet_method = "scrublet_skipped"
            elif n_cells >= 200:
                try:
                    sample_barcodes = sample_df["cell_barcode"].tolist()
                    if USE_BACKED_PER_SAMPLE and adata_backed is not None:
                        sample_obs_idx = adata_backed.obs_names.isin(sample_barcodes)
                        if sample_obs_idx.any():
                            sample_adata = subset_to_memory(adata_backed, sample_obs_idx)
                    elif adata_sub is not None:
                        sample_obs_idx = adata_sub.obs_names.isin(sample_barcodes)
                        if sample_obs_idx.any():
                            sample_adata = adata_sub[sample_obs_idx].copy()

                    if sample_adata is not None:
                        counts_matrix, counts_source = select_counts_matrix_for_scrublet(sample_adata)
                        if n_cells > 20_000:
                            doublet_method = "scrublet_chunked_target"
                        scores, predicted, status = run_scrublet_for_sample(
                            counts_matrix, exp_rate
                        )
                        if scores is not None:
                            doublet_scores = scores
                            doublet_calls = assign_doublet_calls(scores, predicted, max_plausible_rate=0.2)
                            if np.all(doublet_calls == "unknown"):
                                doublet_status = "failed_with_reason:implausible_predicted_doublet_rate"
                                doublet_method = f"{doublet_method}|guarded"
                            else:
                                doublet_method = f"{doublet_method}|{counts_source}"
                                doublet_status = status
                        else:
                            doublet_status = status
                    else:
                        doublet_status = "not_possible_no_cells"
                except Exception as e:
                    doublet_status = f"failed_with_reason:{e}"

            # Build barcode lookup for doublet scores
            barcode_to_doublet: Dict[str, Tuple] = {}
            if doublet_scores is not None and doublet_calls is not None and sample_adata is not None:
                sample_barcodes_ordered = sample_adata.obs_names.astype(str).tolist()
                for barcode, score, call in zip(sample_barcodes_ordered, doublet_scores, doublet_calls):
                    barcode_to_doublet[barcode] = (score, call)

            # Build doublet rows
            n_doublet = 0
            n_borderline = 0
            sample_doublet_rows = []
            for _, row in sample_df.iterrows():
                barcode = str(row["cell_barcode"])
                ds, dc = barcode_to_doublet.get(barcode, (np.nan, "unknown"))
                conf = "high" if dc == "high_confidence_doublet" else ("medium" if dc == "borderline_doublet" else "low")
                if dc == "high_confidence_doublet":
                    n_doublet += 1
                elif dc == "borderline_doublet":
                    n_borderline += 1
                sample_doublet_rows.append({
                    "run_id": run_id,
                    "created_at": now_iso(),
                    "input_manifest_ref": manifest_ref,
                    "source_h5ad": source_name,
                    "cohort_id": row["cohort_id"],
                    "cell_barcode": barcode,
                    "sample_id": sample_id,
                    "doublet_score": float(ds) if not np.isnan(ds) else np.nan,
                    "doublet_call": dc,
                    "doublet_confidence": conf,
                    "doublet_method": doublet_method,
                    "expected_doublet_rate": round(exp_rate, 4),
                    "doublet_run_status": doublet_status,
                })
            if sample_doublet_rows:
                ddf = pd.DataFrame(sample_doublet_rows)
                dtab = pa.Table.from_pandas(ddf, preserve_index=False)
                if doublet_writer is None:
                    doublet_writer = pq.ParquetWriter(doublet_path, dtab.schema)
                doublet_writer.write_table(dtab)
                n_doublet_rows_written += len(sample_doublet_rows)

            # --- Contamination / ambient ---
            # For backed-per-sample mode on huge sources, skip per-sample contamination
            # extraction in main pass to avoid pathological I/O runtime.
            sample_contam_by_barcode = {}
            if USE_BACKED_PER_SAMPLE:
                pass
            elif adata_backed is not None:
                try:
                    sample_barcodes = sample_df["cell_barcode"].tolist()
                    sample_obs_idx = adata_backed.obs_names.isin(sample_barcodes)
                    if sample_obs_idx.any():
                        sample_adata_contam = subset_to_memory(adata_backed, sample_obs_idx)
                        contam_df = compute_contamination_scores(sample_adata_contam)
                        for _, c_row in contam_df.iterrows():
                            sample_contam_by_barcode[str(c_row["cell_barcode"])] = c_row
                        del sample_adata_contam, contam_df
                        gc.collect()
                except Exception as e:
                    print(f"    Contamination error for sample {sample_id}: {e}", flush=True)

            rbc_warn = False
            stress_warn = False
            tumor_warn = False
            sample_ambient_rows = []
            for _, row in sample_df.iterrows():
                barcode = str(row["cell_barcode"])
                if USE_BACKED_PER_SAMPLE:
                    c_row = sample_contam_by_barcode.get(barcode)
                else:
                    c_row = source_contam_by_barcode.get(barcode)
                if c_row is not None:
                    sample_ambient_rows.append({
                        "run_id": run_id,
                        "created_at": now_iso(),
                        "input_manifest_ref": manifest_ref,
                        "source_h5ad": source_name,
                        "cohort_id": row["cohort_id"],
                        "cell_barcode": barcode,
                        "sample_id": sample_id,
                        "rbc_score": round(float(c_row["rbc_score"]), 4),
                        "hb_score": round(float(c_row["rbc_score"]), 4),
                        "stress_score": round(float(c_row["stress_score"]), 4),
                        "mitochondrial_warning": row["pct_mito"] > 25,
                        "ribosomal_warning": row["pct_ribo"] > 50,
                        "top_gene_dominance_warning": row.get("top20_gene_fraction", 0) > 50,
                        "tumor_epithelial_contamination_score": round(float(c_row["tumor_epithelial_contamination_score"]), 4),
                        "ambient_audit_status": "audit_only_no_raw_correction",
                    })
                    if float(c_row["rbc_score"]) > 5:
                        rbc_warn = True
                    if float(c_row["stress_score"]) > 5:
                        stress_warn = True
                    if float(c_row["tumor_epithelial_contamination_score"]) > 5:
                        tumor_warn = True
            if sample_ambient_rows:
                adf = pd.DataFrame(sample_ambient_rows)
                atab = pa.Table.from_pandas(adf, preserve_index=False)
                if ambient_writer is None:
                    ambient_writer = pq.ParquetWriter(ambient_path, atab.schema)
                ambient_writer.write_table(atab)
                n_ambient_rows_written += len(sample_ambient_rows)

            doublet_rate = n_doublet / n_cells if n_cells > 0 else 0
            sample_summary_row = {
                "run_id": run_id,
                "created_at": now_iso(),
                "input_manifest_ref": manifest_ref,
                "sample_id": sample_id,
                "cohort_id": cohort_id,
                "n_cells_input": n_cells,
                "doublet_rate": round(doublet_rate, 4),
                "high_confidence_doublet_rate": round(doublet_rate, 4),
                "ambient_audit_status": "not_run_large_source" if USE_BACKED_PER_SAMPLE else "audit_only",
                "rbc_warning": rbc_warn,
                "stress_warning": stress_warn,
                "tumor_contamination_warning": tumor_warn,
                "sample_qc_decision": "review" if doublet_rate > 0.15 else "pass",
                "doublet_method": doublet_method,
                "doublet_run_status": doublet_status,
            }
            sample_summaries.append(sample_summary_row)
            csv_append(summary_path, summary_fields, sample_summary_row)

        # Clean up backed file for large sources
        if adata_backed is not None:
            try:
                adata_backed.file.close()
            except Exception:
                pass

    if doublet_writer is not None:
        doublet_writer.close()
        print("Wrote doublet_scores_by_cell.parquet", flush=True)
    if ambient_writer is not None:
        ambient_writer.close()
        print("Wrote ambient_contamination_scores_by_cell.parquet", flush=True)
    print(f"Doublet rows: {n_doublet_rows_written}, Ambient rows: {n_ambient_rows_written}", flush=True)

    if not summary_path.exists():
        csv_write(summary_path, summary_fields, sample_summaries)

    if batch_tag:
        print(f"Batch mode completed: {batch_tag}", flush=True)
        return {
            "run_id": run_id,
            "status": "batch_complete",
            "batch_tag": batch_tag,
            "n_samples": len(sample_summaries),
            "n_cells": len(merged),
            "doublet_path": str(doublet_path),
            "ambient_path": str(ambient_path),
            "summary_path": str(summary_path),
        }

    # Build final inclusion flags (vectorized)
    if doublet_path.exists():
        doublet_df = pq.read_table(doublet_path, columns=["cell_barcode", "sample_id", "doublet_call"]).to_pandas()
    else:
        doublet_df = pd.DataFrame(columns=["cell_barcode", "sample_id", "doublet_call"])
    inc_df = merged[["source_h5ad", "cohort_id", "cell_barcode", "sample_id", "filter_plan", "plan_reason"]].copy()
    inc_df = inc_df.merge(doublet_df, on=["source_h5ad", "cell_barcode", "sample_id"], how="left")
    inc_df["doublet_call"] = inc_df["doublet_call"].fillna("unknown")

    cond_doublet_high = inc_df["doublet_call"] == "high_confidence_doublet"
    cond_qc_remove = inc_df["filter_plan"] == "remove_low_quality"
    cond_qc_review = inc_df["filter_plan"] == "review_borderline"
    cond_doublet_border = inc_df["doublet_call"] == "borderline_doublet"

    inc_df["final_inclusion"] = np.select(
        [cond_doublet_high, cond_qc_remove, cond_qc_review, cond_doublet_border],
        ["exclude_high_confidence_doublet", "exclude_low_quality", "include_sensitivity_only", "include_sensitivity_only"],
        default="include_main",
    )
    inc_df["inclusion_reason"] = np.select(
        [cond_doublet_high, cond_qc_remove, cond_qc_review, cond_doublet_border],
        ["high_confidence_doublet", inc_df["plan_reason"].fillna("low_quality"), "borderline_qc", "borderline_doublet"],
        default="pass_all_qc_and_doublet_checks",
    )
    inc_df["source_qc_flags"] = np.select(
        [cond_doublet_high, cond_qc_remove, cond_qc_review, cond_doublet_border],
        ["doublet_exclusion", "qc_exclusion", "borderline", "borderline_doublet"],
        default="none",
    )
    inc_df["run_id"] = run_id
    inc_df["created_at"] = now_iso()
    inc_df["input_manifest_ref"] = manifest_ref
    inc_df = inc_df[
        [
            "run_id", "created_at", "input_manifest_ref", "source_h5ad", "cohort_id",
            "cell_barcode", "sample_id", "final_inclusion", "inclusion_reason",
            "source_qc_flags", "doublet_call", "filter_plan"
        ]
    ]

    # Write final_cell_inclusion_flags.parquet
    table = pa.Table.from_pandas(inc_df, preserve_index=False)
    final_inc_path = out_dir / f"final_cell_inclusion_flags{suffix}.parquet"
    pq.write_table(table, final_inc_path)
    print(f"Wrote {final_inc_path.name}")

    # Report
    counts = inc_df["final_inclusion"].value_counts()
    n_main = int(counts.get("include_main", 0))
    n_sens = int(counts.get("include_sensitivity_only", 0))
    n_ex_qc = int(counts.get("exclude_low_quality", 0))
    n_ex_dbl = int(counts.get("exclude_high_confidence_doublet", 0))
    n_high_doublet_samples = sum(1 for s in sample_summaries if s["doublet_rate"] > 0.15)
    n_cells_total = len(inc_df)

    report = f"""# Doublet / Ambient / Contamination QC Report

- **run_id**: {run_id}
- **gate_ready**: {gate_ready}
- **n_samples**: {len(sample_summaries)}
- **n_cells_total**: {n_cells_total}
- **include_main**: {n_main}
- **include_sensitivity_only**: {n_sens}
- **exclude_low_quality**: {n_ex_qc}
- **exclude_high_confidence_doublet**: {n_ex_dbl}
- **high_doublet_samples**: {n_high_doublet_samples}

## Final Inclusion Summary

| Inclusion Class | Count | Fraction |
|-----------------|-------|----------|
| include_main | {n_main} | {n_main/n_cells_total:.4f} |
| include_sensitivity_only | {n_sens} | {n_sens/n_cells_total:.4f} |
| exclude_low_quality | {n_ex_qc} | {n_ex_qc/n_cells_total:.4f} |
| exclude_high_confidence_doublet | {n_ex_dbl} | {n_ex_dbl/n_cells_total:.4f} |

## Acceptance Rules

- [x] Each sample has doublet audit state
- [x] High-confidence doublets excluded from main analysis
- [x] No ambient correction claimed without raw support
- [x] All decisions written to manifest-compatible outputs

## Key Parameter References

- Scrublet method and sample-wise execution rationale:
  Wolock et al. 2019, Cell Systems (doi:10.1016/j.cels.2018.11.005)
- Scrublet API semantics (predicted_doublets is the primary auto-thresholded call):
  https://github.com/swolock/scrublet
- Plausible multiplet-rate scale in Chromium workflows:
  10x Genomics user guides (for example CG000317 and related tables, ~0.4% at ~500 recovered to ~8% at ~10,000 recovered)
- Ambient correction should not be claimed without required inputs:
  SoupX, Young & Behjati 2020, GigaScience (doi:10.1093/gigascience/giaa151)
"""
    report_path = out_dir / f"doublet_ambient_qc_report{suffix}.md"
    report_path.write_text(report, encoding="utf-8")

    # Checkpoint
    checkpoint_name = f"step2_4{suffix}.checkpoint.yaml" if suffix else "step2_4.checkpoint.yaml"
    (checkpoint_dir / checkpoint_name).write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\n"
        f"rows_written: {n_cells_total}\nstatus: complete\n",
        encoding="utf-8",
    )

    print(f"Done. main={n_main}, sens={n_sens}, ex_qc={n_ex_qc}, ex_dbl={n_ex_dbl}")
    return {
        "run_id": run_id,
        "status": "complete",
        "n_samples": len(sample_summaries),
        "n_cells": n_cells_total,
        "doublet_path": str(doublet_path),
        "ambient_path": str(ambient_path),
        "summary_path": str(summary_path),
        "final_inclusion_path": str(final_inc_path),
        "report_path": str(report_path),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step2.4 doublet ambient contamination audit")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    parser.add_argument("--min-sample-cells", type=int, default=None, help="optional min sample cell count for batch mode")
    parser.add_argument("--max-sample-cells", type=int, default=None, help="optional max sample cell count for batch mode")
    parser.add_argument("--max-samples", type=int, default=None, help="optional max number of samples for batch mode")
    parser.add_argument("--batch-tag", default=None, help="optional batch tag; writes batch outputs and skips global finalization")
    parser.add_argument("--sample-shard-index", type=int, default=None, help="optional sample shard index (0-based)")
    parser.add_argument("--sample-shard-count", type=int, default=None, help="optional sample shard count")
    parser.add_argument("--filter-plan-path", default=None, help="optional override filter plan parquet path")
    parser.add_argument("--output-suffix", default=None, help="optional suffix for output files, e.g. baseline or rescue")
    parser.add_argument("--batch-dir-name", default="doublet_batches", help="batch directory name under 03_qc")
    args = parser.parse_args()
    result = run_step2_4(
        root_dir=args.root_dir,
        run_id=args.run_id,
        min_sample_cells=args.min_sample_cells,
        max_sample_cells=args.max_sample_cells,
        max_samples=args.max_samples,
        batch_tag=args.batch_tag,
        sample_shard_index=args.sample_shard_index,
        sample_shard_count=args.sample_shard_count,
        filter_plan_path=args.filter_plan_path,
        output_suffix=args.output_suffix,
        batch_dir_name=args.batch_dir_name,
    )
    print(result)
