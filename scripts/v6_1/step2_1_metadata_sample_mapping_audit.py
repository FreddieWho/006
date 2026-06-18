from __future__ import annotations

import argparse
import csv
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]

CANDIDATE_SAMPLE_FIELDS = [
    "sample_id",
    "sample_prefix",
    "sample",
    "Sample",
    "patient_id",
    "patient",
    "Patient",
    "donor_id",
    "donor",
    "Donor",
    "orig.ident",
    "library_id",
    "batch",
    "Batch",
]

REQUIRED_OUT_COLS = [
    "run_id",
    "created_at",
    "input_manifest_ref",
    "cell_barcode",
    "source_h5ad",
    "source_obs_table",
    "obs_sample_raw",
    "mapped_sample_id",
    "mapped_patient_id",
    "mapped_cohort_id",
    "mapped_timepoint",
    "mapped_tissue_source",
    "split",
    "mapping_status",
    "mapping_confidence",
    "mapping_reason",
    "alias_used",
    "candidate_sample_ids",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def csv_write(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
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
            if row.get("step_name") == "Step2.1_mapping":
                status = row.get("status")
                fatal = row.get("fatal_inputs_missing", "")
                if status == "ready":
                    return True
                if status == "blocked" and fatal in {"", "mapping_acceptance"}:
                    return True
                return False
    return False


def infer_cohort_from_filename(filename: str) -> str:
    stem = Path(filename).stem.lower()
    # Known mappings for filenames that don't match metadata cohort_id
    mapping = {
        "gse140228_smartseq2": "GSE140228_smartseq2",
        "gse140228_10x": "GSE140228_droplet",
        "gse140228_droplet": "GSE140228_droplet",
        "imbrave150": "IMbrave150",
        "lambrecht_hcc": "lambrecht_hcc",
        "task01": "TASK01",
        "task02": "TASK02",
    }
    if stem in mapping:
        return mapping[stem]
    # GSE prefix: preserve full stem with underscores (e.g., gse123813_bcc)
    if stem.startswith("gse"):
        return stem.upper()
    return stem.upper()


def _scalar_str(series: pd.Series, idx: str) -> str:
    """Safely extract scalar string from Series, handling duplicate indices."""
    val = series[idx]
    if isinstance(val, pd.Series):
        val = val.iloc[0]
    return str(val)


def normalize_id(val: str) -> str:
    return str(val).strip().lower()


def resolve_canonical_cohort_id(raw_cohort: str, sample_meta: pd.DataFrame) -> str:
    if sample_meta.empty or "cohort_id" not in sample_meta.columns:
        return raw_cohort
    cohort_ids = [str(x) for x in sample_meta["cohort_id"].dropna().unique()]
    exact = {cid: cid for cid in cohort_ids}
    lower = {cid.lower(): cid for cid in cohort_ids}
    norm = {normalize_id(cid): cid for cid in cohort_ids}
    if raw_cohort in exact:
        return exact[raw_cohort]
    if raw_cohort.lower() in lower:
        return lower[raw_cohort.lower()]
    return norm.get(normalize_id(raw_cohort), raw_cohort)


def build_sample_lookup(sample_meta: pd.DataFrame) -> Tuple[Dict[Tuple[str, str], pd.Series], Dict[Tuple[str, str], List[pd.Series]]]:
    """Build (cohort, sample_id) -> row and (cohort, patient_id) -> [rows] lookups."""
    sample_lookup: Dict[Tuple[str, str], pd.Series] = {}
    patient_lookup: Dict[Tuple[str, str], List[pd.Series]] = {}
    for _, row in sample_meta.iterrows():
        cohort = str(row.get("cohort_id", ""))
        sample_id = str(row.get("sample_id", ""))
        patient_id = str(row.get("patient_id", ""))
        if cohort and sample_id:
            sample_lookup[(cohort, sample_id)] = row
            sample_lookup[(cohort, normalize_id(sample_id))] = row
        if cohort and patient_id:
            patient_lookup.setdefault((cohort, patient_id), []).append(row)
            patient_lookup.setdefault((cohort, normalize_id(patient_id)), []).append(row)
    return sample_lookup, patient_lookup


def build_gse272734_week_lookup(sample_meta: pd.DataFrame) -> Dict[Tuple[str, str], pd.Series]:
    lookup: Dict[Tuple[str, str], pd.Series] = {}
    if sample_meta.empty:
        return lookup
    for _, row in sample_meta.iterrows():
        cohort = str(row.get("cohort_id", ""))
        if cohort != "GSE272734":
            continue
        patient_id = str(row.get("patient_id", ""))
        sample_id = str(row.get("sample_id", ""))
        match = re.search(r"-wk(\d+)$", sample_id)
        if not patient_id or not match:
            continue
        week_token = f"wk{match.group(1)}"
        lookup[(patient_id, week_token)] = row
        lookup[(normalize_id(patient_id), week_token)] = row
    return lookup


def gse272993_timepoint_to_week(raw_timepoint: str) -> str:
    token = normalize_id(raw_timepoint)
    mapping = {
        "baseline": "wk0",
        "follow up 1": "wk3",
        "follow up 2": "wk6",
        "follow up 3": "wk9",
        "follow up 4": "wk12",
    }
    return mapping.get(token, "")


def discover_sources(root: Path, manifest: Dict[str, str]) -> List[Dict[str, str]]:
    sources = []
    snapshot_ref = manifest.get("data_pool_snapshot_ref", "")
    inventory_path = None
    if snapshot_ref:
        candidate = Path(snapshot_ref).parent / "h5ad_inventory.csv"
        if candidate.exists():
            inventory_path = candidate

    if inventory_path is not None:
        with inventory_path.open("r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                path = row.get("file_path", "")
                cohort = row.get("cohort_id", "")
                if path:
                    sources.append({"path": path, "type": "h5ad", "name": Path(path).name, "cohort": cohort})
        return sources

    h5ad_dir = root / "data" / "processed" / "srt" / "raw"
    if h5ad_dir.exists():
        for f in sorted(h5ad_dir.glob("*.h5ad")):
            cohort = infer_cohort_from_filename(f.name)
            sources.append({"path": str(f), "type": "h5ad", "name": f.name, "cohort": cohort})
    return sources


def extract_obs_fields(path: str) -> List[str]:
    try:
        import h5py
        with h5py.File(path, "r") as f:
            obs = f["obs"]
            cols = [k for k in obs.keys() if not k.startswith("_")]
        return cols
    except Exception:
        return []


def _build_cohort_lookups(
    cohort: str,
    sample_lookup: Dict[Tuple[str, str], pd.Series],
    patient_lookup: Dict[Tuple[str, str], List[pd.Series]],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]], Dict[str, List[pd.Series]], Dict[str, List[pd.Series]]]:
    """Build cohort-specific exact and normalized lookups.
    Sample lookups are per-column dicts for fast vectorized mapping.
    Patient lookups remain row lists for ambiguity handling."""
    cols = ["sample_id", "patient_id", "timepoint", "tissue_source", "split", "split_label"]
    sample_exact: Dict[str, Dict[str, str]] = {c: {} for c in cols}
    sample_norm: Dict[str, Dict[str, str]] = {c: {} for c in cols}
    patient_exact: Dict[str, List[pd.Series]] = {}
    patient_norm: Dict[str, List[pd.Series]] = {}

    for (c, k), v in sample_lookup.items():
        if c == cohort:
            is_norm = k == k.lower()
            maps = sample_norm if is_norm else sample_exact
            for col in cols:
                maps[col][k] = str(v.get(col, ""))

    for (c, k), v in patient_lookup.items():
        if c == cohort:
            if k == k.lower():
                patient_norm[k] = v
            else:
                patient_exact[k] = v

    return sample_exact, sample_norm, patient_exact, patient_norm


def map_cells_from_source(
    source: Dict[str, str],
    sample_lookup: Dict[Tuple[str, str], pd.Series],
    patient_lookup: Dict[Tuple[str, str], List[pd.Series]],
    gse272734_week_lookup: Dict[Tuple[str, str], pd.Series],
    run_id: str,
    manifest_ref: str,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    """Map cells from one source using vectorized pandas operations."""
    path = source["path"]
    cohort = source["cohort"]
    source_name = source["name"]

    try:
        a = ad.read_h5ad(path, backed="r")
        obs = a.obs
        n_cells = a.n_obs
        a.file.close()
    except Exception as e:
        return [], {"source": source_name, "error": str(e), "n_cells": 0, "n_mapped": 0}

    present_candidates = [c for c in CANDIDATE_SAMPLE_FIELDS if c in obs.columns]
    best_field = present_candidates[0] if present_candidates else None

    if "cohort_id" in obs.columns:
        obs_cohorts = obs["cohort_id"].dropna().unique()
        if len(obs_cohorts) == 1:
            cohort = str(obs_cohorts[0])

    if best_field is None or n_cells == 0:
        return [], {
            "source": source_name, "cohort": cohort, "n_cells": n_cells,
            "n_mapped": 0, "n_ambiguous": 0, "n_unmapped": n_cells,
            "candidate_fields": ";".join(present_candidates), "best_field_used": "none",
        }

    # Build cohort-specific lookups
    sample_exact, sample_norm, patient_exact, patient_norm = _build_cohort_lookups(cohort, sample_lookup, patient_lookup)

    rescue_possible = cohort == "GSE272993" and bool(gse272734_week_lookup) and "patient_id" in obs.columns and ("timepoint_raw" in obs.columns or "timepoint" in obs.columns)
    if not rescue_possible and not sample_exact["sample_id"] and not patient_exact and not sample_norm["sample_id"] and not patient_norm:
        return [], {
            "source": source_name,
            "cohort": cohort,
            "n_cells": n_cells,
            "counted_cells_for_acceptance": 0,
            "n_mapped": 0,
            "n_ambiguous": 0,
            "n_unmapped": n_cells,
            "candidate_fields": ";".join(present_candidates),
            "best_field_used": best_field or "none",
            "source_status": "excluded_unregistered_source",
        }

    raw_vals = obs[best_field].astype(str)
    norm_vals = raw_vals.str.strip().str.lower()

    # Initialize result series
    idx = obs.index
    mapped_sample_id = pd.Series("", index=idx)
    mapped_patient_id = pd.Series("", index=idx)
    mapped_timepoint = pd.Series("", index=idx)
    mapped_tissue = pd.Series("", index=idx)
    mapped_cohort_id = pd.Series(cohort, index=idx)
    split = pd.Series("", index=idx)
    status = pd.Series("unmapped", index=idx)
    confidence = pd.Series("low", index=idx)
    reason = pd.Series("no_metadata_match", index=idx)
    candidates = pd.Series("", index=idx)

    # Phase 1: exact sample match
    exact_sample_id = raw_vals.map(sample_exact["sample_id"])
    exact_mask = exact_sample_id.notna()
    if exact_mask.any():
        mapped_sample_id[exact_mask] = exact_sample_id[exact_mask]
        mapped_patient_id[exact_mask] = raw_vals[exact_mask].map(sample_exact["patient_id"])
        mapped_timepoint[exact_mask] = raw_vals[exact_mask].map(sample_exact["timepoint"]).fillna("unknown")
        mapped_tissue[exact_mask] = raw_vals[exact_mask].map(sample_exact["tissue_source"]).fillna("tissue_unknown")
        mapped_cohort_id[exact_mask] = cohort
        s_split = raw_vals[exact_mask].map(sample_exact["split"]).fillna("")
        s_label = raw_vals[exact_mask].map(sample_exact["split_label"]).fillna("")
        split[exact_mask] = np.where(s_split == "", s_label, s_split)
        status[exact_mask] = "mapped"
        confidence[exact_mask] = "high"
        reason[exact_mask] = "exact_sample_match"

    # Phase 2: normalized sample match for unmatched
    unmatched = ~exact_mask
    if unmatched.any():
        norm_sample_id = norm_vals[unmatched].map(sample_norm["sample_id"])
        norm_mask = norm_sample_id.notna()
        if norm_mask.any():
            idx_norm = norm_sample_id[norm_mask].index
            mapped_sample_id[idx_norm] = norm_sample_id[norm_mask]
            mapped_patient_id[idx_norm] = norm_vals[idx_norm].map(sample_norm["patient_id"])
            mapped_timepoint[idx_norm] = norm_vals[idx_norm].map(sample_norm["timepoint"]).fillna("unknown")
            mapped_tissue[idx_norm] = norm_vals[idx_norm].map(sample_norm["tissue_source"]).fillna("tissue_unknown")
            mapped_cohort_id[idx_norm] = cohort
            s_split = norm_vals[idx_norm].map(sample_norm["split"]).fillna("")
            s_label = norm_vals[idx_norm].map(sample_norm["split_label"]).fillna("")
            split[idx_norm] = np.where(s_split == "", s_label, s_split)
            status[idx_norm] = "mapped"
            confidence[idx_norm] = "high"
            reason[idx_norm] = "exact_sample_match"
            full_norm_mask = pd.Series(False, index=idx)
            full_norm_mask[idx_norm] = True
            unmatched = unmatched & ~full_norm_mask

    # Phase 2.5: explicit GSE272993 -> GSE272734 rescue for cross-cohort cells with missing sample field
    if unmatched.any() and cohort == "GSE272993" and gse272734_week_lookup and "patient_id" in obs.columns:
        patient_series = obs["patient_id"].astype(str)
        timepoint_col = "timepoint_raw" if "timepoint_raw" in obs.columns else "timepoint" if "timepoint" in obs.columns else None
        if timepoint_col is not None:
            timepoint_series = obs[timepoint_col].astype(str)
            rescue_mask = unmatched & norm_vals.isin({"", "nan", "na", "none", "null", "n/a"})
            if rescue_mask.any():
                rescue_idx = idx[rescue_mask]
                for cell_id in rescue_idx:
                    patient_id = patient_series[cell_id]
                    week_token = gse272993_timepoint_to_week(timepoint_series[cell_id])
                    meta_row = gse272734_week_lookup.get((patient_id, week_token))
                    if meta_row is None:
                        meta_row = gse272734_week_lookup.get((normalize_id(patient_id), week_token))
                    if meta_row is None:
                        continue
                    mapped_sample_id[cell_id] = str(meta_row.get("sample_id", ""))
                    mapped_patient_id[cell_id] = str(meta_row.get("patient_id", ""))
                    mapped_timepoint[cell_id] = str(meta_row.get("timepoint", "unknown"))
                    mapped_tissue[cell_id] = str(meta_row.get("tissue_source", "tissue_unknown"))
                    mapped_cohort_id[cell_id] = str(meta_row.get("cohort_id", "GSE272734"))
                    split[cell_id] = str(meta_row.get("split", "")) or str(meta_row.get("split_label", ""))
                    status[cell_id] = "mapped"
                    confidence[cell_id] = "medium"
                    reason[cell_id] = "gse272993_cross_cohort_patient_timepoint_rescue"
                unmatched = status != "mapped"

    # Phase 3: exact patient match for remaining unmatched
    if unmatched.any():
        patient_matched = raw_vals[unmatched].map(patient_exact)
        patient_mask = patient_matched.notna()
        if patient_mask.any():
            patient_vals = raw_vals[unmatched][patient_mask]
            patient_full_mask = pd.Series(False, index=idx)
            for pval in patient_vals.unique():
                full_mask = (raw_vals == pval) & unmatched
                patient_full_mask = patient_full_mask | full_mask
                norm_pval = str(pval).strip().lower()
                rows_list = patient_exact.get(pval, [])
                if not rows_list:
                    rows_list = patient_norm.get(norm_pval, [])
                if len(rows_list) == 1:
                    meta_row = rows_list[0]
                    mapped_sample_id[full_mask] = str(meta_row.get("sample_id", ""))
                    mapped_patient_id[full_mask] = str(meta_row.get("patient_id", ""))
                    mapped_timepoint[full_mask] = str(meta_row.get("timepoint", "unknown"))
                    mapped_tissue[full_mask] = str(meta_row.get("tissue_source", "tissue_unknown"))
                    mapped_cohort_id[full_mask] = str(meta_row.get("cohort_id", cohort))
                    split[full_mask] = str(meta_row.get("split", "")) or str(meta_row.get("split_label", ""))
                    status[full_mask] = "mapped"
                    confidence[full_mask] = "medium"
                    reason[full_mask] = "patient_match_single_sample"
                elif len(rows_list) > 1:
                    cands = ";".join(str(r.get("sample_id", "")) for r in rows_list)
                    mapped_patient_id[full_mask] = str(rows_list[0].get("patient_id", ""))
                    candidates[full_mask] = cands
                    status[full_mask] = "ambiguous"
                    confidence[full_mask] = "low"
                    reason[full_mask] = "patient_match_multiple_samples"
            unmatched = unmatched & ~patient_full_mask

    # Phase 4: normalized patient match for remaining unmatched
    if unmatched.any():
        norm_patient_matched = norm_vals[unmatched].map(patient_norm)
        norm_patient_mask = norm_patient_matched.notna()
        if norm_patient_mask.any():
            norm_patient_vals = norm_vals[unmatched][norm_patient_mask]
            norm_patient_full_mask = pd.Series(False, index=idx)
            for pval in norm_patient_vals.unique():
                full_mask = (norm_vals == pval) & unmatched
                norm_patient_full_mask = norm_patient_full_mask | full_mask
                rows_list = patient_norm.get(pval, [])
                if not rows_list:
                    rows_list = patient_exact.get(pval, [])
                if len(rows_list) == 1:
                    meta_row = rows_list[0]
                    mapped_sample_id[full_mask] = str(meta_row.get("sample_id", ""))
                    mapped_patient_id[full_mask] = str(meta_row.get("patient_id", ""))
                    mapped_timepoint[full_mask] = str(meta_row.get("timepoint", "unknown"))
                    mapped_tissue[full_mask] = str(meta_row.get("tissue_source", "tissue_unknown"))
                    mapped_cohort_id[full_mask] = str(meta_row.get("cohort_id", cohort))
                    split[full_mask] = str(meta_row.get("split", "")) or str(meta_row.get("split_label", ""))
                    status[full_mask] = "mapped"
                    confidence[full_mask] = "medium"
                    reason[full_mask] = "patient_match_single_sample"
                elif len(rows_list) > 1:
                    cands = ";".join(str(r.get("sample_id", "")) for r in rows_list)
                    mapped_patient_id[full_mask] = str(rows_list[0].get("patient_id", ""))
                    candidates[full_mask] = cands
                    status[full_mask] = "ambiguous"
                    confidence[full_mask] = "low"
                    reason[full_mask] = "patient_match_multiple_samples"
            unmatched = unmatched & ~norm_patient_full_mask

    # Build output dataframe
    out_df = pd.DataFrame({
        "run_id": run_id,
        "created_at": now_iso(),
        "input_manifest_ref": manifest_ref,
        "cell_barcode": obs.index.astype(str),
        "source_h5ad": source_name,
        "source_obs_table": source_name,
        "obs_sample_raw": raw_vals,
        "mapped_sample_id": mapped_sample_id,
        "mapped_patient_id": mapped_patient_id,
        "mapped_cohort_id": mapped_cohort_id,
        "mapped_timepoint": mapped_timepoint,
        "mapped_tissue_source": mapped_tissue.replace("", "tissue_unknown"),
        "split": split,
        "mapping_status": status,
        "mapping_confidence": confidence,
        "mapping_reason": reason,
        "alias_used": best_field,
        "candidate_sample_ids": candidates,
    })

    for col in REQUIRED_OUT_COLS:
        if col not in out_df.columns:
            out_df[col] = ""
    out_df = out_df[REQUIRED_OUT_COLS]

    rows = out_df.to_dict("records")
    mapped = (status == "mapped").sum()
    ambiguous = (status == "ambiguous").sum()

    audit = {
        "source": source_name,
        "cohort": cohort,
        "n_cells": n_cells,
        "counted_cells_for_acceptance": n_cells,
        "n_mapped": int(mapped),
        "n_ambiguous": int(ambiguous),
        "n_unmapped": n_cells - int(mapped) - int(ambiguous),
        "candidate_fields": ";".join(present_candidates),
        "best_field_used": best_field or "none",
        "source_status": "counted",
    }
    return rows, audit


def rows_to_table(rows: List[Dict[str, object]]) -> pa.Table:
    df = pd.DataFrame(rows)
    for col in REQUIRED_OUT_COLS:
        if col not in df.columns:
            df[col] = ""
    df = df[REQUIRED_OUT_COLS]
    return pa.Table.from_pandas(df, preserve_index=False)


def run_step2_1_mapping_audit(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    out_dir = run_root / "01_mapping_audit"
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_patched_manifest(run_root)
    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    # Gate check
    gate_ready = check_gate(run_root)

    # Load metadata
    sample_meta_path = manifest.get("sample_metadata_ref", "")
    patient_meta_path = manifest.get("patient_metadata_ref", "")
    split_path = manifest.get("patient_split_ref", "")

    sample_meta = pd.read_csv(sample_meta_path) if sample_meta_path and Path(sample_meta_path).exists() else pd.DataFrame()
    patient_meta = pd.read_csv(patient_meta_path) if patient_meta_path and Path(patient_meta_path).exists() else pd.DataFrame()
    split_df = pd.read_csv(split_path) if split_path and Path(split_path).exists() else pd.DataFrame()

    sample_lookup: Dict[Tuple[str, str], pd.Series] = {}
    patient_lookup: Dict[Tuple[str, str], List[pd.Series]] = {}
    if not sample_meta.empty:
        if "cohort_id" in sample_meta.columns:
            sample_meta = sample_meta.copy()
            sample_meta["cohort_id"] = sample_meta["cohort_id"].astype(str)
        sample_lookup, patient_lookup = build_sample_lookup(sample_meta)
    gse272734_week_lookup = build_gse272734_week_lookup(sample_meta)

    # Discover sources
    sources = discover_sources(root, manifest)
    for src in sources:
        src["cohort"] = resolve_canonical_cohort_id(src["cohort"], sample_meta)

    # Phase 1: field audit
    field_audits: List[Dict[str, object]] = []
    for src in sources:
        cols = extract_obs_fields(src["path"])
        candidates = [c for c in CANDIDATE_SAMPLE_FIELDS if c in cols]
        field_audits.append({
            "source_name": src["name"],
            "cohort": src["cohort"],
            "obs_columns": ";".join(cols),
            "candidate_sample_fields": ";".join(candidates),
            "has_sample_id": "sample_id" in cols,
            "has_patient_id": "patient_id" in cols,
            "has_patient": "patient" in cols,
        })

    # Phase 2: cell mapping (only if gate ready)
    mapping_audits: List[Dict[str, object]] = []
    parquet_path = out_dir / "cell_to_sample_mapping_draft.parquet"
    if parquet_path.exists():
        parquet_path.unlink()

    writer = None
    if gate_ready and not sample_meta.empty:
        for idx_src, src in enumerate(sources):
            print(f"[{idx_src+1}/{len(sources)}] Mapping {src['name']} ...", flush=True)
            rows, audit = map_cells_from_source(src, sample_lookup, patient_lookup, gse272734_week_lookup, run_id, manifest_ref)
            if rows:
                table = rows_to_table(rows)
                if writer is None:
                    writer = pq.ParquetWriter(parquet_path, table.schema)
                writer.write_table(table)
            mapping_audits.append(audit)
            print(f"  -> {audit.get('n_cells', 0)} cells, {audit.get('n_mapped', 0)} mapped", flush=True)
        if writer is not None:
            writer.close()
    else:
        mapping_audits.append({"note": "gate_not_ready", "reason": "Step2.1_mapping gate not ready or metadata missing"})

    # Write sample_mapping_audit.csv
    audit_rows = []
    for fa in field_audits:
        ma = next((m for m in mapping_audits if m.get("source") == fa["source_name"]), {})
        audit_rows.append({
            "source_name": fa["source_name"],
            "cohort": fa["cohort"],
            "obs_columns_count": len(fa["obs_columns"].split(";")) if fa["obs_columns"] else 0,
            "candidate_sample_fields": fa["candidate_sample_fields"],
            "has_sample_id": fa["has_sample_id"],
            "has_patient_id": fa["has_patient_id"],
            "has_patient": fa["has_patient"],
            "n_cells": ma.get("n_cells", 0),
            "counted_cells_for_acceptance": ma.get("counted_cells_for_acceptance", ma.get("n_cells", 0)),
            "n_mapped": ma.get("n_mapped", 0),
            "n_ambiguous": ma.get("n_ambiguous", 0),
            "n_unmapped": ma.get("n_unmapped", 0),
            "best_field_used": ma.get("best_field_used", "n/a"),
            "source_status": ma.get("source_status", "counted"),
        })

    audit_path = out_dir / "sample_mapping_audit.csv"
    csv_write(audit_path, [
        "source_name", "cohort", "obs_columns_count", "candidate_sample_fields",
        "has_sample_id", "has_patient_id", "has_patient", "n_cells", "counted_cells_for_acceptance", "n_mapped",
        "n_ambiguous", "n_unmapped", "best_field_used", "source_status",
    ], audit_rows)

    # Compute totals
    total_cells = sum(a.get("counted_cells_for_acceptance", a.get("n_cells", 0)) for a in mapping_audits if "n_cells" in a)
    total_mapped = sum(a.get("n_mapped", 0) for a in mapping_audits if "n_mapped" in a)
    total_ambiguous = sum(a.get("n_ambiguous", 0) for a in mapping_audits if "n_ambiguous" in a)
    mapped_fraction = total_mapped / total_cells if total_cells > 0 else 0.0

    # Unmapped cells report
    unmapped_rows = []
    if gate_ready and parquet_path.exists():
        df = pq.read_table(parquet_path).to_pandas()
        unmapped_df = df[df["mapping_status"] != "mapped"]
        if len(unmapped_df) > 0:
            counted_sources = {row["source_name"] for row in audit_rows if row["source_status"] == "counted"}
            unmapped_df = unmapped_df[unmapped_df["source_h5ad"].isin(counted_sources)]
            summary = unmapped_df.groupby("source_h5ad").size().reset_index(name="n_unmapped")
            for _, row in summary.iterrows():
                unmapped_rows.append({
                    "source_h5ad": row["source_h5ad"],
                    "n_unmapped": row["n_unmapped"],
                    "unmapped_fraction": round(row["n_unmapped"] / total_cells, 6) if total_cells > 0 else 0,
                })

    unmapped_path = out_dir / "unmapped_cells_report.csv"
    csv_write(unmapped_path, ["source_h5ad", "n_unmapped", "unmapped_fraction"], unmapped_rows)

    # Unmapped samples report
    # Samples in metadata but not found in any h5ad mapping
    mapped_samples = set()
    if gate_ready and parquet_path.exists():
        df = pq.read_table(parquet_path).to_pandas()
        mapped_samples = set(df[df["mapped_sample_id"] != ""]["mapped_sample_id"].unique())

    all_meta_samples = set(sample_meta["sample_id"].astype(str)) if not sample_meta.empty and "sample_id" in sample_meta.columns else set()
    unmapped_samples = all_meta_samples - mapped_samples
    unmapped_sample_rows = []
    for s in unmapped_samples:
        rows_meta = sample_meta[sample_meta["sample_id"].astype(str) == s]
        if not rows_meta.empty:
            r = rows_meta.iloc[0]
            unmapped_sample_rows.append({
                "sample_id": s,
                "cohort_id": str(r.get("cohort_id", "")),
                "patient_id": str(r.get("patient_id", "")),
                "reason": "no_cell_mapped_from_any_h5ad",
            })

    unmapped_sample_path = out_dir / "unmapped_samples_report.csv"
    csv_write(unmapped_sample_path, ["sample_id", "cohort_id", "patient_id", "reason"], unmapped_sample_rows)

    # Conflict report
    conflict_rows = []
    if gate_ready and parquet_path.exists():
        df = pq.read_table(parquet_path).to_pandas()
        ambig_df = df[df["mapping_status"] == "ambiguous"]
        if len(ambig_df) > 0:
            summary = ambig_df.groupby("source_h5ad").size().reset_index(name="n_ambiguous")
            for _, row in summary.iterrows():
                conflict_rows.append({
                    "source_h5ad": row["source_h5ad"],
                    "n_ambiguous_cells": row["n_ambiguous"],
                    "conflict_type": "multiple_samples_per_patient",
                })

    conflict_path = out_dir / "mapping_conflict_report.csv"
    csv_write(conflict_path, ["source_h5ad", "n_ambiguous_cells", "conflict_type"], conflict_rows)

    # Patient timepoint consistency
    consistency_rows = []
    if gate_ready and parquet_path.exists():
        df = pq.read_table(parquet_path).to_pandas()
        mapped_df = df[df["mapping_status"] == "mapped"]
        if len(mapped_df) > 0:
            # Group by patient and check timepoint consistency
            patient_timepoints = mapped_df.groupby("mapped_patient_id")["mapped_timepoint"].unique().reset_index()
            for _, row in patient_timepoints.iterrows():
                tps = [t for t in row["mapped_timepoint"] if t and t != "unknown"]
                if len(set(tps)) > 1:
                    consistency_rows.append({
                        "patient_id": row["mapped_patient_id"],
                        "timepoints_found": ";".join(set(tps)),
                        "status": "multiple_timepoints",
                    })

    consistency_path = out_dir / "patient_timepoint_consistency.csv"
    csv_write(consistency_path, ["patient_id", "timepoints_found", "status"], consistency_rows)

    # Blocking issues
    blocking_issues: List[str] = []
    if mapped_fraction < 0.99:
        blocking_issues.append(f"mapped_fraction={mapped_fraction:.4f} < 0.99")
    if total_ambiguous > 0:
        blocking_issues.append(f"ambiguous_cells={total_ambiguous}")

    # Check sample->patient is 1:1 in metadata
    if not sample_meta.empty and "sample_id" in sample_meta.columns and "patient_id" in sample_meta.columns:
        sp_pairs = sample_meta.groupby("sample_id")["patient_id"].nunique()
        multi_patient_samples = sp_pairs[sp_pairs > 1]
        if len(multi_patient_samples) > 0:
            blocking_issues.append(f"sample_to_patient_not_1:1 count={len(multi_patient_samples)}")

    blocking_path = out_dir / "mapping_blocking_issues.md"
    blocking_lines = [
        "# Mapping Blocking Issues",
        "",
        f"- **run_id**: {run_id}",
        f"- **gate_ready**: {gate_ready}",
        f"- **total_cells**: {total_cells}",
        f"- **total_mapped**: {total_mapped}",
        f"- **mapped_fraction**: {mapped_fraction:.6f}",
        f"- **total_ambiguous**: {total_ambiguous}",
        f"- **blocking_issues**: {len(blocking_issues)}",
    ]
    if blocking_issues:
        blocking_lines.extend(["", "## Issues"])
        for issue in blocking_issues:
            blocking_lines.append(f"- {issue}")
    else:
        blocking_lines.append("- No blocking issues")
    blocking_path.write_text("\n".join(blocking_lines) + "\n", encoding="utf-8")

    # Mapping summary
    summary_lines = [
        "# Mapping Summary",
        "",
        f"- **run_id**: {run_id}",
        f"- **gate_ready**: {gate_ready}",
        f"- **n_sources**: {len(sources)}",
        f"- **total_cells**: {total_cells}",
        f"- **total_mapped**: {total_mapped}",
        f"- **mapped_fraction**: {mapped_fraction:.6f}",
        f"- **total_ambiguous**: {total_ambiguous}",
        f"- **unmapped_samples_in_meta**: {len(unmapped_sample_rows)}",
        f"- **patient_timepoint_conflicts**: {len(consistency_rows)}",
        "",
        "## Source Details",
        "",
        "| Source | Cohort | Cells | Counted | Mapped | Fraction | Best Field | Status |",
        "|--------|--------|-------|---------|--------|----------|------------|--------|",
    ]
    for ar in audit_rows:
        denom = ar["counted_cells_for_acceptance"] if ar["counted_cells_for_acceptance"] > 0 else ar["n_cells"]
        frac = round(ar["n_mapped"] / denom, 4) if denom > 0 else 0
        summary_lines.append(
            f"| {ar['source_name']} | {ar['cohort']} | {ar['n_cells']} | {ar['counted_cells_for_acceptance']} | {ar['n_mapped']} | {frac} | {ar['best_field_used']} | {ar['source_status']} |"
        )
    summary_path = out_dir / "mapping_summary.md"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    # Checkpoint
    checkpoint_path = checkpoint_dir / "mapping_audit.checkpoint.csv"
    csv_write(checkpoint_path, ["run_id", "created_at", "source_path", "status", "last_file_written"], [{
        "run_id": run_id,
        "created_at": now_iso(),
        "source_path": str(manifest_ref),
        "status": "completed" if not blocking_issues else "blocking",
        "last_file_written": str(audit_path),
    }])

    return {
        "run_id": run_id,
        "run_root": run_root,
        "gate_ready": gate_ready,
        "mapped_fraction": mapped_fraction,
        "blocking_issues": blocking_issues,
        "mapping_parquet_path": parquet_path if gate_ready else None,
        "audit_path": audit_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Step2.1 metadata sample mapping audit")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    args = parser.parse_args()

    result = run_step2_1_mapping_audit(root_dir=args.root_dir, run_id=args.run_id)
    print(f"run_id: {result['run_id']}")
    print(f"gate_ready: {result['gate_ready']}")
    print(f"mapped_fraction: {result['mapped_fraction']:.4f}")
    if result['blocking_issues']:
        print(f"BLOCKING: {', '.join(result['blocking_issues'])}")
    else:
        print("No blocking issues")


if __name__ == "__main__":
    main()
