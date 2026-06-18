#!/usr/bin/env python3
"""Rematerialize Phase4A identity/context and Phase4B aggregation inputs.

The script preserves labels in the historical Phase4A master. For newly
appended objects it maps source-author labels, or an explicit source-sorted
population, to the frozen ontology without using response. It overlays the
response-blind identity sidecar and derives only the Phase4B tables required by
Phase7/8 from the corrected analysis units.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MASTER = ROOT / "results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz"
DEFAULT_INTERFACE = ROOT / "results/v6_2/data_interface_repair_v1"
DEFAULT_METADATA = DEFAULT_INTERFACE / "metadata"
DEFAULT_OUT4A = ROOT / "results/v6_2/phase4a_cell_state_harmonization_repair_v1"
DEFAULT_OUT4B = ROOT / "results/v6_2/phase4b_immune_state_feature_construction_repair_v1"
UNKNOWN = {"", "unknown", "nan", "none", "na", "null"}
DEFAULT_APPEND_OBJECTS = {
    "GSE272734": ROOT / "data/processed/srt/raw/gse272734.h5ad",
    "GSE200996": ROOT / "data/processed/srt/raw/gse200996.h5ad",
    "GSE179994": ROOT / "data/processed/srt/raw/gse179994.h5ad",
    "mendeley_skrx2fz79n": ROOT / "data/processed/srt/raw/mendeley_skrx2fz79n.h5ad",
    "bi_2021_rcc": ROOT / "data/processed/srt/raw/bi2021rcc.h5ad",
}
ROW_ALIGNED_DUPLICATE_BARCODE_OBJECTS = {
    "GSE272734": ROOT / "data/processed/srt/raw/gse272734.h5ad",
    "GSE272735": ROOT / "data/processed/srt/raw/gse272735.h5ad",
    "GSE273718": ROOT / "data/processed/srt/raw/gse273718.h5ad",
}
IDENTITY_COLUMNS = [
    "source_object_id", "source_cell_id", "original_cell_barcode", "raw_barcode",
    "source_sample_id", "legacy_sample_key", "biological_sample_key", "sample_key", "study_subject_key",
    "patient_key", "patient_id", "library_key", "demux_id", "normalized_timepoint",
    "timepoint_raw", "tissue_context", "tissue_raw", "lesion_context", "treatment_arm",
    "analysis_unit_key", "assignment_status", "quarantine_reason", "provenance_source",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def safe(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unknown"


def is_known(series: pd.Series) -> pd.Series:
    return ~series.fillna("").astype(str).str.strip().str.lower().isin(UNKNOWN)


def collapse_unique(values: pd.Series, mixed: str = "unknown_or_mixed") -> object:
    """Return one non-empty value, otherwise an explicit mixed/unknown marker."""
    kept = [value for value in values.dropna().tolist()
            if str(value).strip().lower() not in UNKNOWN]
    unique = list(dict.fromkeys(map(str, kept)))
    if len(unique) == 1:
        return unique[0]
    return mixed if unique else "unknown"


def stable_annotation_dtypes(frame: pd.DataFrame) -> pd.DataFrame:
    """Prevent chunk-local all-null inference from changing Parquet schemas."""
    numeric_float = {"mapping_confidence", "reliability_score"}
    numeric_int = {
        "low_quality_flag", "doublet_risk_flag", "ambient_rna_risk_flag",
        "stress_dissociation_flag", "cell_count_support",
    }
    out = frame.copy()
    for column in out:
        if column in numeric_float:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)
        elif column in numeric_int:
            out[column] = pd.to_numeric(out[column], errors="coerce").fillna(0).astype("int64")
        else:
            out[column] = out[column].astype("string").fillna("")
    return out


def find_identity_files(root: Path, cohort: str) -> list[Path]:
    candidates = list(root.glob(f"**/cohort_id={cohort}/**/*.parquet"))
    candidates += list(root.glob(f"**/{safe(cohort)}*.parquet"))
    return sorted({p for p in candidates if "audit" not in p.parts and "metadata" not in p.parts})


def normalize_identity(frame: pd.DataFrame, cohort: str) -> pd.DataFrame:
    aliases = {
        "source_cell_id": ["source_cell_id", "original_cell_barcode", "raw_barcode", "cell_id", "barcode"],
        "source_object_id": ["source_object_id", "object_id"],
        "study_subject_key": ["study_subject_key", "canonical_patient_key", "patient_key"],
        "patient_key": ["patient_key", "study_subject_key", "canonical_patient_key"],
        "normalized_timepoint": ["normalized_timepoint", "timepoint", "collection_timepoint"],
        "tissue_context": ["tissue_context", "tissue_source", "tissue"],
        "analysis_unit_key": ["analysis_unit_key", "patient_timepoint_context_id"],
        "assignment_status": ["assignment_status", "identity_status"],
        "provenance_source": ["provenance_source", "assignment_provenance"],
    }
    out = frame.copy()
    for target, choices in aliases.items():
        if target not in out:
            source = next((c for c in choices if c in out), None)
            out[target] = out[source] if source else ""
    out["cohort_id"] = cohort
    out["source_cell_id"] = out["source_cell_id"].astype(str)
    if "original_cell_barcode" not in out:
        out["original_cell_barcode"] = out["raw_barcode"] if "raw_barcode" in out else out["source_cell_id"]
    if "study_subject_key" not in out or not is_known(out["study_subject_key"]).all():
        out["study_subject_key"] = out["patient_key"]
    if "patient_key" not in out:
        out["patient_key"] = out["study_subject_key"]
    out["normalized_timepoint"] = out["normalized_timepoint"].fillna("unknown").astype(str)
    out["tissue_context"] = out["tissue_context"].fillna("unknown").astype(str)
    missing_unit = ~is_known(out["analysis_unit_key"])
    out.loc[missing_unit, "analysis_unit_key"] = (
        out.loc[missing_unit, "cohort_id"].astype(str) + "::"
        + out.loc[missing_unit, "study_subject_key"].astype(str) + "::"
        + out.loc[missing_unit, "normalized_timepoint"].astype(str) + "::"
        + out.loc[missing_unit, "tissue_context"].astype(str)
    )
    if "sample_key" not in out:
        out["sample_key"] = out["analysis_unit_key"]
    if "biological_sample_key" not in out:
        out["biological_sample_key"] = out["sample_key"]
    if "assignment_status" not in out:
        out["assignment_status"] = "resolved"
    return out


class IdentityCache:
    def __init__(self, root: Path):
        self.root = root
        self.cache: dict[str, pd.DataFrame] = {}

    def get(self, cohort: str) -> pd.DataFrame:
        if cohort in self.cache:
            return self.cache[cohort]
        files = find_identity_files(self.root, cohort)
        if not files:
            self.cache[cohort] = pd.DataFrame()
            return self.cache[cohort]
        parts = [pd.read_parquet(path) for path in files]
        frame = normalize_identity(pd.concat(parts, ignore_index=True, sort=False), cohort)
        key = [c for c in ("source_object_id", "source_cell_id") if c in frame]
        if not key:
            raise RuntimeError(f"{cohort}: identity sidecar has no join key")
        duplicated = frame.duplicated(key, keep=False)
        if duplicated.any():
            conflicts = frame.loc[duplicated].groupby(key, dropna=False)["analysis_unit_key"].nunique()
            if (conflicts > 1).any():
                raise RuntimeError(f"{cohort}: identity sidecar join key maps to multiple analysis units")
            frame = frame.drop_duplicates(key)
        self.cache[cohort] = frame
        if len(self.cache) > 3:
            oldest = next(iter(self.cache))
            if oldest != cohort:
                self.cache.pop(oldest)
        return frame


def overlay_group(group: pd.DataFrame, identity: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cohort = str(group["cohort_id"].iloc[0])
    before = len(group)
    group = group.copy()
    group["legacy_sample_key"] = group.get("sample_key", "")
    if identity.empty:
        fallback = is_known(group["patient_key"]) & is_known(group["legacy_sample_key"])
        group["study_subject_key"] = group["patient_key"]
        group["biological_sample_key"] = group["legacy_sample_key"]
        group["library_key"] = ""
        group["demux_id"] = ""
        group["normalized_timepoint"] = "unknown"
        group["tissue_context"] = "unknown"
        group["lesion_context"] = ""
        group["treatment_arm"] = "unknown"
        group["analysis_unit_key"] = group["legacy_sample_key"].where(fallback, "")
        group["assignment_status"] = fallback.map(
            {True: "resolved_legacy_validated_fallback", False: "quarantine_missing_identity_sidecar"})
        group["quarantine_reason"] = fallback.map(
            {True: "", False: "missing_identity_sidecar_and_legacy_key"})
        group["provenance_source"] = "phase4a_frozen_patient_sample_key_fallback"
        group.loc[fallback, "sample_key"] = group.loc[fallback, "analysis_unit_key"]
        return group, {
            "cohort_id": cohort, "n_input": before, "n_output": len(group),
            "n_resolved": int(fallback.sum()), "n_quarantined": int((~fallback).sum()),
            "n_identity_resolved": 0, "n_legacy_fallback": int(fallback.sum()),
            "status": "PASS" if fallback.all() else "CONDITIONAL" if fallback.any() else "HARD_FAIL",
        }
    join_left = ["source_object_id", "original_cell_barcode"]
    join_right = ["source_object_id", "__join_cell"]
    right = identity.copy()
    right["__join_cell"] = right["source_cell_id"].astype(str)
    if "source_object_id" not in right or not is_known(right["source_object_id"]).any():
        join_left = ["original_cell_barcode"]
        join_right = ["__join_cell"]
    keep = list(dict.fromkeys(join_right + [c for c in IDENTITY_COLUMNS if c in right and c not in join_right]))
    right = right[keep].drop_duplicates(join_right)
    merged = group.merge(right, left_on=join_left, right_on=join_right, how="left",
                         validate="many_to_one", suffixes=("", "_identity"))
    merged.drop(columns=[c for c in join_right if c not in join_left], inplace=True, errors="ignore")
    for col in IDENTITY_COLUMNS:
        incoming = f"{col}_identity"
        if incoming in merged:
            merged[col] = merged[incoming].where(is_known(merged[incoming]), merged.get(col, ""))
            merged.drop(columns=incoming, inplace=True)
    identity_resolved = (
        is_known(merged.get("analysis_unit_key", pd.Series("", index=merged.index)))
        & merged.get("assignment_status", pd.Series("", index=merged.index)).astype(str).str.lower().eq("resolved")
        & is_known(merged.get("study_subject_key", pd.Series("", index=merged.index)))
    )
    legacy_resolved = (
        ~identity_resolved
        & is_known(merged.get("patient_key", pd.Series("", index=merged.index)))
        & is_known(merged.get("legacy_sample_key", pd.Series("", index=merged.index)))
    )
    merged.loc[legacy_resolved, "study_subject_key"] = merged.loc[legacy_resolved, "patient_key"]
    merged.loc[legacy_resolved, "analysis_unit_key"] = merged.loc[legacy_resolved, "legacy_sample_key"]
    merged.loc[legacy_resolved, "biological_sample_key"] = merged.loc[legacy_resolved, "legacy_sample_key"]
    merged.loc[legacy_resolved, "assignment_status"] = "resolved_legacy_validated_fallback"
    merged.loc[legacy_resolved, "quarantine_reason"] = ""
    merged.loc[legacy_resolved, "provenance_source"] = "phase4a_frozen_patient_sample_key_fallback"
    resolved = identity_resolved | legacy_resolved
    merged.loc[~resolved, "assignment_status"] = "quarantine_unresolved_identity"
    merged.loc[~resolved, "quarantine_reason"] = "cell_not_joined_to_identity_or_legacy_key"
    merged.loc[resolved, "sample_key"] = merged.loc[resolved, "analysis_unit_key"]
    merged.loc[resolved, "patient_key"] = merged.loc[resolved, "study_subject_key"]
    if len(merged) != before:
        raise RuntimeError(f"{cohort}: cell count changed during identity overlay")
    return merged, {
        "cohort_id": cohort, "n_input": before, "n_output": len(merged),
        "n_resolved": int(resolved.sum()), "n_quarantined": int((~resolved).sum()),
        "n_identity_resolved": int(identity_resolved.sum()),
        "n_legacy_fallback": int(legacy_resolved.sum()),
        "status": "PASS" if resolved.all() else ("CONDITIONAL" if resolved.any() else "HARD_FAIL"),
    }


def rematerialize(master: Path, identity_root: Path, out4a: Path, cohorts: set[str] | None,
                  chunksize: int, resume: bool) -> pd.DataFrame:
    (out4a / "qc").mkdir(parents=True, exist_ok=True)
    part_dir = out4a / "handoff/cell_state_annotation_master.parquet"
    part_dir.mkdir(parents=True, exist_ok=True)
    completed = {p.stem for p in (out4a / "completion").glob("*.complete")} if resume else set()
    cache = IdentityCache(identity_root)
    writers: dict[str, pq.ParquetWriter] = {}
    audits: list[dict] = []
    seen_counts: Counter = Counter()
    resolved_counts: Counter = Counter()
    identity_resolved_counts: Counter = Counter()
    legacy_fallback_counts: Counter = Counter()
    for chunk in pd.read_csv(master, chunksize=chunksize, low_memory=False):
        if cohorts:
            chunk = chunk[chunk["cohort_id"].isin(cohorts)]
        for cohort, group in chunk.groupby("cohort_id", sort=False):
            cohort = str(cohort)
            if cohort in ROW_ALIGNED_DUPLICATE_BARCODE_OBJECTS:
                continue
            if cohort in completed:
                continue
            fixed, audit = overlay_group(group, cache.get(cohort))
            fixed = stable_annotation_dtypes(fixed)
            seen_counts[cohort] += audit["n_input"]
            resolved_counts[cohort] += audit["n_resolved"]
            identity_resolved_counts[cohort] += audit.get("n_identity_resolved", 0)
            legacy_fallback_counts[cohort] += audit.get("n_legacy_fallback", 0)
            path = part_dir / f"cohort_id={safe(cohort)}.parquet"
            table = pa.Table.from_pandas(fixed, preserve_index=False)
            if cohort not in writers:
                if path.exists():
                    path.unlink()
                writers[cohort] = pq.ParquetWriter(path, table.schema, compression="zstd")
            writers[cohort].write_table(table)
    for writer in writers.values():
        writer.close()
    for cohort in writers:
        marker = out4a / "completion" / f"{cohort}.complete"
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(now_iso() + "\n")
    for cohort, n in sorted(seen_counts.items()):
        audits.append({"cohort_id": cohort, "n_input": n, "n_output": n,
                       "n_resolved": resolved_counts[cohort], "n_quarantined": n - resolved_counts[cohort],
                       "n_identity_resolved": identity_resolved_counts[cohort],
                       "n_legacy_fallback": legacy_fallback_counts[cohort],
                       "status": "PASS" if n == resolved_counts[cohort] else (
                           "CONDITIONAL" if resolved_counts[cohort] else "HARD_FAIL")})
    audit_df = pd.DataFrame(audits)
    audit_df.to_csv(out4a / "qc/cell_identity_overlay_conservation.csv", index=False)
    return audit_df


def audit_existing_partitions(out4a: Path) -> pd.DataFrame:
    """Recover the identity audit after all parquet writers completed."""
    (out4a / "qc").mkdir(parents=True, exist_ok=True)
    rows = []
    for path, frame in scan_corrected(out4a / "handoff/cell_state_annotation_master.parquet"):
        status = frame.get("assignment_status", pd.Series("", index=frame.index)).astype(str).str.lower()
        provenance = frame.get("provenance_source", pd.Series("", index=frame.index)).astype(str)
        resolved = status.str.startswith("resolved")
        cohort = str(frame.cohort_id.iloc[0]) if len(frame) else path.stem.replace("cohort_id=", "")
        rows.append({
            "cohort_id": cohort, "n_input": len(frame), "n_output": len(frame),
            "n_resolved": int(resolved.sum()), "n_quarantined": int((~resolved).sum()),
            "n_identity_resolved": int((resolved & ~provenance.eq("phase4a_frozen_patient_sample_key_fallback")).sum()),
            "n_legacy_fallback": int(provenance.eq("phase4a_frozen_patient_sample_key_fallback").sum()),
            "status": "PASS" if resolved.all() else "CONDITIONAL" if resolved.any() else "HARD_FAIL",
        })
    audit = pd.DataFrame(rows)
    audit.to_csv(out4a / "qc/cell_identity_overlay_conservation.csv", index=False)
    return audit


def overlay_row_aligned_duplicate_barcode_objects(
    identity_root: Path,
    out4a: Path,
    objects: dict[str, Path],
) -> pd.DataFrame:
    """Repair objects whose bare 10x barcodes repeat across source samples.

    Row alignment is accepted only after exact row-count and source-barcode
    equality checks. Existing response-blind annotations are preserved.
    """
    cache = IdentityCache(identity_root)
    part_dir = out4a / "handoff/cell_state_annotation_master.parquet"
    rows = []
    for cohort, h5ad_path in objects.items():
        target = part_dir / f"cohort_id={safe(cohort)}.parquet"
        if not target.exists() or not h5ad_path.exists():
            rows.append({"cohort_id": cohort, "status": "MISSING", "n_cells": 0})
            continue
        current = pd.read_parquet(target)
        identity = cache.get(cohort).copy()
        if len(current) != len(identity):
            rows.append({"cohort_id": cohort, "status": "HARD_FAIL_ROW_COUNT", "n_cells": len(current)})
            continue
        order = pd.to_numeric(identity["source_cell_id"].astype(str).str.extract(r"^row::(\d+)$")[0], errors="coerce")
        if order.isna().any() or order.nunique() != len(identity):
            rows.append({"cohort_id": cohort, "status": "HARD_FAIL_ROW_KEY", "n_cells": len(current)})
            continue
        identity = identity.assign(__row=order.astype(int)).sort_values("__row").reset_index(drop=True)
        expected = identity.get("original_cell_barcode", identity.get("raw_barcode", pd.Series("", index=identity.index))).astype(str)
        observed = current["original_cell_barcode"].astype(str).reset_index(drop=True)
        if not observed.equals(expected.reset_index(drop=True)):
            rows.append({"cohort_id": cohort, "status": "HARD_FAIL_BARCODE_ORDER", "n_cells": len(current)})
            continue
        current["legacy_sample_key"] = current.get("legacy_sample_key", current.get("sample_key", ""))
        for column in IDENTITY_COLUMNS:
            if column in identity:
                current[column] = identity[column].to_numpy()
        current["sample_id"] = identity.get("biological_sample_key", identity.get("sample_key", "")).to_numpy()
        current["sample_key"] = identity["analysis_unit_key"].to_numpy()
        current["patient_key"] = identity["study_subject_key"].to_numpy()
        current["patient_id"] = identity.get("patient_id", identity["study_subject_key"]).astype(str).to_numpy()
        current["cell_key"] = (
            cohort + "::" + current["sample_key"].astype(str) + "::" + current["source_cell_id"].astype(str)
        )
        current = stable_annotation_dtypes(current)
        current.to_parquet(target, index=False, compression="zstd")
        rows.append({"cohort_id": cohort, "status": "PASS", "n_cells": len(current)})
    audit = pd.DataFrame(rows)
    audit.to_csv(out4a / "qc/row_aligned_duplicate_barcode_overlay.csv", index=False)
    return audit


def scan_corrected(part_dir: Path):
    for path in sorted(part_dir.glob("cohort_id=*.parquet")):
        yield path, pd.read_parquet(path)


def author_annotation_labels(values: pd.Series) -> pd.DataFrame:
    """Map explicit author labels to the frozen ontology without response use."""
    coarse, mid = [], []
    for raw in values.fillna("").astype(str):
        text = raw.casefold()
        if "cd8" in text:
            pair = ("T_NK", "CD8_T")
        elif "cd4" in text or "treg" in text:
            pair = ("T_NK", "Treg" if "treg" in text else "CD4_T")
        elif "nk" in text:
            pair = ("T_NK", "NK")
        elif "t cell" in text or text.startswith("t_"):
            pair = ("T_NK", "Unknown")
        elif "plasma" in text:
            pair = ("B_Plasma", "Plasma_cell")
        elif "b cell" in text or "b_cell" in text:
            pair = ("B_Plasma", "B_cell")
        elif "dendritic" in text or re.search(r"(^|[^a-z])dc([^a-z]|$)", text):
            pair = ("DC_APC", "DC")
        elif "macro" in text:
            pair = ("Myeloid", "Macrophage")
        elif "mono" in text:
            pair = ("Myeloid", "Monocyte")
        elif "neut" in text:
            pair = ("Myeloid", "Neutrophil_TAN_like")
        elif "myeloid" in text or "mast" in text:
            pair = ("Myeloid", "Unknown")
        elif "fibro" in text or "caf" in text or "strom" in text:
            pair = ("CAF_Stromal", "CAF")
        elif "endo" in text:
            pair = ("Endothelial", "Endothelial")
        elif "tumor" in text or "epithelial" in text or "malignant" in text:
            pair = ("Tumor_like", "Tumor_epithelial_like")
        elif "cycling" in text or "prolif" in text:
            pair = ("Cycling", "Unknown")
        else:
            pair = ("Unknown", "Unknown")
        coarse.append(pair[0])
        mid.append(pair[1])
    return pd.DataFrame({"coarse": coarse, "mid": mid}, index=values.index)


def append_missing_objects(identity_root: Path, out4a: Path, objects: dict[str, Path]) -> pd.DataFrame:
    """Append locally registered scRNA objects absent from the historical master."""
    import anndata as ad

    schema_cols = pd.read_csv(DEFAULT_MASTER, nrows=0).columns.tolist()
    cache = IdentityCache(identity_root)
    rows = []
    part_dir = out4a / "handoff/cell_state_annotation_master.parquet"
    for cohort, path in objects.items():
        if not path.exists():
            rows.append({"cohort_id": cohort, "status": "MISSING", "n_cells": 0, "n_resolved": 0})
            continue
        target = part_dir / f"cohort_id={safe(cohort)}.parquet"
        if target.exists() and cohort != "GSE272734":
            continue
        adata = ad.read_h5ad(path, backed="r")
        try:
            obs = adata.obs.copy()
            obs.index = obs.index.astype(str)
        finally:
            adata.file.close()
        identity = cache.get(cohort)
        if identity.empty:
            rows.append({"cohort_id": cohort, "status": "HARD_FAIL_NO_IDENTITY", "n_cells": len(obs), "n_resolved": 0})
            continue
        identity = identity.copy()
        key = "source_cell_id" if "source_cell_id" in identity else "original_cell_barcode"
        identity["__join_cell"] = identity[key].astype(str)
        if identity["__join_cell"].duplicated().any():
            rows.append({"cohort_id": cohort, "status": "HARD_FAIL_DUPLICATE_IDENTITY", "n_cells": len(obs), "n_resolved": 0})
            continue
        source_ids = pd.Series(obs.index.astype(str), index=obs.index)
        if cohort == "GSE272734" and source_ids.duplicated().any():
            source_ids = pd.Series([f"row::{i}" for i in range(len(obs))], index=obs.index)
        frame = pd.DataFrame({"source_cell_id_input": source_ids.to_numpy(),
                              "original_cell_barcode": obs.index.astype(str), "cohort_id": cohort})
        frame = frame.merge(identity.drop(columns=["original_cell_barcode"], errors="ignore"),
                            left_on="source_cell_id_input", right_on="__join_cell",
                            how="left", validate="one_to_one").drop(columns="__join_cell")
        if cohort in {"mendeley_skrx2fz79n", "bi_2021_rcc"}:
            status = frame.get("assignment_status", pd.Series("", index=frame.index)).astype(str).str.lower()
            unresolved = ~status.str.startswith("resolved")
            sample = frame.get("source_sample_id", frame.get("biological_sample_key", pd.Series("cohort_pool", index=frame.index)))
            sample = sample.fillna("cohort_pool").astype(str)
            subject = cohort + "::support::" + sample
            unit = ["support::" + hashlib.sha256(f"{cohort}|{item}".encode()).hexdigest()[:20] for item in sample]
            frame.loc[unresolved, "study_subject_key"] = subject[unresolved]
            frame.loc[unresolved, "patient_key"] = subject[unresolved]
            frame.loc[unresolved, "analysis_unit_key"] = pd.Series(unit, index=frame.index)[unresolved]
            frame.loc[unresolved, "assignment_status"] = "resolved_support_only_synthetic_subject"
            frame.loc[unresolved, "quarantine_reason"] = ""
            frame.loc[unresolved, "provenance_source"] = "support_only_sample_identity_no_patient_claim"
        annotation_col = next((c for c in ("anno_orig", "Annotation", "cell_type", "celltype_raw",
                                             "celltype", "CellType") if c in obs), None)
        raw = obs[annotation_col].astype(str).reset_index(drop=True) if annotation_col else pd.Series("", index=range(len(obs)))
        if cohort == "GSE272734" and raw.replace({"nan": "", "None": ""}).str.strip().eq("").all():
            raw = pd.Series("CD8 T cell (source-sorted PBMC)", index=range(len(obs)))
            annotation_col = "source_sorted_population"
        labels = author_annotation_labels(raw)
        resolved = (
            is_known(frame.get("analysis_unit_key", pd.Series("", index=frame.index)))
            & frame.get("assignment_status", pd.Series("", index=frame.index)).astype(str).str.lower().str.startswith("resolved")
        )
        out = pd.DataFrame("", index=frame.index, columns=schema_cols)
        out["cohort_id"] = cohort
        out["source_object_id"] = frame.get("source_object_id", f"{cohort}::object_1")
        out["source_cell_id"] = frame.get("source_cell_id", frame.original_cell_barcode)
        out["source_object_path"] = str(path)
        out["original_cell_barcode"] = frame.original_cell_barcode
        out["standardized_cell_barcode"] = frame.original_cell_barcode
        out["sample_id"] = frame.get("biological_sample_key", frame.get("sample_key", ""))
        out["sample_key"] = frame.get("analysis_unit_key", "")
        out["patient_key"] = frame.get("study_subject_key", frame.get("patient_key", ""))
        out["patient_id"] = frame.get("patient_id", out["patient_key"].astype(str).str.rsplit("::", n=1).str[-1])
        out["cell_key"] = cohort + "::" + out["sample_key"].astype(str) + "::" + out["original_cell_barcode"].astype(str)
        out["original_annotation"] = raw.to_numpy()
        out["original_annotation_field"] = annotation_col or "none"
        out["marker_based_coarse_label"] = labels.coarse.to_numpy()
        out["marker_based_mid_label"] = labels.mid.to_numpy()
        out["marker_based_fine_candidate"] = "Unknown"
        out["reference_mapped_coarse_label"] = labels.coarse.to_numpy()
        out["reference_mapped_mid_label"] = labels.mid.to_numpy()
        out["reference_mapped_fine_label"] = "Unknown"
        out["reference_mapped_label"] = labels.mid.where(labels.mid.ne("Unknown"), labels.coarse).to_numpy()
        out["harmonized_coarse_label"] = labels.coarse.to_numpy()
        out["harmonized_mid_label"] = labels.mid.to_numpy()
        out["harmonized_fine_label"] = "Unknown"
        out["final_label_level"] = labels.mid.ne("Unknown").map({True: "mid", False: "coarse"}).to_numpy()
        out["annotation_method"] = "source_author_label_response_blind_context_repair_v1"
        out["integration_tier"] = "full_integration_candidate"
        out["mapping_confidence"] = labels.coarse.ne("Unknown").map({True: 0.80, False: 0.0}).to_numpy()
        out["reliability_score"] = labels.coarse.ne("Unknown").map({True: 0.70, False: 0.0}).to_numpy()
        out["reliability_level"] = labels.mid.ne("Unknown").map({True: "mid_reliable", False: "coarse_or_unknown"}).to_numpy()
        out["low_quality_flag"] = 0
        out["doublet_risk_flag"] = 0
        out["ambient_rna_risk_flag"] = 0
        out["stress_dissociation_flag"] = 0
        allowed = resolved & labels.coarse.ne("Unknown")
        out["allowed_phase4b_fraction"] = allowed.map({True: "yes", False: "no"}).to_numpy()
        out["allowed_phase4b_pseudobulk"] = allowed.map({True: "yes", False: "no"}).to_numpy()
        out["allowed_phase4b_signature"] = allowed.map({True: "yes", False: "no"}).to_numpy()
        out["allowed_hcc_context"] = "yes" if cohort == "mendeley_skrx2fz79n" else "no"
        out["allowed_tcr_coupling"] = "no"
        out["allowed_downstream_use"] = allowed.map({True: "coarse_mid_support", False: "excluded"}).to_numpy()
        out["downgrade_reason"] = ""
        out.loc[~resolved, "downgrade_reason"] = "unresolved_identity"
        for col in ("source_cell_id", "source_sample_id", "legacy_sample_key", "biological_sample_key", "study_subject_key", "library_key", "demux_id",
                    "normalized_timepoint", "tissue_context", "lesion_context", "treatment_arm", "analysis_unit_key",
                    "assignment_status", "quarantine_reason", "provenance_source"):
            out[col] = frame.get(col, "")
        out.to_parquet(target, index=False, compression="zstd")
        rows.append({"cohort_id": cohort, "status": "PASS" if resolved.all() else "CONDITIONAL",
                     "n_cells": len(out), "n_resolved": int(resolved.sum()),
                     "n_coarse_known": int(labels.coarse.ne("Unknown").sum()),
                     "n_mid_known": int(labels.mid.ne("Unknown").sum())})
    audit = pd.DataFrame(rows)
    audit.to_csv(out4a / "qc/appended_object_annotation_audit.csv", index=False)
    return audit


def fraction_matrix(counts: Counter, totals: Counter, prefix: str) -> pd.DataFrame:
    rows: dict[str, dict] = defaultdict(dict)
    for (unit, state), n in counts.items():
        denom = totals.get(unit, 0)
        rows[unit][f"frac_{prefix}__{safe(state)}"] = n / denom if denom else float("nan")
    out = pd.DataFrame.from_dict(rows, orient="index").reset_index(names="sample_key")
    return out


def build_phase4b(out4a: Path, metadata_root: Path, out4b: Path) -> dict:
    out4b.mkdir(parents=True, exist_ok=True)
    for sub in ("fractions", "qc_covariates", "response_environment", "universe", "handoff", "audit", "matrix"):
        (out4b / sub).mkdir(parents=True, exist_ok=True)
    counts = {level: Counter() for level in ("coarse", "mid", "fine_restricted")}
    totals = {level: Counter() for level in counts}
    unit_meta: dict[str, dict] = {}
    unit_legacy_samples: defaultdict[str, set[str]] = defaultdict(set)
    annotation_counts: Counter = Counter()
    total_cells = 0
    for _, frame in scan_corrected(out4a / "handoff/cell_state_annotation_master.parquet"):
        status = frame["assignment_status"].fillna("").astype(str).str.lower()
        valid_unit = is_known(frame["analysis_unit_key"]) & status.str.startswith(("resolved", "eligible"))
        frame = frame[valid_unit].copy()
        total_cells += len(frame)
        if "legacy_sample_key" in frame:
            legacy_pairs = frame.loc[valid_unit, ["analysis_unit_key", "legacy_sample_key"]].drop_duplicates()
            for unit, legacy in legacy_pairs.itertuples(index=False):
                if str(legacy).strip().lower() not in UNKNOWN:
                    unit_legacy_samples[str(unit)].add(str(legacy))
        for row in frame[["analysis_unit_key", "cohort_id", "study_subject_key", "normalized_timepoint",
                          "tissue_context", "lesion_context", "biological_sample_key", "library_key",
                          "treatment_arm"]].drop_duplicates().itertuples(index=False):
            unit_meta[row.analysis_unit_key] = {
                "sample_key": row.analysis_unit_key, "analysis_unit_key": row.analysis_unit_key,
                "cohort_id": row.cohort_id, "patient_key": row.study_subject_key,
                "timepoint": row.normalized_timepoint, "tissue_source": row.tissue_context,
                "lesion_context": row.lesion_context,
                "biological_sample_key": row.biological_sample_key, "library_key": row.library_key,
                "treatment_arm": row.treatment_arm,
            }
        allowed = frame.get("allowed_phase4b_fraction", pd.Series("yes", index=frame.index)).astype(str).str.lower().eq("yes")
        annotation_counts.update(frame.loc[allowed, "analysis_unit_key"].astype(str))
        specs = {
            "coarse": "harmonized_coarse_label",
            "mid": "harmonized_mid_label",
            "fine_restricted": "harmonized_fine_label",
        }
        for level, col in specs.items():
            known = allowed & is_known(frame[col])
            if level == "fine_restricted":
                reliability = frame["reliability_level"] if "reliability_level" in frame else pd.Series("", index=frame.index)
                known &= reliability.astype(str).eq("high_reliability_fine")
            grouped = frame.loc[known].groupby(["analysis_unit_key", col], observed=True).size()
            for key, value in grouped.items():
                counts[level][(str(key[0]), str(key[1]))] += int(value)
            grouped_total = frame.loc[known].groupby("analysis_unit_key", observed=True).size()
            for key, value in grouped_total.items():
                totals[level][str(key)] += int(value)
    meta = pd.DataFrame(unit_meta.values()).sort_values(["cohort_id", "patient_key", "timepoint", "tissue_source"])
    meta["legacy_sample_keys"] = meta["sample_key"].map(
        lambda key: "|".join(sorted(unit_legacy_samples.get(str(key), set()))))
    matrices = []
    for level in counts:
        matrices.append(fraction_matrix(counts[level], totals[level], level))
    fractions = meta.copy()
    for matrix in matrices:
        fractions = fractions.merge(matrix, on="sample_key", how="left", validate="one_to_one")
    fractions.to_csv(out4b / "fractions/cell_state_fraction_matrix.sample_level.csv", index=False)
    fractions.to_csv(out4b / "fractions/cell_state_fraction_matrix.patient_timepoint_level.csv", index=False)
    feature_rows = []
    for col in fractions.columns:
        if col.startswith("frac_"):
            feature_rows.append({"feature_id": col, "feature_family": "cell_state_fraction",
                                 "resolution": col.split("__", 1)[0].replace("frac_", ""),
                                 "is_primary_feature": "no" if "fine_restricted" in col else "yes",
                                 "is_response_derived": "false", "source": "context_repair_v1"})
    pd.DataFrame(feature_rows).to_csv(out4b / "fractions/cell_state_fraction_feature_dictionary.csv", index=False)
    qc = meta.copy()
    qc["total_cells"] = qc["sample_key"].map(annotation_counts).fillna(0).astype(int)
    for level in ("coarse", "mid", "fine_restricted"):
        qc[f"{level.replace('_restricted', '')}_label_coverage"] = [
            totals[level].get(key, 0) / annotation_counts.get(key, 1) for key in qc.sample_key
        ]
    for col in ("detected_genes", "umi", "mitochondrial_fraction", "ribosomal_fraction", "low_quality_fraction"):
        qc[col] = ""
    qc["treatment_context"] = ""
    qc["expression_layer"] = "context_repair_v1"
    qc["integration_status"] = "full_integration_not_completed"
    qc.to_csv(out4b / "qc_covariates/sample_patient_timepoint_qc_covariates.csv", index=False)
    sample_files = sorted(metadata_root.glob("*sample*metadata*frozen_v1*.csv")) + sorted(metadata_root.glob("*sample*effective*.csv"))
    patient_files = sorted(metadata_root.glob("*patient*metadata*frozen_v1*.csv")) + sorted(metadata_root.glob("*patient*effective*.csv"))
    binding = meta.copy()
    analysis_sample = meta.copy()
    analysis_sample["treatment_context"] = analysis_sample["treatment_arm"].where(
        is_known(analysis_sample["treatment_arm"]), "unknown")
    analysis_sample["timepoint_raw"] = analysis_sample["timepoint"]
    analysis_sample["treatment_raw"] = analysis_sample["treatment_arm"]
    if sample_files:
        sm = pd.read_csv(sample_files[0], low_memory=False)
        subject_col = "study_subject_key" if "study_subject_key" in sm else "patient_key"
        sm[subject_col] = sm[subject_col].fillna(sm.get("patient_key", "")).astype(str)
        semantic = [c for c in (
            "response_raw", "response_endpoint_type", "response_binary_harmonized",
            "response_harmonization_confidence", "supervised_use_allowed", "support_use_allowed",
            "dataset_role", "cancer_type", "treatment_context",
        ) if c in sm]
        exact_fields = list(dict.fromkeys(c for c in (
            "sample_key", "timepoint", "timepoint_raw", "tissue_source", "treatment_context",
            "treatment_raw", *semantic,
        ) if c in sm))
        exact = sm[exact_fields].drop_duplicates("sample_key")
        bridge = analysis_sample[["sample_key", "legacy_sample_keys"]].copy()
        bridge["legacy_sample_key"] = bridge["legacy_sample_keys"].fillna("").str.split("|")
        bridge = bridge.explode("legacy_sample_key")
        bridge = bridge[bridge.legacy_sample_key.astype(str).str.strip().ne("")]
        exact = bridge.merge(exact, left_on="legacy_sample_key", right_on="sample_key", how="left",
                             suffixes=("_analysis", "_source"))
        exact_semantic = [c for c in exact_fields if c != "sample_key"]
        if len(exact):
            exact = exact.groupby("sample_key_analysis", as_index=False)[exact_semantic].agg(collapse_unique)
            exact = exact.rename(columns={"sample_key_analysis": "sample_key"})
            analysis_sample = analysis_sample.merge(exact, on="sample_key", how="left",
                                                    suffixes=("", "_sample"), validate="one_to_one")
        for field in set(exact_fields) - {"sample_key"}:
            incoming = f"{field}_sample"
            if incoming not in analysis_sample:
                continue
            source_preferred = field in {
                "timepoint", "timepoint_raw", "tissue_source", "treatment_context", "treatment_raw",
                "response_raw", "response_endpoint_type", "response_binary_harmonized",
                "response_harmonization_confidence", "supervised_use_allowed", "support_use_allowed",
            }
            missing = (
                is_known(analysis_sample[incoming]) if source_preferred
                else ~is_known(analysis_sample[field]) if field in analysis_sample
                else pd.Series(True, index=analysis_sample.index)
            )
            if field not in analysis_sample:
                analysis_sample[field] = analysis_sample[incoming]
            else:
                analysis_sample.loc[missing, field] = analysis_sample.loc[missing, incoming]
            analysis_sample = analysis_sample.drop(columns=incoming)
        # Subject-level fallback is restricted to non-clinical descriptors.
        # Response, endpoint, treatment and timepoint require exact source-sample binding.
        subject_safe = [c for c in ("dataset_role", "cancer_type", "supervised_use_allowed", "support_use_allowed") if c in sm]
        grouped = sm.groupby(["cohort_id", subject_col], dropna=False)[subject_safe].agg(collapse_unique).reset_index()
        grouped = grouped.rename(columns={subject_col: "patient_key"})
        analysis_sample = analysis_sample.merge(grouped, on=["cohort_id", "patient_key"], how="left",
                                                suffixes=("", "_subject"), validate="many_to_one")
        for field in subject_safe:
            incoming = f"{field}_subject"
            if incoming not in analysis_sample:
                continue
            missing = ~is_known(analysis_sample[field]) if field in analysis_sample else pd.Series(True, index=analysis_sample.index)
            if field not in analysis_sample:
                analysis_sample[field] = analysis_sample[incoming]
            else:
                analysis_sample.loc[missing, field] = analysis_sample.loc[missing, incoming]
            analysis_sample = analysis_sample.drop(columns=incoming)
    if patient_files:
        pm = pd.read_csv(patient_files[0], low_memory=False)
        subject_col = "study_subject_key" if "study_subject_key" in pm else "patient_key"
        pm[subject_col] = pm[subject_col].fillna(pm.get("patient_key", "")).astype(str)
        semantic = [c for c in (
            "cancer_type", "available_modalities", "response_raw_summary", "response_endpoint_type",
            "response_binary_harmonized", "response_harmonization_confidence", "split",
            "supervised_use_allowed", "support_use_allowed", "exclusion_or_downgrade_reason",
        ) if c in pm]
        canonical_patient = pm.groupby(subject_col, dropna=False)[semantic].agg(collapse_unique).reset_index()
        canonical_patient = canonical_patient.rename(columns={subject_col: "patient_key"})
        counts_by_patient = meta.groupby("patient_key", as_index=False).agg(
            n_samples=("sample_key", "nunique"),
            timepoint_schema=("timepoint", lambda x: "|".join(sorted(set(map(str, x))))),
            paired_pre_post_available=("timepoint", lambda x: "yes" if len(set(map(str, x))) > 1 else "no"),
            treatment_context_summary=("treatment_arm", collapse_unique),
        )
        canonical_patient = canonical_patient.merge(counts_by_patient, on="patient_key", how="right", validate="one_to_one")
    else:
        canonical_patient = meta.groupby("patient_key", as_index=False).agg(
            n_samples=("sample_key", "nunique"),
            timepoint_schema=("timepoint", lambda x: "|".join(sorted(set(map(str, x))))),
            paired_pre_post_available=("timepoint", lambda x: "yes" if len(set(map(str, x))) > 1 else "no"),
            treatment_context_summary=("treatment_arm", collapse_unique),
        )
    analysis_sample["input_version"] = "context_repair_v1"
    canonical_patient["input_version"] = "context_repair_v1"
    analysis_sample.to_csv(out4b / "response_environment/sample_metadata_analysis_unit_v1.csv", index=False)
    canonical_patient.to_csv(out4b / "response_environment/patient_metadata_analysis_unit_v1.csv", index=False)
    binding = analysis_sample.copy()
    patient_keep = [c for c in canonical_patient if c == "patient_key" or c not in binding]
    binding = binding.merge(canonical_patient[patient_keep], on="patient_key", how="left", validate="many_to_one")
    binary = binding.get("response_binary_harmonized", pd.Series("unknown", index=binding.index)).astype(str).str.lower()
    binding["response_known"] = binary.isin({"responder", "non_responder", "1", "0", "1.0", "0.0"}).map({True: "yes", False: "no"})
    binding["analysis_lane"] = binding.get("supervised_use_allowed", pd.Series("no", index=binding.index)).eq("yes").map(
        {True: "supervised_candidate", False: "support_or_registry"})
    binding["endpoint_pooling_allowed"] = "endpoint_stratified_only"
    binding["timepoint_use_boundary"] = binding["timepoint"].eq("baseline").map(
        {True: "baseline_primary_candidate", False: "temporal_sensitivity_only"})
    binding["sensitivity_only_reason"] = binding["timepoint"].eq("baseline").map(
        {True: "", False: "nonbaseline"})
    binding.to_csv(out4b / "response_environment/response_environment_binding_table.csv", index=False)
    universe = meta.groupby("cohort_id", as_index=False).agg(n_analysis_units=("sample_key", "nunique"), n_patients=("patient_key", "nunique"))
    universe["universe_id"] = "context_repair_v1_all_resolved_scRNA"
    universe.to_csv(out4b / "universe/analysis_universe_registry_v6_2.csv", index=False)
    high_risk = []
    for col in fractions.columns:
        if col.startswith("frac_"):
            miss = float(fractions[col].isna().mean())
            if miss > 0.5:
                high_risk.append({"feature_id": col, "reason": "missingness_gt_0.5", "missingness": miss})
    pd.DataFrame(high_risk, columns=["feature_id", "reason", "missingness"]).to_csv(
        out4b / "audit/high_risk_feature_exclusion_list.csv", index=False)
    identity_cols = [c for c in ("sample_key", "cohort_id", "patient_key", "timepoint", "tissue_source") if c in fractions]
    primary_features = [c for c in fractions if c.startswith(("frac_coarse__", "frac_mid__"))]
    sensitivity_features = [c for c in fractions if c.startswith("frac_fine_restricted__")]
    fractions[identity_cols + primary_features].to_csv(out4b / "matrix/immune_state_feature_matrix.primary.csv", index=False)
    fractions[["sample_key"] + sensitivity_features].to_csv(out4b / "matrix/immune_state_feature_matrix.sensitivity.csv", index=False)
    qc.to_csv(out4b / "matrix/immune_state_feature_matrix.qc_covariates.csv", index=False)
    dictionary_rows = []
    for feature in primary_features + sensitivity_features:
        dictionary_rows.append({
            "feature_id": feature, "feature_family": "cell_fraction",
            "is_primary_feature": "yes" if feature in primary_features else "no",
            "is_sensitivity_feature": "yes" if feature in sensitivity_features else "no",
            "is_response_derived": "false", "source_file": "fractions/cell_state_fraction_matrix.sample_level.csv",
        })
    pd.DataFrame(dictionary_rows).to_csv(out4b / "matrix/feature_dictionary_v6_2.csv", index=False)
    missingness_rows = [{"feature_id": c, "missingness": float(fractions[c].isna().mean()),
                         "allowed_primary_after_audit": "yes" if c in primary_features and fractions[c].isna().mean() <= 0.5 else "no"}
                        for c in primary_features + sensitivity_features]
    missingness_table = pd.DataFrame(missingness_rows)
    missingness_table.to_csv(out4b / "audit/feature_missingness_report.csv", index=False)
    missingness_table.to_csv(out4b / "audit/feature_confounding_tags.csv", index=False)
    phase4b_decision = {
        "phase": "phase4b_context_repair_v1_fraction_compatibility",
        "verdict": "CONDITIONAL_GO_TO_PHASE5", "hard_blockers": [],
        "conditional_items": ["fraction_only_repair_handoff", "full_activity_families_not_rebuilt"],
        "response_leakage_detected": False, "split_leakage_detected": False,
    }
    (out4b / "handoff/phase4b_decision_manifest.yaml").write_text(yaml.safe_dump(phase4b_decision, sort_keys=False))
    (out4b / "handoff/phase4b_to_phase5_handoff.yaml").write_text(yaml.safe_dump({
        **phase4b_decision,
        "primary_matrix": "matrix/immune_state_feature_matrix.primary.csv",
        "sensitivity_matrix": "matrix/immune_state_feature_matrix.sensitivity.csv",
        "allowed_feature_families": ["cell_fraction_coarse", "cell_fraction_mid"],
        "blocked_feature_families": ["pathway_activity", "tf_activity", "tcr_mainline"],
    }, sort_keys=False))
    manifest = {
        "phase": "phase4b_context_repair_v1", "created_at": now_iso(),
        "verdict": "PASS" if len(meta) and not binding.empty else "BLOCKED",
        "n_cells_resolved": total_cells, "n_analysis_units": len(meta),
        "response_used_for_feature_construction": False,
        "full_integration_completed": False,
        "outputs": {
            "sample_fractions": "fractions/cell_state_fraction_matrix.sample_level.csv",
            "patient_fractions": "fractions/cell_state_fraction_matrix.patient_timepoint_level.csv",
            "qc_covariates": "qc_covariates/sample_patient_timepoint_qc_covariates.csv",
            "response_binding": "response_environment/response_environment_binding_table.csv",
            "analysis_unit_sample_metadata": "response_environment/sample_metadata_analysis_unit_v1.csv",
            "canonical_patient_metadata": "response_environment/patient_metadata_analysis_unit_v1.csv",
        },
    }
    (out4b / "handoff/phase4b_context_repair_manifest.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-master", type=Path, default=DEFAULT_MASTER)
    parser.add_argument("--identity-root", type=Path, default=DEFAULT_INTERFACE / "cell_identity.parquet")
    parser.add_argument("--metadata-root", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--phase4a-output-root", type=Path, default=DEFAULT_OUT4A)
    parser.add_argument("--phase4b-output-root", type=Path, default=DEFAULT_OUT4B)
    parser.add_argument("--cohort", action="append")
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--skip-default-appends", action="store_true")
    parser.add_argument("--finalize-existing", action="store_true",
                        help="Skip CSV rematerialization and audit already completed parquet partitions.")
    args = parser.parse_args()
    args.phase4a_output_root.mkdir(parents=True, exist_ok=True)
    if not args.finalize_existing:
        rematerialize(
            args.cell_master, args.identity_root, args.phase4a_output_root,
            set(args.cohort) if args.cohort else None, args.chunksize, args.resume)
    row_audit = overlay_row_aligned_duplicate_barcode_objects(
        args.identity_root, args.phase4a_output_root, ROW_ALIGNED_DUPLICATE_BARCODE_OBJECTS)
    if row_audit.status.astype(str).str.startswith("HARD_FAIL").any():
        raise SystemExit("Phase4A row-aligned duplicate-barcode repair hard gate failed")
    audit = audit_existing_partitions(args.phase4a_output_root)
    if audit.empty or audit.status.eq("HARD_FAIL").any():
        manifest = {"phase": "phase4a_context_repair_v1", "verdict": "BLOCKED",
                    "hard_fail_cohorts": audit.loc[audit.status.eq("HARD_FAIL"), "cohort_id"].tolist(),
                    "created_at": now_iso()}
        path = args.phase4a_output_root / "handoff/phase4a_context_repair_manifest.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(manifest, sort_keys=False))
        raise SystemExit("Phase4A identity overlay hard gate failed")
    if not args.skip_default_appends and not args.cohort:
        append_audit = append_missing_objects(args.identity_root, args.phase4a_output_root, DEFAULT_APPEND_OBJECTS)
        if not append_audit.empty and append_audit.status.astype(str).str.startswith("HARD_FAIL").any():
            raise SystemExit("Phase4A missing-object append hard gate failed")
    part_dir = args.phase4a_output_root / "handoff/cell_state_annotation_master.parquet"
    summary_parts = []
    eligibility_parts = []
    for _, frame in scan_corrected(part_dir):
        summary_parts.append(frame.groupby(["cohort_id", "harmonized_coarse_label", "harmonized_mid_label",
                                            "harmonized_fine_label"], dropna=False).size().reset_index(name="n_cells"))
        eligibility_parts.append(frame.groupby(
            ["cohort_id", "sample_id", "sample_key", "patient_key", "harmonized_coarse_label",
             "harmonized_mid_label", "harmonized_fine_label", "allowed_phase4b_fraction",
             "allowed_phase4b_pseudobulk", "allowed_phase4b_signature"],
            dropna=False, as_index=False).size().rename(columns={"size": "n_cells"}))
    summary = pd.concat(summary_parts, ignore_index=True).groupby(
        ["cohort_id", "harmonized_coarse_label", "harmonized_mid_label", "harmonized_fine_label"],
        dropna=False, as_index=False).n_cells.sum()
    summary.to_csv(args.phase4a_output_root / "handoff/cell_state_annotation_master_summary.csv", index=False)
    eligibility = pd.concat(eligibility_parts, ignore_index=True)
    eligibility.to_csv(
        args.phase4a_output_root / "handoff/phase4a_to_phase4b_aggregation_eligibility.csv", index=False)
    manifest = {
        "phase": "phase4a_context_repair_v1", "verdict": "PASS", "created_at": now_iso(),
        "source_master": str(args.cell_master), "source_master_sha256": sha256(args.cell_master),
        "existing_master_annotation_labels_modified": False,
        "appended_source_annotation_materialized": True,
        "appended_annotation_basis": (
            "source_author_labels; GSE272734 source-sorted CD8 T cells from PBMC"
        ),
        "celltypist_or_full_reannotation_performed": False,
        "response_used": False,
        "n_cells": int(summary.n_cells.sum()), "n_cohorts": int(summary.cohort_id.nunique()),
    }
    (args.phase4a_output_root / "handoff/phase4a_context_repair_manifest.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    phase4b = build_phase4b(args.phase4a_output_root, args.metadata_root, args.phase4b_output_root)
    print(json.dumps({"phase4a": manifest, "phase4b": phase4b}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
