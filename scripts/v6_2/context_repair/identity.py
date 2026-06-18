"""Obs-first identity construction for context_repair_v1.

This module never opens an expression matrix or consumes response labels.  It
only converts a source object's ``obs`` table and cell index into stable,
auditable identity sidecars.
"""

from __future__ import annotations

import hashlib
import re
from typing import Callable

import pandas as pd


IDENTITY_VERSION = "context_repair_v1"
UNKNOWN = "unknown"
CELL_COLUMNS = [
    "source_object_id", "source_cell_id", "raw_barcode", "cohort_id",
    "source_sample_id", "patient_id", "treatment_arm",
    "library_key", "demux_id", "biological_sample_key", "study_subject_key",
    "normalized_timepoint", "tissue_context", "lesion_context",
    "analysis_unit_key", "assignment_status", "assignment_provenance",
    "quarantine_reason",
]


def token(value: object, default: str = UNKNOWN) -> str:
    if value is None or pd.isna(value):
        return default
    value = str(value).strip()
    return value if value and value.lower() not in {"nan", "none", "na"} else default


def value(frame: pd.DataFrame, *columns: str, default: str = UNKNOWN) -> pd.Series:
    out = pd.Series(default, index=frame.index, dtype="object")
    for column in columns:
        if column in frame:
            candidate = frame[column].map(token)
            out = out.mask(out.eq(default) & candidate.ne(default), candidate)
    return out


def key(prefix: str, *parts: object) -> str:
    payload = "|".join(token(part) for part in parts)
    return f"{prefix}::{hashlib.sha256(payload.encode()).hexdigest()[:20]}"


def normalize_timepoint(raw: object) -> str:
    text = token(raw).lower().replace("_", "-").replace(" ", "-")
    if text in {"pre", "baseline", "pre-tx", "pretreatment", "pre-treatment"}:
        return "baseline"
    if text in {"post", "post-treatment", "post-tx"}:
        return "post_treatment"
    if text in {"on-treatment", "on-treatment", "on-tx", "treatment"} or text.startswith("w"):
        return "on_treatment"
    follow = re.search(r"follow-?up-?(\d+)", text)
    if follow:
        return f"on_treatment_followup_{follow.group(1)}"
    return token(raw)


def canonical_272993(patient: str) -> str:
    return re.sub(r"_(?:actla-?4|apd1)$", "", patient, flags=re.I)


def parse_200996(barcode: str) -> tuple[str, str]:
    match = re.search(r"_(P[^_]+)_([^_]+)$", barcode)
    return (match.group(1), normalize_timepoint(match.group(2))) if match else (UNKNOWN, UNKNOWN)


def cohort_from_obs(path_stem: str, obs: pd.DataFrame) -> str:
    cohort = value(obs, "cohort_id")
    known = cohort[cohort.ne(UNKNOWN)].unique()
    if len(known) == 1:
        return str(known[0])
    aliases = {
        "mendeley_skrx2fz79n": "mendeley_skrx2fz79n",
        "bi2021rcc": "bi_2021_rcc",
    }
    return aliases.get(path_stem.lower(), path_stem.upper())


def _base(obs: pd.DataFrame, cohort_id: str, source_object_id: str) -> pd.DataFrame:
    frame = obs.copy()
    frame.index = frame.index.astype(str)
    out = pd.DataFrame(index=frame.index)
    out["source_object_id"] = source_object_id
    out["source_cell_id"] = frame.index
    original = value(frame, "_source_barcode", default=UNKNOWN)
    out["raw_barcode"] = original.where(original.ne(UNKNOWN), pd.Series(frame.index, index=frame.index))
    out["cohort_id"] = cohort_id
    # Several legacy objects use donor_id/Patient for the source-proven
    # individual identifier. They are identity fields, never response labels.
    out["_patient"] = value(frame, "patient_id", "patient", "patient_alias", "donor_id", "Patient")
    out["_timepoint"] = value(frame, "timepoint", "timepoint_raw", "patient_timepoint_sort")
    out["_tissue"] = value(frame, "tissue_source", "tissue", "organ", "Tissue")
    out["_sample"] = value(frame, "sample", "sample_id", "sample_prefix", "orig.ident")
    out["_library"] = value(frame, "library_id", "library", "sample_prefix", "orig.ident", "seq_run")
    out["_demux"] = value(frame, "new_hash_ID", "hash")
    out["_treatment"] = value(frame, "treatment", "treatment_raw", "treatment_regimen_raw", "Cohort")
    out["_provenance"] = "obs_generic"
    out["source_sample_id"] = out["_sample"]
    out["patient_id"] = out["_patient"]
    out["treatment_arm"] = out["_treatment"]
    return out


def _gse123813(out: pd.DataFrame, lesion: str) -> pd.DataFrame:
    treatment = out["_timepoint"].map(normalize_timepoint)
    patient = out["_patient"]
    out["study_subject_key"] = [f"GSE123813::{p}" for p in patient]
    out["biological_sample_key"] = [key("specimen", "GSE123813", p, t, lesion) for p, t in zip(patient, treatment)]
    out["library_key"] = [key("library", "GSE123813", p, t, lesion) for p, t in zip(patient, treatment)]
    out["demux_id"] = UNKNOWN
    out["normalized_timepoint"] = treatment
    out["tissue_context"] = "tumor"
    out["lesion_context"] = lesion
    out["_provenance"] = "obs:patient_id+timepoint_raw+cohort_lesion"
    return out


def _gse232240(out: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    patient, timepoint, tissue = out["_patient"], out["_timepoint"].map(normalize_timepoint), out["_tissue"]
    seq_run, orig = value(raw, "seq_run"), value(raw, "orig.ident")
    out["library_key"] = [key("library", "GSE232240", a, b) for a, b in zip(seq_run, orig)]
    out["demux_id"] = value(raw, "new_hash_ID")
    out["study_subject_key"] = [f"GSE232240::{p}" for p in patient]
    out["biological_sample_key"] = [key("specimen", "GSE232240", p, t, d) for p, t, d in zip(patient, timepoint, out.demux_id)]
    out["normalized_timepoint"], out["tissue_context"], out["lesion_context"] = timepoint, tissue, "tumor"
    out["_provenance"] = "obs:seq_run+orig.ident;new_hash_ID;patient_id+timepoint"
    return out


def _gse176021(out: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    prefix, patient = value(raw, "sample_prefix"), out["_patient"]
    condition = pd.Series([
        p[len(subject) + 1:].rsplit("_", 1)[0] if p.startswith(subject + "_") else p
        for p, subject in zip(prefix, patient)
    ], index=out.index)
    lower = condition.str.lower()
    tissue = pd.Series("pbmc", index=out.index, dtype="object")
    tissue.loc[lower.eq("tumor")] = "tumor"
    tissue.loc[lower.eq("normal")] = "adjacent_normal"
    tissue.loc[lower.eq("ln")] = "tdln"
    tissue.loc[lower.eq("mettumor")] = "metastatic_tumor"
    timepoint = out["_timepoint"].map(normalize_timepoint)
    timepoint.loc[lower.eq("w2")] = "week_2"
    timepoint.loc[lower.eq("w4")] = "week_4"
    timepoint.loc[lower.eq("m3")] = "month_3"
    out["library_key"] = [key("library", "GSE176021", p) for p in prefix]
    out["demux_id"] = UNKNOWN
    out["study_subject_key"] = [f"GSE176021::{p}" for p in patient]
    out["biological_sample_key"] = [key("specimen", "GSE176021", p, s) for p, s in zip(patient, prefix)]
    out["normalized_timepoint"] = timepoint
    out["tissue_context"] = tissue
    out["lesion_context"] = condition
    out["_provenance"] = "obs:sample_prefix+patient_id;unknown_tokens_preserved"
    return out


def _gse272993(out: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    patient = value(raw, "patient_id", "patient_alias").map(canonical_272993)
    out["library_key"] = [key("library", "GSE272993", x, y) for x, y in zip(value(raw, "orig.ident"), value(raw, "hash"))]
    out["demux_id"] = value(raw, "hash")
    out["study_subject_key"] = [f"GSE272993::{p}" for p in patient]
    out["patient_id"] = patient
    out["treatment_arm"] = value(raw, "treatment")
    out["normalized_timepoint"] = out["_timepoint"].map(normalize_timepoint)
    out["tissue_context"] = out["_tissue"]
    out["lesion_context"] = UNKNOWN
    out["biological_sample_key"] = [
        key("specimen", "GSE272993", p, t, s, arm, source)
        for p, t, s, arm, source in zip(
            patient, out.normalized_timepoint, out.tissue_context,
            out.treatment_arm, out.source_sample_id)
    ]
    out["_provenance"] = "obs:patient_id+timepoint+tissue+treatment_arm+source_sample;14-1189_cross_arm_canonical"
    return out


def _gse272734(out: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    patient = value(raw, "patient_id")
    sample = value(raw, "sample_prefix")
    week = sample.str.extract(r"-wk(\d+)$", expand=False)
    timepoint = week.map(lambda x: "baseline" if str(x) == "0" else f"week_{int(x)}" if str(x).isdigit() else "unknown")
    out["study_subject_key"] = [f"GSE272734::{p}" for p in patient]
    out["patient_id"] = patient
    out["source_sample_id"] = sample
    out["library_key"] = [key("library", "GSE272734", s) for s in sample]
    out["demux_id"] = UNKNOWN
    out["biological_sample_key"] = [key("specimen", "GSE272734", p, s) for p, s in zip(patient, sample)]
    out["normalized_timepoint"] = timepoint
    out["tissue_context"] = "PBMC"
    out["lesion_context"] = "sorted_CD8_T"
    out["treatment_arm"] = "anti_PD1"
    out["_provenance"] = "obs:patient_id+sample_prefix_week;study_design:PBMC_sorted_CD8_anti_PD1"
    return out


def _gse200996(out: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    patient = value(raw, "Patient_ID", "patient_id", "donor_id")
    stage = value(raw, "Stage", "stage_raw", "timepoint").map(normalize_timepoint)
    tissue = value(raw, "tissue", "tissue_source", "organ")
    source = value(raw, "source_matrix", "orig.ident", "lib")
    out["study_subject_key"] = [f"GSE200996::{p}" for p in patient]
    out["patient_id"] = patient
    out["treatment_arm"] = value(raw, "Cohort")
    out["normalized_timepoint"] = stage
    out["tissue_context"] = tissue
    out["lesion_context"] = UNKNOWN
    out["library_key"] = [key("library", "GSE200996", item) for item in source]
    out["demux_id"] = UNKNOWN
    out["biological_sample_key"] = [key("specimen", "GSE200996", p, t, u, s) for p, t, u, s in zip(patient, stage, tissue, source)]
    out["_provenance"] = "obs:Patient_ID+Stage+tissue+source_matrix_or_lib"
    return out


def _generic(out: pd.DataFrame) -> pd.DataFrame:
    out["study_subject_key"] = [f"{c}::{p}" for c, p in zip(out.cohort_id, out._patient)]
    out["normalized_timepoint"] = out["_timepoint"].map(normalize_timepoint)
    out["tissue_context"] = out["_tissue"]
    out["lesion_context"] = UNKNOWN
    out["tissue_context"] = out["_tissue"]
    out["lesion_context"] = UNKNOWN
    out["library_key"] = [key("library", c, x) for c, x in zip(out.cohort_id, out._library)]
    out["demux_id"] = out["_demux"]
    out["biological_sample_key"] = [key("specimen", c, p, t, tissue, sample, library) for c, p, t, tissue, sample, library in zip(out.cohort_id, out._patient, out.normalized_timepoint, out.tissue_context, out._sample, out._library)]
    return out


def build_cell_identity(obs: pd.DataFrame, cohort_id: str, source_object_id: str) -> pd.DataFrame:
    out = _base(obs, cohort_id, source_object_id)
    raw = obs.copy(); raw.index = raw.index.astype(str)
    handlers: dict[str, Callable[[], pd.DataFrame]] = {
        "GSE123813_bcc": lambda: _gse123813(out, "BCC"),
        "GSE123813_scc": lambda: _gse123813(out, "SCC"),
        "GSE232240": lambda: _gse232240(out, raw),
        "GSE176021": lambda: _gse176021(out, raw),
        "GSE272993": lambda: _gse272993(out, raw),
        "GSE272734": lambda: _gse272734(out, raw),
        "GSE200996": lambda: _gse200996(out, raw),
    }
    out = handlers.get(cohort_id, lambda: _generic(out))()
    out["analysis_unit_key"] = [
        key("analysis", c, subject, t, tissue, lesion, arm)
        for c, subject, t, tissue, lesion, arm in zip(
            out.cohort_id, out.study_subject_key, out.normalized_timepoint,
            out.tissue_context, out.lesion_context, out.treatment_arm)
    ]
    missing = out.study_subject_key.str.endswith(f"::{UNKNOWN}") | out.library_key.eq(UNKNOWN) | out.biological_sample_key.eq(UNKNOWN)
    out["quarantine_reason"] = missing.map(lambda bad: "missing_required_identity_component" if bad else "")
    out["assignment_status"] = out.quarantine_reason.map(lambda reason: "quarantined" if reason else "resolved")
    out["assignment_provenance"] = out._provenance + ";metadata_only_no_response"
    return out.loc[:, CELL_COLUMNS].reset_index(drop=True)


def apply_conflict_quarantine(cells: pd.DataFrame) -> pd.DataFrame:
    duplicate = cells.duplicated(["source_object_id", "source_cell_id"], keep=False)
    # analysis_unit_key is a hash of precisely the four context fields below;
    # a distinct payload cannot share a key absent a hash collision.  Retain a
    # cheap explicit collision check on unique context rows, avoiding a
    # million-row groupby for large source objects.
    context = cells[["analysis_unit_key", "study_subject_key", "normalized_timepoint", "tissue_context", "lesion_context"]].drop_duplicates()
    conflict_keys = context.groupby("analysis_unit_key").size().loc[lambda s: s.gt(1)].index
    conflict = cells.analysis_unit_key.isin(conflict_keys)
    for mask, reason in ((duplicate, "duplicate_source_cell_key"), (conflict, "conflicting_analysis_unit_key")):
        cells.loc[mask, "quarantine_reason"] = reason
        cells.loc[mask, "assignment_status"] = "quarantined"
    return cells


def cardinality_audit(cells: pd.DataFrame) -> pd.DataFrame:
    return cells.groupby("cohort_id", dropna=False).agg(
        n_cells=("source_cell_id", "size"), n_libraries=("library_key", "nunique"),
        n_biological_samples=("biological_sample_key", "nunique"), n_subjects=("study_subject_key", "nunique"),
        n_analysis_units=("analysis_unit_key", "nunique"), resolved_cells=("assignment_status", lambda x: int(x.eq("resolved").sum())),
        quarantined_cells=("assignment_status", lambda x: int(x.eq("quarantined").sum())),
    ).reset_index().assign(identity_version=IDENTITY_VERSION)
