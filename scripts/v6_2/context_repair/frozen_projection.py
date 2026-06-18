"""Streaming, immutable Phase6 module projection for repaired cell contexts.

The only operation here is ``count pseudobulk -> frozen membership projection``.
There is intentionally no fitting, feature selection, or membership mutation.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MEMBERSHIP = ROOT / "results/v6_2/phase6_module_algorithm_benchmark_strengthened/consensus/module_membership.frozen_v1.csv"
PHASE7_COLUMNS = ["expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "timepoint", "tissue_context", "cell_state_level", "cell_state", "n_cells_used", "library_size"]


class ProjectionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _path(value: str | Path, base: Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else (base / p).resolve()


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet" or path.is_dir():
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t" if ".tsv" in path.suffixes else ",")


def _col(frame: pd.DataFrame, requested: str | None, options: tuple[str, ...], label: str) -> str:
    if requested:
        if requested not in frame:
            raise ProjectionError(f"missing {label} column {requested!r}")
        return requested
    for name in options:
        if name in frame:
            return name
    raise ProjectionError(f"cannot infer {label}; tried {options}")


@dataclass(frozen=True)
class ObjectInput:
    cohort_id: str
    object_id: str
    h5ad: Path
    count_layer: str
    identity_sidecar: Path
    h5ad_cell_id_column: str | None
    h5ad_gene_column: str | None
    h5ad_gene_mapping: Path | None


@dataclass(frozen=True)
class Inputs:
    objects: list[ObjectInput]
    annotations: Path
    membership: Path
    annotation_object_column: str | None
    annotation_cell_id_column: str | None
    identity_cell_id_column: str | None
    identity_status_column: str | None
    analysis_unit_column: str
    coarse_state_column: str | None
    mid_state_column: str | None
    gene_column: str | None
    expression_gene_reference: Path | None
    chunk_size: int
    minimum_membership_coverage: float


def _objects(raw: dict[str, Any], base: Path) -> list[ObjectInput]:
    entries = raw.get("objects")
    if not entries:
        # Keep the first version's single-object config executable.
        entries = [raw]
    result = []
    for i, entry in enumerate(entries):
        missing = [k for k in ("cohort_id", "h5ad", "count_layer", "identity_sidecar") if not entry.get(k)]
        if missing:
            raise ProjectionError(f"objects[{i}] missing {', '.join(missing)}")
        item = ObjectInput(str(entry["cohort_id"]), str(entry.get("object_id", entry["cohort_id"])), _path(entry["h5ad"], base),
                           str(entry["count_layer"]), _path(entry["identity_sidecar"], base),
                           entry.get("h5ad_cell_id_column"), entry.get("h5ad_gene_column"),
                           _path(entry["h5ad_gene_mapping"], base) if entry.get("h5ad_gene_mapping") else None)
        if not item.h5ad.exists() or not item.identity_sidecar.exists():
            raise ProjectionError(f"object input missing: {item.h5ad} or {item.identity_sidecar}")
        if item.h5ad_gene_mapping and not item.h5ad_gene_mapping.exists():
            raise ProjectionError(f"gene mapping input missing: {item.h5ad_gene_mapping}")
        result.append(item)
    if len({x.object_id for x in result}) != len(result):
        raise ProjectionError("object_id must be unique")
    return result


def load_config(path: Path) -> tuple[dict[str, Any], Inputs]:
    config = yaml.safe_load(path.read_text()) or {}
    raw = config.get("inputs", config)
    if not raw.get("phase4a_annotations"):
        raise ProjectionError("missing inputs.phase4a_annotations")
    membership = _path(raw.get("membership", DEFAULT_MEMBERSHIP), path.parent)
    annotations = _path(raw["phase4a_annotations"], path.parent)
    if not membership.exists() or not annotations.exists():
        raise ProjectionError("membership or Phase4A annotation input does not exist")
    inputs = Inputs(
        objects=_objects(raw, path.parent), annotations=annotations, membership=membership,
        annotation_object_column=raw.get("annotation_object_column"), annotation_cell_id_column=raw.get("annotation_cell_id_column"),
        identity_cell_id_column=raw.get("identity_cell_id_column"), identity_status_column=raw.get("identity_status_column"),
        analysis_unit_column=str(raw.get("analysis_unit_column", "analysis_unit_key")),
        coarse_state_column=raw.get("coarse_state_column"), mid_state_column=raw.get("mid_state_column"),
        gene_column=raw.get("gene_column"),
        expression_gene_reference=_path(raw["expression_gene_reference"], path.parent)
        if raw.get("expression_gene_reference") else None,
        chunk_size=int(raw.get("chunk_size", 50_000)),
        minimum_membership_coverage=float(raw.get("minimum_membership_coverage", 0.8)),
    )
    if inputs.chunk_size < 1:
        raise ProjectionError("chunk_size must be positive")
    if inputs.expression_gene_reference and not inputs.expression_gene_reference.exists():
        raise ProjectionError("expression_gene_reference does not exist")
    return config, inputs


def _load_annotations(inputs: Inputs, obj: ObjectInput) -> pd.DataFrame:
    """Read one annotation partition; parquet directories are supported."""
    if inputs.annotations.is_dir():
        direct = inputs.annotations / f"cohort_id={obj.cohort_id}.parquet"
        if direct.exists():
            frame = pd.read_parquet(direct)
            object_col = inputs.annotation_object_column or "source_object_id"
            return frame.loc[frame[object_col].astype(str).eq(obj.object_id)].copy() if object_col in frame else frame
        try:
            import pyarrow.dataset as ds
            dataset = ds.dataset(inputs.annotations, format="parquet", partitioning="hive")
            object_col = inputs.annotation_object_column or "source_object_id"
            names = set(dataset.schema.names)
            if object_col in names:
                return dataset.to_table(filter=ds.field(object_col) == obj.object_id).to_pandas()
        except ImportError:
            pass
    frame = read_table(inputs.annotations)
    object_col = inputs.annotation_object_column or next((x for x in ("source_object_id", "object_id") if x in frame), None)
    if object_col:
        frame = frame.loc[frame[object_col].astype(str).eq(obj.object_id)].copy()
    return frame


def _membership(genes: pd.Index, path: Path, gene_column: str | None) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame]:
    frozen = pd.read_csv(path)
    module = _col(frozen, None, ("frozen_module_id", "module_id"), "membership module")
    gene = _col(frozen, gene_column, ("gene", "gene_symbol", "feature"), "membership gene")
    # relative_weight is scale-invariant across the original LDA representations.
    weight = "relative_weight" if "relative_weight" in frozen else _col(frozen, None, ("membership_weight", "weight", "consensus_weight"), "membership weight")
    x = frozen[[module, gene, weight]].copy()
    x[module], x[gene] = x[module].astype(str), x[gene].astype(str)
    x[weight] = pd.to_numeric(x[weight], errors="coerce")
    x = x.dropna(subset=[weight])
    totals = x.groupby(module, as_index=False)[weight].sum().rename(columns={weight: "total_weight"})
    index = pd.Series(np.arange(len(genes)), index=genes.astype(str)).groupby(level=0).first()
    x["gene_index"] = x[gene].map(index)
    overlap = x.dropna(subset=["gene_index"]).copy()
    if overlap.empty:
        raise ProjectionError("zero frozen-membership gene overlap with h5ad var_names")
    covered = overlap.groupby(module, as_index=False)[weight].sum().rename(columns={weight: "covered_weight"})
    coverage = totals.merge(covered, on=module, how="left").fillna({"covered_weight": 0})
    coverage["membership_coverage"] = coverage.covered_weight / coverage.total_weight.replace(0, np.nan)
    coverage = coverage.rename(columns={module: "module_id"})
    overlap = overlap.merge(covered, on=module, how="left")
    overlap["normalized_weight"] = overlap[weight] / overlap.covered_weight
    modules = sorted(coverage["module_id"].astype(str).unique())
    module_index = {name: i for i, name in enumerate(modules)}
    matrix = sparse.coo_matrix((overlap.normalized_weight, (overlap.gene_index.astype(int), overlap[module].map(module_index))),
                               shape=(len(genes), len(modules))).tocsr()
    return matrix, modules, coverage


def _object_gene_space(adata: Any, obj: ObjectInput) -> tuple[pd.Index, sparse.csr_matrix | None]:
    raw = pd.Index(
        adata.var[obj.h5ad_gene_column].astype(str)
        if obj.h5ad_gene_column
        else adata.var_names.astype(str)
    )
    if obj.h5ad_gene_mapping:
        mapping = pd.read_csv(obj.h5ad_gene_mapping)
        mapping = mapping.sort_values(["ENSEMBL", "mapping_rank"]).drop_duplicates("ENSEMBL")
        lookup = mapping.set_index("ENSEMBL")["SYMBOL"]
        mapped = pd.Series(raw.astype(str)).map(lookup)
        raw = pd.Index(mapped.where(mapped.notna(), pd.Series(raw.astype(str))).astype(str))
    if len(raw) != adata.n_vars:
        raise ProjectionError(f"{obj.object_id}: selected gene namespace is incomplete")
    genes = pd.Index(pd.unique(raw))
    if len(genes) == len(raw):
        return genes, None
    codes = genes.get_indexer(raw)
    collapse = sparse.coo_matrix(
        (np.ones(len(raw)), (np.arange(len(raw)), codes)),
        shape=(len(raw), len(genes)),
    ).tocsr()
    return genes, collapse


def _identity_table(obj: ObjectInput, inputs: Inputs, obs: pd.DataFrame, annotations: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    sidecar = read_table(obj.identity_sidecar)
    # A sidecar may be shared by several objects.  Its join identity is the
    # pair (source_object_id, source_cell_id), never a bare barcode globally.
    side_object = next((x for x in ("source_object_id", "object_id") if x in sidecar), None)
    if side_object:
        sidecar = sidecar.loc[sidecar[side_object].astype(str).eq(obj.object_id)].copy()
    side_id = _col(sidecar, inputs.identity_cell_id_column, ("source_cell_id", "original_cell_barcode", "cell_id", "cell.id", "barcode"), "identity cell id")
    ann_id = _col(annotations, inputs.annotation_cell_id_column, ("source_cell_id", "original_cell_barcode", side_id, "cell_id", "cell.id", "barcode"), "annotation cell id")
    unit = _col(sidecar, inputs.analysis_unit_column, ("analysis_unit_key",), "analysis-unit")
    status = _col(sidecar, inputs.identity_status_column, ("assignment_status", "identity_status"), "identity status")
    h5ad_id = obj.h5ad_cell_id_column
    source_cells = (
        pd.Series([f"row::{i}" for i in range(len(obs))], index=obs.index)
        if h5ad_id == "__row_index__"
        else obs[h5ad_id].astype(str) if h5ad_id
        else pd.Series(obs.index.astype(str), index=obs.index)
    )
    table = pd.DataFrame({"_row": np.arange(len(obs)), "_cell": source_cells.to_numpy()})
    for frame, key, name in ((sidecar, side_id, "sidecar"), (annotations, ann_id, "annotation")):
        if key not in frame or frame[key].duplicated().any():
            raise ProjectionError(f"{name} must contain one row per source cell")
    side_keep = [side_id, unit, status] + [x for x in (
        "study_subject_key", "patient_key", "patient_id", "normalized_timepoint", "timepoint",
        "tissue_context", "tissue_source", "lesion_context", "treatment_arm",
    ) if x in sidecar]
    table = table.merge(sidecar[side_keep].rename(columns={side_id: "_cell"}), on="_cell", how="left", validate="one_to_one")
    coarse = _col(annotations, inputs.coarse_state_column, ("cell_state_coarse", "coarse_label", "coarse_state"), "coarse cell-state")
    mid = _col(annotations, inputs.mid_state_column, ("cell_state_mid", "mid_label", "cell_state"), "mid cell-state")
    ann_extra = [x for x in ("analysis_unit_key", "patient_key", "study_subject_key", "normalized_timepoint",
                              "tissue_context", "assignment_status") if x in annotations]
    ann_keep = list(dict.fromkeys([ann_id, coarse, mid, *ann_extra]))
    ann = annotations[ann_keep].rename(columns={ann_id: "_cell", **{x: f"__ann_{x}" for x in ann_extra}})
    table = table.merge(ann, on="_cell", how="left", validate="one_to_one")
    ann_status = table.get("__ann_assignment_status", pd.Series("", index=table.index)).fillna("").astype(str)
    annotation_resolved = ann_status.str.lower().str.startswith("resolved")
    replacements = {
        unit: "analysis_unit_key", "study_subject_key": "study_subject_key", "patient_key": "patient_key",
        "normalized_timepoint": "normalized_timepoint", "tissue_context": "tissue_context",
    }
    for target, source_name in replacements.items():
        incoming = f"__ann_{source_name}"
        if incoming in table:
            table.loc[annotation_resolved, target] = table.loc[annotation_resolved, incoming]
    if "__ann_assignment_status" in table:
        table.loc[annotation_resolved, status] = table.loc[annotation_resolved, "__ann_assignment_status"]
    patient = table.get("patient_key", table.get("study_subject_key", table.get("patient_id", pd.Series(index=table.index, dtype=object))))
    if "study_subject_key" in table:
        patient = patient.where(patient.fillna("").astype(str).str.strip().ne(""), table["study_subject_key"])
    table["patient_key"] = patient.fillna("").astype(str)
    table["timepoint"] = table.get("normalized_timepoint", table.get("timepoint", pd.Series("unknown", index=table.index))).fillna("unknown").astype(str)
    table["tissue_context"] = table.get("tissue_context", table.get("tissue_source", pd.Series("unknown", index=table.index))).fillna("unknown").astype(str)
    table["_resolved"] = table[status].astype(str).str.lower().str.startswith(("resolved", "eligible"))
    counts = {"total_cells": len(table), "unresolved_identity": int((~table._resolved).sum()),
              "missing_analysis_unit": int(table[unit].isna().sum()), "missing_annotation": int(table[[coarse, mid]].isna().all(axis=1).sum())}
    table = table.rename(columns={unit: "analysis_unit_key", coarse: "_coarse", mid: "_mid"})
    return table, counts


def _project_object(obj: ObjectInput, inputs: Inputs, membership: sparse.csr_matrix, modules: list[str], coverage: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame]:
    import anndata as ad
    adata = ad.read_h5ad(obj.h5ad, backed="r")
    if obj.count_layer != "X" and obj.count_layer not in adata.layers:
        raise ProjectionError(f"{obj.object_id}: missing count layer {obj.count_layer!r}")
    annotations = _load_annotations(inputs, obj)
    identity, stats = _identity_table(obj, inputs, adata.obs, annotations)
    rows, expression_rows, conservation = [], [], []
    genes, collapse = _object_gene_space(adata, obj)
    for level, state_col in (("coarse", "_coarse"), ("mid", "_mid")):
        eligible = identity._resolved & identity.analysis_unit_key.notna() & identity[state_col].notna()
        eligible &= identity.analysis_unit_key.astype(str).str.strip().ne("") & identity[state_col].astype(str).str.strip().ne("")
        selected = identity.loc[eligible].copy()
        if selected.empty:
            continue
        key_cols = ["analysis_unit_key", "patient_key", "timepoint", "tissue_context", state_col]
        keys = selected[key_cols].astype(str).drop_duplicates().sort_values(key_cols).reset_index(drop=True)
        key_map = {tuple(row): i for i, row in keys.iterrows()}
        codes = np.fromiter((key_map[tuple(row)] for _, row in selected[key_cols].astype(str).iterrows()), int, len(selected))
        total = sparse.csr_matrix((len(keys), len(genes)), dtype=np.float64)
        eligible_total = 0.0
        selected["_code"] = codes
        for start in range(0, len(selected), inputs.chunk_size):
            chunk = selected.iloc[start:start + inputs.chunk_size]
            source = adata.X if obj.count_layer == "X" else adata.layers[obj.count_layer]
            matrix = source[chunk._row.to_numpy()]
            matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
            # AnnData count layers are often float32 even when values are
            # integer counts. Summing large cohorts in float32 can lose a few
            # units and trigger a false conservation failure.
            matrix = matrix.astype(np.float64, copy=False)
            if collapse is not None:
                matrix = matrix @ collapse
            eligible_total += float(matrix.sum())
            selector = sparse.csr_matrix((np.ones(len(chunk)), (chunk._code.to_numpy(), np.arange(len(chunk)))), shape=(len(keys), len(chunk)))
            total += selector @ matrix
        score = total @ membership
        registry = pd.DataFrame({"analysis_unit_key": keys.analysis_unit_key, "patient_key": keys.patient_key, "timepoint": keys.timepoint,
                                 "tissue_context": keys.tissue_context, "cell_state": keys[state_col], "cell_state_level": level})
        registry["cohort_id"], registry["object_id"], registry["sample_key"] = obj.cohort_id, obj.object_id, registry.analysis_unit_key
        registry["expression_unit_id"] = [f"{obj.cohort_id}::{obj.object_id}::{u}::{level}::{s}::count_pseudobulk" for u, s in zip(registry.analysis_unit_key, registry.cell_state)]
        registry["n_cells_used"] = np.bincount(codes, minlength=len(keys))
        registry["library_size"] = np.asarray(total.sum(axis=1)).ravel()
        registry["n_genes"] = np.asarray((total > 0).sum(axis=1)).ravel()
        registry["inclusion_status"] = np.where(registry.library_size > 0, "included", "zero_library")
        registry["layer_family"] = "raw_count"
        score_frame = registry[PHASE7_COLUMNS].copy()
        score_frame[modules] = score.toarray() if sparse.issparse(score) else np.asarray(score)
        low = set(coverage.loc[coverage.membership_coverage < inputs.minimum_membership_coverage, "module_id"])
        for module in low:
            if module in score_frame:
                score_frame[module] = np.nan
        # This is deliberately dense only *after* streaming aggregation: its
        # shape is expression-units x genes, never cells x genes.
        library = np.asarray(total.sum(axis=1)).ravel()
        scale = np.divide(1_000_000.0, library, out=np.zeros_like(library, dtype=float), where=library > 0)
        expression_indices = np.arange(len(genes))
        expression_genes = genes
        if inputs.expression_gene_reference:
            import pyarrow.parquet as pq
            reference_columns = pq.ParquetFile(inputs.expression_gene_reference).schema.names
            identifier_columns = set(PHASE7_COLUMNS) | {"n_genes", "inclusion_status", "layer_family"}
            requested = [x for x in reference_columns if x not in identifier_columns]
            position = pd.Series(np.arange(len(genes)), index=genes.astype(str)).groupby(level=0).first()
            available = [gene for gene in requested if gene in position.index]
            expression_indices = position.reindex(available).to_numpy(dtype=int)
            expression_genes = pd.Index(available)
        logcpm = total[:, expression_indices].multiply(scale[:, None]).tocsr()
        logcpm.data = np.log1p(logcpm.data)
        expression = pd.concat(
            [
                registry[PHASE7_COLUMNS].reset_index(drop=True),
                pd.DataFrame(logcpm.toarray(), columns=expression_genes),
            ],
            axis=1,
        )
        rows.append((registry, score_frame))
        expression_rows.append(expression)
        conservation.append({"cohort_id": obj.cohort_id, "object_id": obj.object_id, "cell_state_level": level,
                             "eligible_cell_total": eligible_total,
                             "pseudobulk_total": float(total.sum()), "n_cells_used": len(selected)})
    if not rows:
        raise ProjectionError(f"{obj.object_id}: no resolved cells with Phase4A coarse or mid annotation")
    registry, scores = (pd.concat([x[i] for x in rows], ignore_index=True) for i in (0, 1))
    expression = pd.concat(expression_rows, ignore_index=True)
    conservation = pd.DataFrame(conservation)
    conservation["eligible_equals_pseudobulk"] = np.isclose(conservation.eligible_cell_total, conservation.pseudobulk_total, rtol=0, atol=0)
    conservation_rows = conservation.assign(
        check="counts_conservation",
        status=np.where(conservation.eligible_equals_pseudobulk, "PASS", "HARD_FAIL"),
        evidence=conservation.apply(lambda row: json.dumps(row.to_dict(), default=str), axis=1),
    )
    blocks = pd.DataFrame([
        {"cohort_id": obj.cohort_id, "object_id": obj.object_id, "check": "identity_status", "status": "BLOCK" if stats["unresolved_identity"] else "PASS", "evidence": json.dumps(stats)},
    ] + [{"cohort_id": obj.cohort_id, "object_id": obj.object_id, "check": f"membership_coverage::{r.module_id}",
           "status": "LOW_COVERAGE" if r.membership_coverage < inputs.minimum_membership_coverage else "PASS", "evidence": str(r.membership_coverage)} for r in coverage.itertuples()])
    return registry, scores, expression, {**stats, "cohort_id": obj.cohort_id, "object_id": obj.object_id}, pd.concat([conservation_rows, blocks], ignore_index=True, sort=False)


def run_projection(config: dict[str, Any], inputs: Inputs, output_root: Path, cohort: str | None = None) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=True)
    before = sha256(inputs.membership)
    registries, scores, expressions, coverage_rows, audits = [], [], [], [], []
    selected = [x for x in inputs.objects if not cohort or x.cohort_id == cohort]
    if not selected:
        raise ProjectionError(f"--cohort {cohort!r} matches no configured object")
    # Membership is made per object because h5ad feature namespaces can differ.
    for obj in selected:
        import anndata as ad
        probe = ad.read_h5ad(obj.h5ad, backed="r")
        genes, _ = _object_gene_space(probe, obj)
        membership, modules, gene_coverage = _membership(genes, inputs.membership, inputs.gene_column)
        reg, score, expression, stats, audit = _project_object(obj, inputs, membership, modules, gene_coverage)
        registries.append(reg); scores.append(score); expressions.append(expression); coverage_rows.append(gene_coverage.assign(cohort_id=obj.cohort_id, object_id=obj.object_id)); audits.append(audit)
    after = sha256(inputs.membership)
    if before != after:
        raise ProjectionError("frozen membership SHA256 changed during projection")
    registry, score = pd.concat(registries, ignore_index=True), pd.concat(scores, ignore_index=True)
    expression = pd.concat(expressions, ignore_index=True, sort=False)
    coverage = pd.concat(coverage_rows, ignore_index=True)
    audit = pd.concat(audits, ignore_index=True, sort=False)
    audit.to_csv(output_root / "counts_conservation_and_block_audit.csv", index=False)
    blocking = audit.get("status", pd.Series(dtype=str)).isin({"HARD_FAIL", "BLOCK"})
    if blocking.any():
        blocked_manifest = {
            "verdict": "BLOCKED", "membership_sha256_before": before,
            "membership_sha256_after": after, "membership_path": str(inputs.membership),
            "blocking_rows": int(blocking.sum()), "formal_outputs_published": False,
        }
        (output_root / "projection_manifest.yaml").write_text(yaml.safe_dump(blocked_manifest, sort_keys=False))
        raise ProjectionError("identity or count-conservation block; formal projection not published")
    registry.to_csv(output_root / "unified_expression_unit_registry.csv", index=False)
    score.to_parquet(output_root / "frozen_module_score_matrix.parquet", index=False)
    score.to_csv(output_root / "frozen_module_score_matrix.csv", index=False)
    expression.to_parquet(output_root / "repaired_cell_state_pseudobulk.logcpm.parquet", index=False)
    coverage.to_csv(output_root / "gene_overlap_coverage_audit.csv", index=False)
    manifest = {"verdict": "PASS", "membership_sha256_before": before, "membership_sha256_after": after, "membership_path": str(inputs.membership),
                "prohibitions": ["module_fit", "module_retraining", "gene_selection", "membership_mutation"],
                "objects": [x.object_id for x in selected], "n_expression_units": len(registry), "n_score_rows": len(score),
                "phase7_logcpm_parquet": str(output_root / "repaired_cell_state_pseudobulk.logcpm.parquet")}
    (output_root / "projection_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    return manifest
