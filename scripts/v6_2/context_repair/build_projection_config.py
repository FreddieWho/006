#!/usr/bin/env python3
"""Build the immutable frozen-module projection config from repaired assets."""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
REGISTRY = ROOT / "results/v6_2/phase3_5_single_cell_processing_qc_gate/analysis_object_registry.frozen_v0.csv"
MEMBERSHIP = ROOT / "results/v6_2/phase6_module_algorithm_benchmark_strengthened/consensus/module_membership.frozen_v1.csv"
EXPRESSION_REFERENCE = ROOT / "results/v6_2/phase6_response_blind_module_discovery/foundation/logcpm_cellstate_pseudobulk.v0.parquet"


def _sample_values(node: h5py.Group | h5py.Dataset, limit: int = 10000) -> np.ndarray:
    if isinstance(node, h5py.Group):
        if "data" not in node:
            return np.array([])
        values = node["data"][:limit]
    else:
        if node.ndim == 1:
            values = node[:limit]
        else:
            rows = min(node.shape[0], 100)
            cols = min(node.shape[1], max(1, limit // max(rows, 1)))
            values = node[:rows, :cols]
    return np.asarray(values, dtype=float).ravel()


def infer_layer(path: Path, registered_decision: str | None) -> tuple[str | None, str, float | None]:
    """Admit only a source-proven or integer-audited count representation."""
    with h5py.File(path, "r") as handle:
        if "layers" in handle and "counts" in handle["layers"]:
            values = _sample_values(handle["layers/counts"])
            layer = "counts"
            source = "source_counts_layer"
        elif registered_decision == "raw_counts_in_X":
            values = _sample_values(handle["X"])
            layer = "X"
            source = "phase3_5_raw_counts_in_X"
        elif registered_decision in {"processed_or_unknown_X", "normalized_expression"}:
            return None, "registered_noncount_or_unknown_X", None
        else:
            values = _sample_values(handle["X"])
            layer = "X"
            source = "integer_audited_X"
    finite = values[np.isfinite(values)]
    noninteger = float(np.mean(np.abs(finite - np.round(finite)) > 1e-6)) if finite.size else 1.0
    if noninteger > 0.01 or (finite.size and np.nanmin(finite) < 0):
        return None, f"{source}_failed_integer_audit", noninteger
    return layer, source, noninteger


def infer_gene_column(path: Path) -> tuple[str | None, str]:
    """Select a source-proven gene symbol column when var_names are row numbers."""
    with h5py.File(path, "r") as handle:
        var = handle["var"]
        index_key = var.attrs.get("_index", "_index")
        index_key = index_key.decode() if isinstance(index_key, bytes) else str(index_key)
        index = var[index_key][: min(200, len(var[index_key]))]
        decoded = [x.decode() if isinstance(x, bytes) else str(x) for x in index]
        numeric_fraction = float(np.mean([value.isdigit() for value in decoded])) if decoded else 0.0
        ensembl_fraction = float(np.mean([value.upper().startswith("ENSG") for value in decoded])) if decoded else 0.0
        if (numeric_fraction > 0.9 or ensembl_fraction > 0.5) and "gene_symbol" in var:
            basis = "numeric" if numeric_fraction > 0.9 else "ensembl"
            return "gene_symbol", f"var_gene_symbol_replaces_{basis}_var_index"
        if ensembl_fraction > 0.5:
            return None, "ensembl_var_index_external_mapping"
    return None, "var_names"


def build(args: argparse.Namespace) -> dict:
    objects = pd.read_csv(args.interface_root / "object_metadata_audit.csv")
    registry = pd.read_csv(args.registry)
    by_path = {
        str(Path(row.object_path).resolve()): row
        for row in registry.itertuples()
        if getattr(row, "object_type", "") == "h5ad"
    }
    entries, audit = [], []
    for row in objects.itertuples():
        path = (ROOT / row.h5ad_path).resolve()
        registered = by_path.get(str(path))
        decision = str(getattr(registered, "layer_decision", "")) if registered else None
        layer, basis, noninteger = infer_layer(path, decision)
        gene_column, gene_basis = infer_gene_column(path)
        status = "included_raw_count_projection" if layer else "retained_support_not_count_projection"
        audit.append({
            "cohort_id": row.cohort_id,
            "object_id": row.source_object_id,
            "h5ad": str(path),
            "registered_layer_decision": decision or "not_in_frozen_v0_registry",
            "selected_count_layer": layer or "",
            "layer_evidence": basis,
            "selected_gene_column": gene_column or "var_names",
            "gene_namespace_evidence": gene_basis,
            "sample_noninteger_fraction": noninteger,
            "projection_status": status,
        })
        if not layer:
            continue
        entries.append({
            "cohort_id": str(row.cohort_id),
            "object_id": str(row.source_object_id),
            "h5ad": str(path),
            "count_layer": layer,
            "identity_sidecar": str(args.interface_root / "cell_identity.parquet" / f"cohort_id={row.cohort_id}"),
            **({"h5ad_gene_column": gene_column} if gene_column else {}),
            **(
                {"h5ad_gene_mapping": str(args.ensembl_map)}
                if not gene_column and gene_basis == "ensembl_var_index_external_mapping"
                else {}
            ),
            **(
                {"h5ad_cell_id_column": "__row_index__"}
                if str(row.cohort_id) in {"GSE272734", "GSE272735", "GSE273718"}
                else {}
            ),
        })
    if not entries:
        raise RuntimeError("No source-proven count objects are eligible for frozen projection")
    config = {
        "version": "context_repair_v1_frozen_projection",
        "inputs": {
            "objects": entries,
            "phase4a_annotations": str(args.phase4a_root / "handoff/cell_state_annotation_master.parquet"),
            "membership": str(args.membership),
            "expression_gene_reference": str(args.expression_reference),
            "annotation_object_column": "source_object_id",
            "annotation_cell_id_column": "source_cell_id",
            "identity_cell_id_column": "source_cell_id",
            "identity_status_column": "assignment_status",
            "analysis_unit_column": "analysis_unit_key",
            "coarse_state_column": "harmonized_coarse_label",
            "mid_state_column": "harmonized_mid_label",
            "chunk_size": args.chunk_size,
            "minimum_membership_coverage": args.minimum_membership_coverage,
        },
        "prohibitions": ["module_fit", "module_retraining", "response_guided_gene_selection", "membership_mutation"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True))
    pd.DataFrame(audit).to_csv(args.audit_output, index=False)
    return {"config": str(args.output), "n_included": len(entries), "n_support_only": len(audit) - len(entries)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interface-root", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1")
    parser.add_argument("--phase4a-root", type=Path, default=ROOT / "results/v6_2/phase4a_cell_state_harmonization_repair_v1")
    parser.add_argument("--registry", type=Path, default=REGISTRY)
    parser.add_argument("--membership", type=Path, default=MEMBERSHIP)
    parser.add_argument("--expression-reference", type=Path, default=EXPRESSION_REFERENCE)
    parser.add_argument("--ensembl-map", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1/gene_reference/ensembl_to_symbol_org_hs_eg_db.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1/configs/frozen_projection_v1.yaml")
    parser.add_argument("--audit-output", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1/configs/projection_object_eligibility_audit.csv")
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--minimum-membership-coverage", type=float, default=0.8)
    args = parser.parse_args()
    for name in ("interface_root", "phase4a_root", "registry", "membership", "expression_reference", "ensembl_map", "output", "audit_output"):
        setattr(args, name, getattr(args, name).resolve())
    print(yaml.safe_dump(build(args), sort_keys=False))


if __name__ == "__main__":
    main()
