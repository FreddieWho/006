#!/usr/bin/env python3
"""Run v7 Stage 2 response-blind ontology and measurement."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from src.v7.ontology.aggregate import aggregate_object, load_object_specs, load_shards, shard_input_fingerprint
from src.v7.ontology.contracts import (
    Stage2Error,
    architecture_surrogate_contract,
    content_manifest,
    load_config,
    prohibited_columns,
    require_inputs,
    sha256,
    write_json,
    write_tsv,
    write_yaml,
)
from src.v7.ontology.diagnostics import (
    ambient_diagnostics,
    carrier_diagnostics,
    confounding_diagnostics,
    cross_modal_mapping,
    gse193736_repeatability_diagnostics,
    gse207422_pair_vector_diagnostics,
    gse207422_paired_descriptive_diagnostics,
    leave_family_out_state_profiles,
    matched_null_coherence,
    method_agreement_diagnostics,
    patient_timepoint_scores,
    robust_standardize_scores,
    summarize_reliability,
    technical_panel_reliability,
)
from src.v7.ontology.pipeline import (
    build_gene_resolver,
    compatibility_audit,
    compile_vocabulary,
    family_state_gene_profiles,
    load_written_vocabulary,
    score_bridges,
    score_gse207422_whole_sample,
    score_pseudobulk_shards,
    write_vocabulary,
)
from src.v7.ontology.technical_panel import run_cell_level_technical_panel


def _preflight(config: dict[str, Any], output_root: Path) -> None:
    diagnostics = output_root / "diagnostics"
    diagnostics.mkdir(parents=True, exist_ok=True)
    inputs = require_inputs(config)
    write_tsv(inputs, diagnostics / "input_audit.tsv", ["input_role"])
    specs = load_object_specs(config["inputs"]["v6_projection_config"])
    expected = int(config["execution"]["required_count_compatible_objects"])
    if len(specs) != expected:
        raise Stage2Error(f"STAGE2_BLOCKED_OBJECT_COUNT: expected={expected}, observed={len(specs)}")
    leakage_rows = []
    for role in ("stage1_registry", "stage1_roles"):
        path = config["inputs"][role]
        columns = pd.read_csv(path, sep="\t", nrows=0).columns
        hits = prohibited_columns(columns, config["prohibited_tokens"])
        leakage_rows.append(
            {
                "source": role,
                "prohibited_columns_present_in_source": "|".join(hits),
                "n_prohibited_source_columns": len(hits),
                "read_into_measurement": False,
                "status": "PASS_DROPPED_AT_SOURCE" if hits else "PASS",
            }
        )
    leakage_rows.extend(
        [
            {
                "source": "GSE207422_metadata",
                "prohibited_columns_present_in_source": "not_loaded_usecols_Sample_Patient_Resource_only",
                "n_prohibited_source_columns": "not_read",
                "read_into_measurement": False,
                "status": "PASS_COLUMN_WHITELIST",
            },
            {
                "source": "score_outputs",
                "prohibited_columns_present_in_source": "",
                "n_prohibited_source_columns": 0,
                "read_into_measurement": False,
                "status": "PASS_PENDING_OUTPUT_AUDIT",
            },
        ]
    )
    write_tsv(pd.DataFrame(leakage_rows), diagnostics / "leakage_preflight.tsv", ["source"])


def _vocabulary(config: dict[str, Any], output_root: Path) -> None:
    vocabulary = compile_vocabulary(config)
    write_vocabulary(vocabulary, output_root)


def _aggregation_semantic_context(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "algorithm_contract": "v7_count_pseudobulk_identity_bound_v2",
        "stage_version": str(config["version"]),
        "primary_levels": list(config["execution"]["primary_levels"]),
        "gene_mapping_sha256": {
            role: sha256(config["inputs"][role])
            for role in ("gencode_gene_map", "org_hs_gene_map")
        },
    }


def _previous_hash_manifest(output_root: Path) -> dict[str, dict[str, Any]]:
    path = output_root / "diagnostics" / "raw_object_input_hashes.json"
    if not path.exists():
        return {}
    try:
        records = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return {
        str(record["object_id"]): record
        for record in records
        if isinstance(record, dict) and record.get("object_id")
    }


def _migrate_legacy_shard_fingerprints(
    spec: Any,
    output_root: Path,
    expected_fingerprint: dict[str, Any],
    previous: dict[str, Any],
) -> None:
    """Refresh an old digest when the shard already carries source hashes.

    An audit with no input fingerprint is intentionally not migrated: its
    source-to-shard binding is unknowable, so ``_shard_ready`` must force a
    fresh aggregation rather than relying on a later hash manifest.
    """

    safe = spec.object_id.replace("::", "__").replace("/", "_")
    for level in ("coarse", "mid"):
        audit_path = output_root / "cache" / "pseudobulk" / f"{safe}__{level}.audit.json"
        if not audit_path.exists():
            return
        record = json.loads(audit_path.read_text(encoding="utf-8"))
        if (
            record.get("h5ad_bytes") != spec.h5ad.stat().st_size
            or not record.get("counts_conserved")
        ):
            return
        observed = record.get("input_fingerprint") or {}
        if not observed:
            return
        # The fingerprint digest deliberately excludes provenance labels.
        # Refresh an old audit when every source/semantic field already
        # matches, while refusing to bless a genuinely changed input.
        comparable = {
            key: value
            for key, value in expected_fingerprint.items()
            if key not in {"fingerprint_sha256", "hash_provenance"}
        }
        observed_comparable = {
            key: value
            for key, value in observed.items()
            if key not in {"fingerprint_sha256", "hash_provenance"}
        }
        if observed_comparable != comparable:
            return
        if observed.get("fingerprint_sha256") == expected_fingerprint.get("fingerprint_sha256"):
            continue
        record["input_fingerprint"] = expected_fingerprint
        write_json(record, audit_path)


def _shard_ready(
    spec: Any,
    output_root: Path,
    expected_fingerprint: dict[str, Any],
) -> bool:
    safe = spec.object_id.replace("::", "__").replace("/", "_")
    for level in ("coarse", "mid"):
        base = output_root / "cache" / "pseudobulk" / f"{safe}__{level}"
        audit = base.with_suffix(".audit.json")
        required = [
            base.with_suffix(".counts.npz"),
            base.with_suffix(".units.parquet"),
            base.with_suffix(".genes.json"),
            audit,
        ]
        if not all(path.exists() for path in required):
            return False
        record = json.loads(audit.read_text(encoding="utf-8"))
        observed_fingerprint = record.get("input_fingerprint") or {}
        if (
            record.get("h5ad_bytes") != spec.h5ad.stat().st_size
            or not record.get("counts_conserved")
            or observed_fingerprint.get("fingerprint_sha256")
            != expected_fingerprint.get("fingerprint_sha256")
        ):
            return False
    return True


def _aggregate(
    config: dict[str, Any],
    output_root: Path,
    *,
    resume: bool,
    cohort: str | None,
    workers: int | None,
) -> None:
    specs = load_object_specs(config["inputs"]["v6_projection_config"])
    if cohort:
        specs = [spec for spec in specs if spec.cohort_id == cohort]
        if not specs:
            raise Stage2Error(f"STAGE2_BLOCKED_UNKNOWN_COHORT: {cohort}")
    overrides = config.get("study_family_overrides") or {}
    semantic_context = _aggregation_semantic_context(config)
    previous_manifest = _previous_hash_manifest(output_root)
    fingerprints = {
        spec.object_id: shard_input_fingerprint(
            spec,
            study_family=overrides.get(spec.cohort_id, spec.cohort_id),
            semantic_context=semantic_context,
            previous=previous_manifest.get(spec.object_id),
        )
        for spec in specs
    }
    if resume:
        for spec in specs:
            _migrate_legacy_shard_fingerprints(
                spec,
                output_root,
                fingerprints[spec.object_id],
                previous_manifest.get(spec.object_id, {}),
            )
    pending = [
        spec
        for spec in specs
        if not (
            resume
            and _shard_ready(spec, output_root, fingerprints[spec.object_id])
        )
    ]
    chunk_size = int(config["execution"]["chunk_size"])
    max_workers = int(workers or config["execution"].get("workers", 1))
    if max_workers < 1:
        raise Stage2Error("STAGE2_BLOCKED_WORKERS: workers must be positive")
    if max_workers == 1:
        for spec in pending:
            aggregate_object(
                spec,
                output_root,
                chunk_size,
                overrides.get(spec.cohort_id, spec.cohort_id),
                fingerprints[spec.object_id],
            )
            print(json.dumps({"event": "aggregate_complete", "cohort_id": spec.cohort_id}), flush=True)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    aggregate_object,
                    spec,
                    output_root,
                    chunk_size,
                    overrides.get(spec.cohort_id, spec.cohort_id),
                    fingerprints[spec.object_id],
                ): spec
                for spec in pending
            }
            for future in as_completed(futures):
                spec = futures[future]
                future.result()
                print(json.dumps({"event": "aggregate_complete", "cohort_id": spec.cohort_id}), flush=True)


def _score(config: dict[str, Any], output_root: Path) -> None:
    vocabulary = load_written_vocabulary(output_root)
    resolver = build_gene_resolver(config)
    sc_scores, sc_mapping = score_pseudobulk_shards(
        output_root,
        vocabulary,
        resolver,
        workers=int(config["execution"].get("workers", 1)),
    )
    bridge_scores, bridge_mapping, crosswalk, design = score_bridges(config, vocabulary, resolver)
    scores = robust_standardize_scores(pd.concat([sc_scores, bridge_scores], ignore_index=True))
    score_root = output_root / "scores"
    score_root.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(score_root / "expression_unit_scores.parquet", index=False)
    patients = patient_timepoint_scores(sc_scores.loc[sc_scores.modality.eq("scRNA")])
    patients = robust_standardize_scores(
        patients,
        group_columns=("cohort_id", "modality", "aggregation_level", "aggregation_method", "cell_state_level", "feature_id"),
    )
    patients.to_parquet(score_root / "patient_timepoint_scores.parquet", index=False)
    write_tsv(crosswalk, output_root / "gse207422_patient_timepoint_crosswalk.tsv", ["patient_id"])
    write_tsv(design, output_root / "gse193736_perturbation_design.tsv", ["lineage", "perturbation", "rest_stim", "replicate"])
    mapping = pd.concat([sc_mapping, bridge_mapping], ignore_index=True, sort=False)
    mapping = (
        mapping.groupby(
            ["input_id", "normalized_id", "resolved_symbol", "mapping_source", "mapping_status", "candidates"],
            as_index=False,
            dropna=False,
        )
        .agg(n_objects=("object_id", "nunique"), n_cohorts=("cohort_id", "nunique"))
    )
    write_tsv(mapping, output_root / "gene_mapping_audit.tsv", ["mapping_status", "input_id"])
    compatibility = compatibility_audit(sc_scores, config["inputs"]["legacy_score_reference"])
    write_tsv(compatibility, output_root / "diagnostics" / "legacy_fm_compatibility.tsv", ["feature_id"])


def _technical_panel(config: dict[str, Any], output_root: Path) -> None:
    vocabulary = load_written_vocabulary(output_root)
    resolver = build_gene_resolver(config)
    scores = pd.read_parquet(output_root / "scores" / "expression_unit_scores.parquet")
    specs = load_object_specs(config["inputs"]["v6_projection_config"])
    summary, comparison = run_cell_level_technical_panel(
        specs,
        config["execution"]["cell_level_technical_panel"],
        vocabulary,
        resolver,
        scores,
        chunk_size=int(config["execution"]["chunk_size"]),
    )
    write_tsv(
        summary,
        output_root / "diagnostics" / "cell_level_technical_panel.tsv",
        ["cohort_id", "object_id", "cell_state_level", "feature_family", "feature_id"],
    )
    cache_path = output_root / "cache" / "cell_level_technical_panel_units.parquet"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_parquet(cache_path, index=False)


def _parent_reliability(reliability: pd.DataFrame, vocabulary: Any) -> pd.DataFrame:
    rows = []
    parents = vocabulary.module_dictionary.loc[vocabulary.module_dictionary.feature_level.eq("parent_axis")]
    for parent in parents.itertuples(index=False):
        components = vocabulary.module_dictionary.loc[
            vocabulary.module_dictionary.parent_axis_id.eq(parent.feature_id), "feature_id"
        ]
        evidence = reliability.loc[reliability.feature_id.isin(components)]
        if evidence.empty or evidence.technical_status.eq("valid").sum() == 0:
            technical, d2, reason = "reject", "annotation_only_or_reject", "no_scoreable_components"
        elif evidence.D2_class.eq("measurable").all():
            technical, d2, reason = "valid", "measurable", "all_explicit_components_measurable"
        else:
            technical, d2, reason = "valid", "joint_representation_with_uncertainty", "vector_components_have_mixed_measurement_evidence"
        rows.append(
            {
                "feature_id": parent.feature_id,
                "feature_family": "mechanism_axis",
                "n_expression_units": int(evidence.n_expression_units.max()) if not evidence.empty else 0,
                "n_study_families": int(evidence.n_study_families.max()) if not evidence.empty else 0,
                "median_coverage": float(evidence.median_coverage.median()) if not evidence.empty else np.nan,
                "technical_status": technical,
                "D2_class": d2,
                "D2_reason": reason,
                "n_components": len(components),
                "n_measurable_components": int(evidence.D2_class.eq("measurable").sum()),
            }
        )
    return pd.DataFrame(rows)


def _gse207422_paired_score_table(
    scores: pd.DataFrame,
    whole_sample_sc: pd.DataFrame,
    crosswalk: pd.DataFrame,
) -> pd.DataFrame:
    """Bind whole-sample scRNA and bulk scores using the audited crosswalk."""

    def patient_id(values: pd.Series) -> pd.Series:
        return values.astype(str).str.rsplit("::", n=1).str[-1]

    def normalized_timepoint(value: Any) -> str:
        text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        if text in {"baseline", "pre", "pre_treatment", "pretreatment"}:
            return "pre_treatment"
        if text in {"post", "post_treatment", "posttreatment"}:
            return "post_treatment"
        return text

    sc = whole_sample_sc.copy()
    sc["patient_id"] = patient_id(sc.patient_key)
    bulk = scores.loc[
        scores.cohort_id.eq("GSE207422_bulk") & scores.modality.eq("bulkRNA")
    ].copy()
    bulk["patient_id"] = patient_id(bulk.patient_key)
    rows: list[pd.DataFrame] = []
    for bridge in crosswalk.itertuples(index=False):
        sc_pair = sc.loc[sc.patient_id.eq(str(bridge.patient_id))].copy()
        bulk_pair = bulk.loc[bulk.patient_id.eq(str(bridge.patient_id))].copy()
        if sc_pair.empty or bulk_pair.empty:
            raise Stage2Error(
                f"STAGE2_BLOCKED_GSE207422_PAIR_SCORE_MISSING: {bridge.patient_id}"
            )
        if (
            sc_pair.groupby("feature_id").size().gt(1).any()
            or bulk_pair.groupby("feature_id").size().gt(1).any()
        ):
            raise Stage2Error(
                f"STAGE2_BLOCKED_GSE207422_PAIR_SCORE_DUPLICATE: {bridge.patient_id}"
            )
        observed_sc_time = {
            normalized_timepoint(value) for value in sc_pair.timepoint.unique()
        }
        observed_bulk_time = {
            normalized_timepoint(value) for value in bulk_pair.timepoint.unique()
        }
        if observed_sc_time != {normalized_timepoint(bridge.sc_timepoint)}:
            raise Stage2Error(
                f"STAGE2_BLOCKED_GSE207422_SC_TIMEPOINT: {bridge.patient_id}"
            )
        if observed_bulk_time != {normalized_timepoint(bridge.bulk_timepoint)}:
            raise Stage2Error(
                f"STAGE2_BLOCKED_GSE207422_BULK_TIMEPOINT: {bridge.patient_id}"
            )
        for modality, frame in (("scRNA", sc_pair), ("bulkRNA", bulk_pair)):
            material = frame[
                ["feature_id", "score_standardized", "coverage", "technical_status"]
            ].copy()
            material["patient_id"] = str(bridge.patient_id)
            material["pair_id"] = str(bridge.patient_id)
            material["modality"] = modality
            material["match_status"] = str(bridge.match_status)
            material["bridge_role"] = str(bridge.bridge_role)
            material["sc_timepoint"] = str(bridge.sc_timepoint)
            material["bulk_timepoint"] = str(bridge.bulk_timepoint)
            rows.append(material)
    result = pd.concat(rows, ignore_index=True)
    if result.duplicated(["pair_id", "feature_id", "modality"]).any():
        raise Stage2Error("STAGE2_BLOCKED_GSE207422_PAIRED_TABLE_DUPLICATE")
    return result


def _diagnose(config: dict[str, Any], output_root: Path) -> None:
    vocabulary = load_written_vocabulary(output_root)
    scores = pd.read_parquet(output_root / "scores" / "expression_unit_scores.parquet")
    patients = pd.read_parquet(output_root / "scores" / "patient_timepoint_scores.parquet")
    sc = scores.loc[scores.modality.eq("scRNA")].copy()
    lfo = leave_family_out_state_profiles(sc)
    carrier = carrier_diagnostics(patients)
    ambient = ambient_diagnostics(sc)
    confounding = confounding_diagnostics(sc)
    method = method_agreement_diagnostics(
        sc,
        seed=int(config["seed"]),
        bootstrap_replicates=int(config["execution"]["bootstrap_replicates"]),
    )
    technical_panel = technical_panel_reliability(
        pd.read_csv(
            output_root / "diagnostics" / "cell_level_technical_panel.tsv",
            sep="\t",
        ),
        seed=int(config["seed"]),
        bootstrap_replicates=int(config["execution"]["bootstrap_replicates"]),
    )
    resolver = build_gene_resolver(config)
    profile_meta, profile_matrix, profile_symbols, profile_observed = family_state_gene_profiles(output_root, resolver)
    null_summary, null_records = matched_null_coherence(
        profile_matrix,
        profile_symbols,
        vocabulary.module_membership,
        vocabulary.module_dictionary,
        replicates=int(config["execution"]["null_replicates"]),
        seed=int(config["seed"]),
    )
    reliability = summarize_reliability(
        sc,
        lfo,
        carrier,
        ambient,
        confounding,
        method,
        null_intervals=null_summary,
        technical_panel=technical_panel,
        seed=int(config["seed"]),
        bootstrap_replicates=int(config["execution"]["bootstrap_replicates"]),
    )
    reliability = pd.concat([reliability, _parent_reliability(reliability, vocabulary)], ignore_index=True, sort=False)
    reliability["evidence_modality"] = "scRNA"
    reliability["D2_scope"] = "scRNA_measurement_only"
    write_tsv(reliability, output_root / "measurement_reliability.tsv", ["feature_family", "feature_id"])
    whole_sample_sc = score_gse207422_whole_sample(
        output_root, vocabulary, resolver
    )
    crosswalk = pd.read_csv(
        output_root / "gse207422_patient_timepoint_crosswalk.tsv", sep="\t"
    )
    paired_scores = _gse207422_paired_score_table(
        scores, whole_sample_sc, crosswalk
    )
    paired_bulk = gse207422_paired_descriptive_diagnostics(
        paired_scores, value_column="score_standardized"
    )
    paired_vectors = gse207422_pair_vector_diagnostics(
        paired_scores, value_column="score_standardized"
    )
    perturbation_design = pd.read_csv(
        output_root / "gse193736_perturbation_design.tsv", sep="\t"
    )
    perturb_repeatability = gse193736_repeatability_diagnostics(
        scores.loc[scores.modality.eq("bulk_perturbation")],
        perturbation_design,
    )
    mapping = cross_modal_mapping(
        scores,
        reliability,
        paired_bulk=paired_bulk,
        perturb_repeatability=perturb_repeatability,
    )
    spatial_rows = []
    for feature in vocabulary.module_dictionary.itertuples(index=False):
        for level in config["spatial_contract"]["levels"]:
            spatial_rows.append(
                {
                    "feature_family": feature.feature_family,
                    "feature_id": feature.feature_id,
                    "modality": "spatial",
                    "aggregation_level": level,
                    "n_observations": 0,
                    "n_scoreable": 0,
                    "n_datasets": 0,
                    "n_patients": 0,
                    "median_coverage": np.nan,
                    "technical_status": "not_run",
                    "D2_class": "pending_stage3",
                    "D2_reason": "NOT_RUN_PENDING_STAGE3",
                }
            )
    mapping = pd.concat([mapping, pd.DataFrame(spatial_rows)], ignore_index=True, sort=False)
    write_tsv(mapping, output_root / "cross_modal_mapping.tsv", ["feature_family", "feature_id", "modality", "aggregation_level"])
    diagnostics = output_root / "diagnostics"
    write_tsv(lfo, diagnostics / "leave_study_family_out.tsv", ["feature_id", "heldout_study_family"])
    write_tsv(carrier, diagnostics / "carrier_dependency.tsv", ["feature_id"])
    write_tsv(ambient, diagnostics / "ambient_observed_profile.tsv", ["feature_id"])
    write_tsv(confounding, diagnostics / "confounding_diagnostics.tsv", ["feature_id"])
    write_tsv(method, diagnostics / "method_agreement.tsv", ["feature_id"])
    write_tsv(null_summary, diagnostics / "matched_expression_null_summary.tsv", ["feature_id"])
    write_tsv(null_records, diagnostics / "matched_expression_null_replicates.tsv", ["feature_id", "null_id"])
    write_tsv(
        paired_scores,
        diagnostics / "gse207422_paired_scores.tsv",
        ["pair_id", "feature_id", "modality"],
    )
    write_tsv(
        paired_bulk,
        diagnostics / "gse207422_paired_feature_diagnostics.tsv",
        ["feature_id"],
    )
    write_tsv(
        paired_vectors,
        diagnostics / "gse207422_pair_vector_diagnostics.tsv",
        ["comparison_type", "comparison_id"],
    )
    write_tsv(
        perturb_repeatability,
        diagnostics / "gse193736_repeatability.tsv",
        ["feature_id"],
    )
    audit_rows = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted((output_root / "cache" / "pseudobulk").glob("*.audit.json"))
    ]
    write_tsv(
        pd.DataFrame(audit_rows),
        diagnostics / "aggregation_audit.tsv",
        ["cohort_id", "object_id", "cell_state_level"],
    )
    profile_meta.to_parquet(output_root / "cache" / "family_state_gene_profile_metadata.parquet", index=False)
    np.savez_compressed(
        output_root / "cache" / "family_state_gene_profile_observed.npz",
        observed=profile_observed,
    )
    write_json(
        {"genes": profile_symbols},
        output_root / "cache" / "family_state_gene_profile_symbols.json",
    )
    write_yaml(architecture_surrogate_contract(config["version"]), output_root / "architecture_surrogate_contract.yaml")
    stage3_contract = {
        "contract": "v7_stage3_measurement_input",
        "version": config["version"],
        "spatial_execution_status": "NOT_RUN_PENDING_STAGE3",
        "D2_scope": "scRNA_measurement_only",
        "modality_evidence_contract": {
            "scRNA": "feature_level_D2",
            "bulkRNA": "GSE207422_two_exact_pairs_descriptive_only",
            "bulk_perturbation": "GSE193736_repeatability_interface_only_no_directional_validation",
            "spatial": "pending_real_projection_in_stage3",
        },
        "required_vocabulary": {
            "module_dictionary": "module_dictionary.tsv",
            "module_membership": "module_membership.tsv",
            "cell_state_ontology": "cell_state_ontology.tsv",
            "measurement_reliability": "measurement_reliability.tsv",
        },
        "levels": ["spot", "bin", "region", "patient"],
        "required_identity": ["patient_id", "block_id", "section_id", "observation_id", "leakage_group_id"],
        "required_measurement_fields": [
            "feature_id",
            "score_native",
            "score_standardized",
            "expected_genes",
            "observed_genes",
            "coverage",
            "technical_status",
            "D2_class",
            "qc_flags",
        ],
        "forbidden_stage2_claims": [
            "cross_modality_measurement_equivalence",
            "bulk_feature_validation_from_two_pairs",
            "perturbation_direction",
            "spatial_barrier",
            "TLS_structure",
            "NET_structure",
            "immune_exclusion",
            "PD1_failure",
        ],
    }
    write_yaml(stage3_contract, output_root / "STAGE3_MEASUREMENT_CONTRACT.yaml")


def _support_tables(config: dict[str, Any], output_root: Path) -> None:
    specs = load_object_specs(config["inputs"]["v6_projection_config"])
    overrides = config.get("study_family_overrides") or {}
    family = pd.DataFrame(
        [
            {
                "cohort_id": spec.cohort_id,
                "study_family": overrides.get(spec.cohort_id, spec.cohort_id),
                "object_id": spec.object_id,
                "modality": "scRNA",
                "family_rule": "explicit_override" if spec.cohort_id in overrides else "cohort_identity",
            }
            for spec in specs
        ]
        + [
            {"cohort_id": "GSE207422_bulk", "study_family": "GSE207422", "object_id": "GSE207422_bulk::matrix", "modality": "bulkRNA", "family_rule": "paired_source"},
            {"cohort_id": "GSE193736", "study_family": "GSE193736", "object_id": "GSE193736::bulkRNAseq_counts", "modality": "bulk_perturbation", "family_rule": "cohort_identity"},
        ]
    )
    write_tsv(family, output_root / "cohort_family_map.tsv", ["study_family", "cohort_id", "object_id"])
    exposure = family[["cohort_id", "study_family", "object_id", "modality"]].copy()
    exposure["stage2_claim"] = "response_blind_cross_modal_measurement"
    exposure["selection_exposure"] = "used_for_stage2_measurement_development_and_reliability"
    exposure["independent_validation_for_same_claim"] = False
    exposure["outcome_fields_read"] = False
    write_tsv(exposure, output_root / "selection_exposure.tsv", ["stage2_claim", "study_family", "cohort_id"])


def _validate(config: dict[str, Any], output_root: Path, *, require_complete: bool) -> dict[str, Any]:
    required = [
        "module_dictionary.tsv",
        "module_membership.tsv",
        "cell_state_ontology.tsv",
        "cell_state_label_mapping.tsv",
        "cross_modal_mapping.tsv",
        "measurement_reliability.tsv",
        "architecture_surrogate_contract.yaml",
        "STAGE3_MEASUREMENT_CONTRACT.yaml",
        "scores/expression_unit_scores.parquet",
        "scores/patient_timepoint_scores.parquet",
        "diagnostics/cell_level_technical_panel.tsv",
        "diagnostics/aggregation_audit.tsv",
        "diagnostics/legacy_fm_compatibility.tsv",
        "diagnostics/leave_study_family_out.tsv",
        "diagnostics/carrier_dependency.tsv",
        "diagnostics/ambient_observed_profile.tsv",
        "diagnostics/confounding_diagnostics.tsv",
        "diagnostics/method_agreement.tsv",
        "diagnostics/matched_expression_null_summary.tsv",
        "diagnostics/matched_expression_null_replicates.tsv",
        "diagnostics/raw_object_input_hashes.json",
        "diagnostics/gse207422_paired_scores.tsv",
        "diagnostics/gse207422_paired_feature_diagnostics.tsv",
        "diagnostics/gse207422_pair_vector_diagnostics.tsv",
        "diagnostics/gse193736_repeatability.tsv",
        "cache/family_state_gene_profile_observed.npz",
        "cache/family_state_gene_profile_symbols.json",
        "gene_mapping_audit.tsv",
        "gse207422_patient_timepoint_crosswalk.tsv",
        "gse193736_perturbation_design.tsv",
    ]
    missing = [name for name in required if not (output_root / name).exists()]
    if missing:
        raise Stage2Error("STAGE2_BLOCKED_MISSING_OUTPUT: " + ",".join(missing))
    audits = [json.loads(path.read_text(encoding="utf-8")) for path in sorted((output_root / "cache" / "pseudobulk").glob("*.audit.json"))]
    object_count = len({row["object_id"] for row in audits})
    if require_complete and object_count != int(config["execution"]["required_count_compatible_objects"]):
        raise Stage2Error(f"STAGE2_BLOCKED_PARTIAL_OBJECTS: {object_count}")
    expected_shards = int(config["execution"]["required_count_compatible_objects"]) * len(
        config["execution"]["primary_levels"]
    )
    if require_complete and len(audits) != expected_shards:
        raise Stage2Error(
            f"STAGE2_BLOCKED_PARTIAL_SHARDS: expected={expected_shards}, observed={len(audits)}"
        )
    observed_levels = {
        (row["object_id"], row["cell_state_level"]) for row in audits
    }
    expected_levels = {
        (object_id, level)
        for object_id in {row["object_id"] for row in audits}
        for level in config["execution"]["primary_levels"]
    }
    if require_complete and observed_levels != expected_levels:
        raise Stage2Error("STAGE2_BLOCKED_OBJECT_LEVEL_COVERAGE")
    if any(not row.get("counts_conserved") for row in audits):
        raise Stage2Error("STAGE2_BLOCKED_COUNTS_CONSERVATION")
    fingerprint_manifest = json.loads(
        (output_root / "diagnostics" / "raw_object_input_hashes.json").read_text(
            encoding="utf-8"
        )
    )
    expected_fingerprints = {
        str(row["object_id"]): str(row.get("fingerprint_sha256", ""))
        for row in fingerprint_manifest
    }
    observed_fingerprints: dict[str, set[str]] = {}
    for row in audits:
        fingerprint = row.get("input_fingerprint") or {}
        value = str(fingerprint.get("fingerprint_sha256", ""))
        observed_fingerprints.setdefault(str(row["object_id"]), set()).add(value)
    if require_complete and (
        set(observed_fingerprints) != set(expected_fingerprints)
        or any(len(values) != 1 or "" in values for values in observed_fingerprints.values())
        or any(
            next(iter(observed_fingerprints[object_id]))
            != expected_fingerprints[object_id]
            for object_id in expected_fingerprints
        )
    ):
        raise Stage2Error("STAGE2_BLOCKED_SHARD_INPUT_FINGERPRINT_MISMATCH")
    scores = pd.read_parquet(output_root / "scores" / "expression_unit_scores.parquet")
    patients = pd.read_parquet(output_root / "scores" / "patient_timepoint_scores.parquet")
    forbidden = prohibited_columns(scores.columns, config["prohibited_tokens"])
    forbidden_patients = prohibited_columns(patients.columns, config["prohibited_tokens"])
    if forbidden or forbidden_patients:
        raise Stage2Error(
            "STAGE2_BLOCKED_OUTPUT_LEAKAGE_COLUMNS: "
            + ",".join(sorted(set(forbidden + forbidden_patients)))
        )
    if scores.duplicated(["expression_unit_id", "feature_id"]).any():
        raise Stage2Error("STAGE2_BLOCKED_DUPLICATE_SCORE_KEY")
    reliability = pd.read_csv(output_root / "measurement_reliability.tsv", sep="\t")
    vocabulary = pd.read_csv(output_root / "module_dictionary.tsv", sep="\t")
    scoreable = set(
        vocabulary.loc[
            vocabulary.scoreable.astype(str).str.lower().isin({"true", "1", "yes"}),
            "feature_id",
        ].astype(str)
    )
    if set(scores.feature_id.astype(str)) != scoreable:
        raise Stage2Error("STAGE2_BLOCKED_SCORE_VOCABULARY_MISMATCH")
    score_counts = scores.groupby("expression_unit_id", sort=False).feature_id.nunique()
    if not score_counts.eq(len(scoreable)).all():
        raise Stage2Error("STAGE2_BLOCKED_INCOMPLETE_EXPRESSION_UNIT_FEATURES")
    patient_key = [
        "patient_timepoint_id",
        "cell_state_level",
        "feature_id",
        "aggregation_method",
    ]
    if patients.duplicated(patient_key).any():
        raise Stage2Error("STAGE2_BLOCKED_DUPLICATE_PATIENT_SCORE_KEY")
    if set(reliability.feature_id) != set(vocabulary.feature_id):
        raise Stage2Error("STAGE2_BLOCKED_RELIABILITY_VOCABULARY_MISMATCH")
    compatibility = pd.read_csv(output_root / "diagnostics" / "legacy_fm_compatibility.tsv", sep="\t")
    compatibility_pass = compatibility.status.astype(str).str.startswith("PASS_").all()
    formula_replay = compatibility.formula_replay_pass.astype(str).str.lower().isin(
        {"true", "1", "yes"}
    ).all()
    replay_columns = [
        "topic_mass_replay_pass",
        "native_score_replay_pass",
        "top_gene_sensitivity_replay_pass",
        "library_size_replay_pass",
    ]
    complete_replay = all(
        column in compatibility
        and compatibility[column].astype(str).str.lower().isin(
            {"true", "1", "yes"}
        ).all()
        for column in replay_columns
    )
    context_mapping = (
        "context_mapping_status" in compatibility
        and compatibility.context_mapping_status.eq(
            "PASS_COMMON_CONTEXT_PRESENT"
        ).all()
    )
    if require_complete and not (
        compatibility_pass and formula_replay and complete_replay and context_mapping
    ):
        raise Stage2Error("STAGE2_BLOCKED_LEGACY_NUMERIC_COMPATIBILITY")
    technical_panel = pd.read_csv(
        output_root / "diagnostics" / "cell_level_technical_panel.tsv", sep="\t"
    )
    expected_panel = set(map(str, config["execution"]["cell_level_technical_panel"]))
    observed_panel = set(technical_panel.cohort_id.astype(str))
    if require_complete and observed_panel != expected_panel:
        raise Stage2Error(
            "STAGE2_BLOCKED_TECHNICAL_PANEL_INCOMPLETE: "
            f"expected={sorted(expected_panel)}, observed={sorted(observed_panel)}"
        )
    expected_panel_rows = len(expected_panel) * len(config["execution"]["primary_levels"]) * len(
        scoreable
    )
    panel_key = ["cohort_id", "cell_state_level", "feature_id"]
    if require_complete and (
        len(technical_panel) != expected_panel_rows
        or technical_panel.duplicated(panel_key).any()
    ):
        raise Stage2Error(
            f"STAGE2_BLOCKED_TECHNICAL_PANEL_ROWS: expected={expected_panel_rows}, observed={len(technical_panel)}"
        )
    nulls = pd.read_csv(
        output_root / "diagnostics" / "matched_expression_null_summary.tsv", sep="\t"
    )
    if require_complete and (
        set(nulls.feature_id.astype(str)) != scoreable
        or not nulls.diagnostic_status.astype(str).str.startswith(
            ("PASS_", "NOT_TESTABLE_")
        ).all()
    ):
        raise Stage2Error("STAGE2_BLOCKED_MATCHED_NULL_INCOMPLETE")
    null_records = pd.read_csv(
        output_root / "diagnostics" / "matched_expression_null_replicates.tsv",
        sep="\t",
    )
    if not null_records.empty:
        mask_match = null_records.observation_mask_match.astype(str).str.lower().isin(
            {"true", "1", "yes"}
        )
        equal_gene_count = pd.to_numeric(
            null_records.n_target_genes, errors="coerce"
        ).eq(pd.to_numeric(null_records.n_null_genes, errors="coerce"))
        if require_complete and not (mask_match & equal_gene_count).all():
            raise Stage2Error("STAGE2_BLOCKED_MATCHED_NULL_MASK_MISMATCH")
    leave_out = pd.read_csv(
        output_root / "diagnostics" / "leave_study_family_out.tsv", sep="\t"
    )
    finite_lfo = (
        leave_out.assign(
            finite=pd.to_numeric(
                leave_out.state_profile_spearman, errors="coerce"
            ).notna()
        )
        .groupby("feature_id", sort=False)
        .finite.any()
    )
    if require_complete and (
        set(finite_lfo.index.astype(str)) != scoreable or not finite_lfo.all()
    ):
        raise Stage2Error("STAGE2_BLOCKED_LEAVE_FAMILY_OUT_INCOMPLETE")
    allowed_d2 = {
        "measurable",
        "joint_representation_with_uncertainty",
        "annotation_only_or_reject",
    }
    if not set(reliability.D2_class.astype(str)).issubset(allowed_d2):
        raise Stage2Error("STAGE2_BLOCKED_INVALID_D2_CLASS")
    if (
        not reliability.evidence_modality.eq("scRNA").all()
        or not reliability.D2_scope.eq("scRNA_measurement_only").all()
    ):
        raise Stage2Error("STAGE2_BLOCKED_D2_SCOPE")
    paired = pd.read_csv(
        output_root / "diagnostics" / "gse207422_paired_feature_diagnostics.tsv",
        sep="\t",
    )
    if require_complete and (
        set(paired.feature_id.astype(str)) != scoreable
        or not paired.n_exact_pairs.eq(2).all()
        or not paired.primary_pair_ids.eq("P05|P08").all()
        or not paired.paired_diagnostic_status.eq(
            "paired_descriptive_insufficient_n"
        ).all()
    ):
        raise Stage2Error("STAGE2_BLOCKED_GSE207422_PAIRED_DIAGNOSTIC")
    repeatability = pd.read_csv(
        output_root / "diagnostics" / "gse193736_repeatability.tsv", sep="\t"
    )
    if require_complete and set(repeatability.feature_id.astype(str)) != scoreable:
        raise Stage2Error("STAGE2_BLOCKED_GSE193736_REPEATABILITY_INCOMPLETE")
    cross_modal = pd.read_csv(output_root / "cross_modal_mapping.tsv", sep="\t")
    direct_mapping = cross_modal.loc[
        cross_modal.feature_id.astype(str).isin(scoreable)
        & cross_modal.modality.isin(["scRNA", "bulkRNA", "bulk_perturbation"])
    ]
    expected_modal_keys = {
        (feature_id, modality)
        for feature_id in scoreable
        for modality in ("scRNA", "bulkRNA", "bulk_perturbation")
    }
    observed_modal_keys = set(
        zip(
            direct_mapping.feature_id.astype(str),
            direct_mapping.modality.astype(str),
        )
    )
    bulk_classes = set(
        direct_mapping.loc[direct_mapping.modality.eq("bulkRNA"), "D2_class"]
    )
    perturb_classes = set(
        direct_mapping.loc[
            direct_mapping.modality.eq("bulk_perturbation"), "D2_class"
        ]
    )
    if require_complete and (
        observed_modal_keys != expected_modal_keys
        or bulk_classes != {"paired_descriptive_insufficient_n"}
        or perturb_classes != {"interface_only_no_directional_validation"}
    ):
        raise Stage2Error("STAGE2_BLOCKED_MODALITY_SPECIFIC_D2")
    return {
        "status": "VALID" if (not require_complete or object_count == int(config["execution"]["required_count_compatible_objects"])) else "PARTIAL",
        "n_objects": object_count,
        "n_shards": len(audits),
        "n_score_rows": len(scores),
        "n_features": len(vocabulary),
        "n_reliability_rows": len(reliability),
        "legacy_compatibility": bool(compatibility_pass),
        "technical_panel_cohorts": len(observed_panel),
    }


def _finalize(config: dict[str, Any], output_root: Path) -> dict[str, Any]:
    specs = load_object_specs(config["inputs"]["v6_projection_config"])
    overrides = config.get("study_family_overrides") or {}
    semantic_context = _aggregation_semantic_context(config)
    previous_manifest = _previous_hash_manifest(output_root)
    fingerprints = [
        shard_input_fingerprint(
            spec,
            study_family=overrides.get(spec.cohort_id, spec.cohort_id),
            semantic_context=semantic_context,
            previous=previous_manifest.get(spec.object_id),
        )
        for spec in specs
    ]
    write_json(
        fingerprints,
        output_root / "diagnostics" / "raw_object_input_hashes.json",
    )
    _support_tables(config, output_root)
    validation = _validate(config, output_root, require_complete=True)
    reliability = pd.read_csv(output_root / "measurement_reliability.tsv", sep="\t")
    counts = reliability.D2_class.value_counts().to_dict()
    scoreable = reliability.loc[reliability.feature_family.ne("mechanism_axis")]
    parents = reliability.loc[reliability.feature_family.eq("mechanism_axis")]
    scoreable_counts = scoreable.D2_class.value_counts().to_dict()
    parent_counts = parents.D2_class.value_counts().to_dict()
    technical_panel = pd.read_csv(
        output_root / "diagnostics" / "cell_level_technical_panel.tsv", sep="\t"
    )
    panel_median = float(
        pd.to_numeric(
            technical_panel.cell_mean_vs_pseudobulk_spearman, errors="coerce"
        ).median()
    )
    null_summary = pd.read_csv(
        output_root / "diagnostics" / "matched_expression_null_summary.tsv",
        sep="\t",
    )
    null_status_counts = null_summary.diagnostic_status.value_counts().to_dict()
    paired_feature = pd.read_csv(
        output_root / "diagnostics" / "gse207422_paired_feature_diagnostics.tsv",
        sep="\t",
    )
    paired_vectors = pd.read_csv(
        output_root / "diagnostics" / "gse207422_pair_vector_diagnostics.tsv",
        sep="\t",
    )
    repeatability = pd.read_csv(
        output_root / "diagnostics" / "gse193736_repeatability.tsv",
        sep="\t",
    )
    report = f"""# v7 Stage 2：共同生物学词汇与跨模态测量

## 技术路线

- 固定继承 8 个 legacy FM，并并行建立 10 个机制父轴、31 个显式组件。
- 对 44 个 count-compatible 单细胞对象从 raw counts 重建患者/样本×时间点×coarse/mid state pseudobulk。
- 在 GSE207422 建立 response-blind scRNA–bulk 桥，在 GSE193736 建立扰动 bulk 技术桥。
- 在 TASK01、GSE207422_sc、GSE301741 对全部合格细胞执行逐细胞技术计分，并与 count pseudobulk 的排序一致性对照；不导出全细胞分数表。
- 通过 exact-observation-mask matched-expression null、leave-study-family-out、coverage、native ambient proxy、carrier、method agreement 和 confounding 诊断确定 scRNA 测量等级。null 候选不足时保留 `NOT_TESTABLE`，不把未测基因当作 0。
- 空间 spot/bin/region 仅输出 Stage 3 接口，本阶段未执行真实空间投影。

## 结果 / 证据

- 完整处理对象：{validation['n_objects']}；pseudobulk 分片：{validation['n_shards']}；计数守恒全部通过。
- 可执行测量特征：39 个（8 FM + 31 机制组件）；另有 10 个机制父轴按组件向量解释。
- 39 个可执行特征的 scRNA D2 分布：{json.dumps(scoreable_counts, ensure_ascii=False, sort_keys=True)}；10 个机制父轴：{json.dumps(parent_counts, ensure_ascii=False, sort_keys=True)}。D2 仅表示 scRNA 测量证据，不外推到其他模态。
- exact-mask null 诊断状态：{json.dumps(null_status_counts, ensure_ascii=False, sort_keys=True)}；具体不可估计的特征保留为不确定性，不强行使用。
- 8 个 legacy FM 均通过当前输入上的冻结公式重放。旧 v6 行级数值因上游 cell-to-context 分配更新而不再要求逐行相等，已作为 historical context drift 明示记录。
- 逐细胞技术面板完成队列：{validation['technical_panel_cohorts']}；234 个“队列×层级×特征”比较的 rank concordance 中位数为 {panel_median:.3f}，低或负一致性的组件保留为 D2 uncertainty。
- GSE207422 精确患者—时间点配对仅 P05、P08；P07 为时间点不一致的支持性桥接。P05/P08 的 39 维路线 Spearman 分别为 {paired_vectors.loc[paired_vectors.comparison_id.eq('P05'), 'feature_vector_spearman'].iloc[0]:.3f} 和 {paired_vectors.loc[paired_vectors.comparison_id.eq('P08'), 'feature_vector_spearman'].iloc[0]:.3f}，只作接口描述；39 个特征的患者级配对校准都是 `paired_descriptive_insufficient_n`。
- GSE193736 的 24 个扰动条件/重复样本均可使用同一测量接口；39 个特征的重复谱 median Spearman 为 {float(pd.to_numeric(repeatability.replicate_profile_median_spearman, errors='coerce').median()):.3f}，但不据此作扰动方向或疗效结论。
- ambient 诊断使用 low-quality/ambient 标签作观察性 proxy，不等同于已证明污染因果；platform 效应在本阶段单一 scRNA platform 下不可估计。

## 结论

Stage 2 已建立可追溯、response-blind 的共同生物学词汇。scRNA 的 D2 测量基础已可供后续空间与患者级分析；bulk 仅有两位精确配对的描述性支持，perturbation 仅有重复性接口证据，不是已验证的等价模态。空间屏障、PD-1 失败和 PD1+X 修复仍需后续阶段的数据与分析。

## 风险与限制

- 机制父轴是多组件向量，不把方向不同的生物学强行压成单值。
- D2 uncertainty 是测量证据的一部分；不能把带不确定性的程序写成已确认机制。matched-null 的区间是算法稳定/置换包络，不是生物效应置信区间。
- 逐细胞先归一化与 counts 先汇总再归一化是不同估计量；技术面板只检验排序稳健性，不要求数值相等。
- GSE207422 的精确配对患者数很少，只支持接口可行性。
- 真实空间投影状态为 `NOT_RUN_PENDING_STAGE3`。
"""
    (output_root / "ONTOLOGY_AND_MEASUREMENT_REPORT.md").write_text(report, encoding="utf-8")
    tracked_outputs = [
        output_root / name
        for name in [
            "module_dictionary.tsv",
            "module_membership.tsv",
            "cell_state_ontology.tsv",
            "cell_state_label_mapping.tsv",
            "cross_modal_mapping.tsv",
            "measurement_reliability.tsv",
            "architecture_surrogate_contract.yaml",
            "STAGE3_MEASUREMENT_CONTRACT.yaml",
            "cohort_family_map.tsv",
            "selection_exposure.tsv",
            "gene_mapping_audit.tsv",
            "ONTOLOGY_AND_MEASUREMENT_REPORT.md",
            "diagnostics/cell_level_technical_panel.tsv",
            "diagnostics/aggregation_audit.tsv",
            "diagnostics/input_audit.tsv",
            "diagnostics/leakage_preflight.tsv",
            "diagnostics/legacy_fm_compatibility.tsv",
            "diagnostics/leave_study_family_out.tsv",
            "diagnostics/carrier_dependency.tsv",
            "diagnostics/ambient_observed_profile.tsv",
            "diagnostics/confounding_diagnostics.tsv",
            "diagnostics/method_agreement.tsv",
            "diagnostics/matched_expression_null_summary.tsv",
            "diagnostics/matched_expression_null_replicates.tsv",
            "diagnostics/gse207422_paired_scores.tsv",
            "diagnostics/gse207422_paired_feature_diagnostics.tsv",
            "diagnostics/gse207422_pair_vector_diagnostics.tsv",
            "diagnostics/gse193736_repeatability.tsv",
            "diagnostics/raw_object_input_hashes.json",
            "cache/family_state_gene_profile_observed.npz",
            "cache/family_state_gene_profile_symbols.json",
            "gse207422_patient_timepoint_crosswalk.tsv",
            "gse193736_perturbation_design.tsv",
            "scores/expression_unit_scores.parquet",
            "scores/patient_timepoint_scores.parquet",
        ]
    ]
    manifest = {
        "stage": "v7_stage2_ontology_and_measurement",
        "version": config["version"],
        "status": "STAGE2_COMPLETE",
        "response_blind": True,
        "spatial_execution_status": "NOT_RUN_PENDING_STAGE3",
        "D2_scope": "scRNA_measurement_only",
        "modality_evidence_scope": {
            "scRNA": "feature_level_D2",
            "bulkRNA": "paired_descriptive_only_two_exact_pairs",
            "bulk_perturbation": "repeatability_interface_only_no_directional_validation",
            "spatial": "pending_stage3",
        },
        "validation": validation,
        "D2_counts": counts,
        "raw_object_hash_manifest": "diagnostics/raw_object_input_hashes.json",
        "outputs": content_manifest(tracked_outputs),
        "prohibitions": [
            "response_guided_vocabulary",
            "legacy_FM_retraining",
            "missing_gene_as_zero",
            "spatial_barrier_claim",
            "PD1_failure_claim",
        ],
    }
    write_yaml(manifest, output_root / "STAGE2_RUN_MANIFEST.yaml")
    return validation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/v7/stage2.yaml")
    parser.add_argument(
        "--step",
        choices=[
            "preflight",
            "vocabulary",
            "aggregate",
            "score",
            "technical-panel",
            "diagnose",
            "finalize",
            "full",
        ],
        default="full",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--cohort")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    try:
        config = load_config(args.config)
        output_root = config["output_root"]
        output_root.mkdir(parents=True, exist_ok=True)
        if args.validate_only:
            print(json.dumps(_validate(config, output_root, require_complete=True), sort_keys=True))
            return 0
        steps = (
            ["preflight", "vocabulary", "aggregate", "score", "technical-panel", "diagnose", "finalize"]
            if args.step == "full"
            else [args.step]
        )
        for step in steps:
            print(json.dumps({"event": "stage2_step_start", "step": step}), flush=True)
            if step == "preflight":
                _preflight(config, output_root)
            elif step == "vocabulary":
                _vocabulary(config, output_root)
            elif step == "aggregate":
                _aggregate(config, output_root, resume=args.resume, cohort=args.cohort, workers=args.workers)
            elif step == "score":
                _score(config, output_root)
            elif step == "technical-panel":
                _technical_panel(config, output_root)
            elif step == "diagnose":
                _diagnose(config, output_root)
            elif step == "finalize":
                result = _finalize(config, output_root)
                print(json.dumps(result, sort_keys=True), flush=True)
            print(json.dumps({"event": "stage2_step_complete", "step": step}), flush=True)
        return 0
    except (Stage2Error, ValueError, OSError) as error:
        print(f"STAGE2_BLOCKED: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
