"""Pure, response-blind compilers for the v7 Stage 2 vocabulary.

The functions in this module only read their declared inputs and return
deterministically ordered data frames.  Writing and orchestration belong to the
Stage 2 runner.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd
import yaml

from .contracts import Stage2Error, sha256


FROZEN_FM_MEMBERSHIP_SHA256 = (
    "285c11de4b071c45b0ec43464fb963650340c29bb2de3037570578b29ee44036"
)
FROZEN_FM_IDS = tuple(f"FM{index:02d}" for index in range(1, 9))


MODULE_DICTIONARY_COLUMNS = [
    "feature_id",
    "feature_family",
    "feature_level",
    "parent_axis_id",
    "feature_name",
    "display_name",
    "component_role",
    "scoring_semantics",
    "signed",
    "scoreable",
    "n_genes",
    "source_resource",
    "source_set",
    "source_sha256",
    "provenance_status",
    "legacy_primary_method",
    "legacy_primary_run",
    "legacy_primary_module_id",
    "support_tier",
    "n_supporting_methods",
    "bootstrap_median_top40_jaccard",
    "leakage_detected",
    "default_for_phase7",
    "notes",
]

MODULE_MEMBERSHIP_COLUMNS = [
    "feature_id",
    "feature_family",
    "parent_axis_id",
    "gene_symbol",
    "direction",
    "weight",
    "relative_weight",
    "is_top_gene",
    "source_resource",
    "source_set",
    "source_record",
    "source_sha256",
    "provenance_status",
]

ONTOLOGY_COLUMNS = [
    "ontology_id",
    "level",
    "node_kind",
    "label",
    "parent_label",
    "parent_ontology_id",
    "definition",
    "reliability_rule",
    "allowed_downstream_use",
    "is_exclusive",
    "is_synthetic",
    "legacy_annotation_role",
    "source_version",
    "source_sha256",
    "provenance_status",
]

LABEL_MAPPING_COLUMNS = [
    "mapping_id",
    "cohort_id",
    "annotation_field",
    "original_label",
    "harmonized_coarse",
    "harmonized_mid",
    "harmonized_fine",
    "coarse_ontology_id",
    "mid_ontology_id",
    "fine_ontology_id",
    "mapping_confidence",
    "mapping_status",
    "mapping_resolution_status",
    "shadow_candidate_level",
    "shadow_candidate_label",
    "shadow_candidate_ontology_id",
    "n_shadow_observations",
    "mapping_source",
    "source_sha256",
]


@dataclass(frozen=True)
class Stage2Vocabulary:
    """All Stage 2 vocabulary tables, held in memory."""

    module_dictionary: pd.DataFrame
    module_membership: pd.DataFrame
    cell_state_ontology: pd.DataFrame
    label_mapping: pd.DataFrame


@dataclass(frozen=True)
class _ComponentSpec:
    axis_slug: str
    component_slug: str
    source_kind: str
    source_set: str
    role: str = "primary"


_AXES: tuple[tuple[str, str, str], ...] = (
    (
        "angiogenesis_endothelial",
        "VEGF / angiogenesis / endothelial",
        "Molecular angiogenesis and endothelial activation; not spatial vascular architecture.",
    ),
    (
        "caf_ecm",
        "CAF / ECM",
        "Stromal and extracellular-matrix programs; not spatial immune exclusion by itself.",
    ),
    ("hypoxia", "Hypoxia", "Cellular hypoxia programs."),
    (
        "suppressive_myeloid",
        "Suppressive myeloid",
        "Suppressive-myeloid-like expression; not a causal suppressor-cell claim.",
    ),
    (
        "apc_ifn_apm",
        "APC / IFN / antigen presentation",
        "Antigen presentation and interferon components retained separately.",
    ),
    (
        "cytotoxic_t_nk",
        "Cytotoxic T / NK",
        "T- and NK-cell cytotoxic effector programs.",
    ),
    (
        "exhaustion",
        "Progenitor / terminal exhaustion",
        "Progenitor, terminal and broad exhaustion are separate non-exclusive components.",
    ),
    (
        "tls_b",
        "TLS / B-like",
        "TLS/B-like molecular program; Stage 2 does not establish a TLS structure.",
    ),
    (
        "neutrophil_tan_net",
        "Neutrophil / TAN / NET-like",
        "Neutrophil and NET-like expression; Stage 2 does not establish NET topology.",
    ),
    (
        "wnt_beta_catenin",
        "WNT / beta-catenin",
        "WNT/beta-catenin expression; Stage 2 does not establish spatial exclusion.",
    ),
)


_COMPONENTS: tuple[_ComponentSpec, ...] = (
    _ComponentSpec("angiogenesis_endothelial", "curated_angiogenesis", "curated", "angiogenesis"),
    _ComponentSpec("angiogenesis_endothelial", "hallmark_angiogenesis", "hallmark", "HALLMARK_ANGIOGENESIS"),
    _ComponentSpec("angiogenesis_endothelial", "reactome_vegf", "reactome", "REACTOME_SIGNALING_BY_VEGF"),
    _ComponentSpec("caf_ecm", "curated_stromal_tgfb", "curated", "stromal_TGFb"),
    _ComponentSpec("caf_ecm", "hallmark_emt", "hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
    _ComponentSpec("caf_ecm", "reactome_ecm", "reactome", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"),
    _ComponentSpec("hypoxia", "curated_hypoxia", "curated", "hypoxia"),
    _ComponentSpec("hypoxia", "hallmark_hypoxia", "hallmark", "HALLMARK_HYPOXIA"),
    _ComponentSpec("suppressive_myeloid", "curated_myeloid_suppression", "curated", "myeloid_suppression"),
    _ComponentSpec("apc_ifn_apm", "curated_apc_dc_activation", "curated", "APC_DC_activation"),
    _ComponentSpec("apc_ifn_apm", "curated_ifn_response", "curated", "IFN_response"),
    _ComponentSpec("apc_ifn_apm", "curated_antigen_presentation", "curated", "antigen_presentation"),
    _ComponentSpec("apc_ifn_apm", "hallmark_ifn_alpha", "hallmark", "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
    _ComponentSpec("apc_ifn_apm", "hallmark_ifn_gamma", "hallmark", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    _ComponentSpec("apc_ifn_apm", "reactome_mhc_i", "reactome", "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION"),
    _ComponentSpec("apc_ifn_apm", "reactome_mhc_ii", "reactome", "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"),
    _ComponentSpec("apc_ifn_apm", "reactome_ifn_alpha_beta", "reactome", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"),
    _ComponentSpec("apc_ifn_apm", "reactome_ifn_gamma", "reactome", "REACTOME_INTERFERON_GAMMA_SIGNALING"),
    _ComponentSpec("cytotoxic_t_nk", "curated_cytotoxicity", "curated", "cytotoxicity"),
    _ComponentSpec("cytotoxic_t_nk", "curated_nk_cytotoxicity", "curated", "NK_cytotoxicity"),
    _ComponentSpec("cytotoxic_t_nk", "hallmark_allograft_rejection", "hallmark", "HALLMARK_ALLOGRAFT_REJECTION", "support"),
    _ComponentSpec("exhaustion", "curated_broad_exhaustion", "curated", "exhaustion"),
    _ComponentSpec("exhaustion", "sentinel_progenitor_exhaustion", "sentinel", "CD8_PROGENITOR_EXHAUSTION"),
    _ComponentSpec("exhaustion", "sentinel_terminal_exhaustion", "sentinel", "CD8_TERMINAL_EXHAUSTION"),
    _ComponentSpec("tls_b", "sentinel_tls_b", "sentinel", "TLS_B"),
    _ComponentSpec("tls_b", "reactome_bcr", "reactome", "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR", "support"),
    _ComponentSpec("neutrophil_tan_net", "sentinel_neutrophil_net", "sentinel", "NEUTROPHIL_NET"),
    _ComponentSpec("neutrophil_tan_net", "reactome_neutrophil_degranulation", "reactome", "REACTOME_NEUTROPHIL_DEGRANULATION", "support"),
    _ComponentSpec("wnt_beta_catenin", "sentinel_tumor_wnt", "sentinel", "TUMOR_WNT_EXCLUSION"),
    _ComponentSpec("wnt_beta_catenin", "hallmark_wnt_beta_catenin", "hallmark", "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
    _ComponentSpec("wnt_beta_catenin", "reactome_wnt", "reactome", "REACTOME_SIGNALING_BY_WNT"),
)


def _require_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise Stage2Error(f"STAGE2_BLOCKED_SCHEMA: {label} missing columns {','.join(missing)}")


def _stable(frame: pd.DataFrame, columns: Sequence[str], sort_by: Sequence[str]) -> pd.DataFrame:
    material = frame.reindex(columns=columns).copy()
    if not material.empty:
        material = material.sort_values(list(sort_by), kind="mergesort", na_position="last")
    return material.reset_index(drop=True)


def _boolean(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _gene(value: object) -> str:
    gene = _text(value)
    if not gene:
        raise Stage2Error("STAGE2_BLOCKED_VOCABULARY: empty gene symbol")
    return gene


def _mapping_id(parts: Iterable[object]) -> str:
    payload = "\x1f".join(_text(part) for part in parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:20]


def compile_legacy_modules(
    membership_path: str | Path,
    dictionary_path: str | Path,
    expected_sha256: str = FROZEN_FM_MEMBERSHIP_SHA256,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compile and verify the immutable v6 FM01--FM08 baseline."""

    membership_file = Path(membership_path)
    dictionary_file = Path(dictionary_path)
    observed_sha = sha256(membership_file)
    if observed_sha != expected_sha256:
        raise Stage2Error(
            "STAGE2_BLOCKED_FROZEN_FM_HASH: "
            f"expected={expected_sha256} observed={observed_sha}"
        )

    membership = pd.read_csv(membership_file)
    dictionary = pd.read_csv(dictionary_file)
    _require_columns(
        membership,
        [
            "frozen_module_id",
            "primary_method",
            "primary_run",
            "primary_module_id",
            "gene",
            "membership_weight",
            "relative_weight",
            "is_top_gene",
        ],
        "legacy FM membership",
    )
    _require_columns(
        dictionary,
        [
            "frozen_module_id",
            "primary_method",
            "primary_run",
            "primary_module_id",
            "support_tier",
            "n_supporting_methods",
            "bootstrap_median_top40_jaccard",
            "leakage_detected",
            "default_for_phase7",
        ],
        "legacy FM dictionary",
    )

    observed_ids = tuple(sorted(membership["frozen_module_id"].astype(str).unique()))
    dictionary_ids = tuple(sorted(dictionary["frozen_module_id"].astype(str).unique()))
    if observed_ids != FROZEN_FM_IDS or dictionary_ids != FROZEN_FM_IDS:
        raise Stage2Error(
            "STAGE2_BLOCKED_FROZEN_FM_IDS: expected FM01-FM08 in membership and dictionary"
        )
    if dictionary["frozen_module_id"].duplicated().any():
        raise Stage2Error("STAGE2_BLOCKED_FROZEN_FM_DICTIONARY: duplicate module id")
    if membership.duplicated(["frozen_module_id", "gene"]).any():
        raise Stage2Error("STAGE2_BLOCKED_FROZEN_FM_MEMBERSHIP: duplicate module-gene row")

    for column in ("membership_weight", "relative_weight"):
        membership[column] = pd.to_numeric(membership[column], errors="raise")
        if membership[column].isna().any() or (membership[column] < 0).any():
            raise Stage2Error(f"STAGE2_BLOCKED_FROZEN_FM_MEMBERSHIP: invalid {column}")
    if dictionary["leakage_detected"].map(_boolean).any():
        raise Stage2Error("STAGE2_BLOCKED_FROZEN_FM_LEAKAGE: frozen dictionary reports leakage")

    membership_records: list[dict[str, object]] = []
    for row in membership.to_dict("records"):
        membership_records.append(
            {
                "feature_id": _text(row["frozen_module_id"]),
                "feature_family": "legacy_fm",
                "parent_axis_id": "",
                "gene_symbol": _gene(row["gene"]),
                "direction": 1,
                "weight": float(row["membership_weight"]),
                "relative_weight": float(row["relative_weight"]),
                "is_top_gene": _boolean(row["is_top_gene"]),
                "source_resource": "v6_frozen_lda_topic_membership",
                "source_set": _text(row["primary_run"]),
                "source_record": _text(row["primary_module_id"]),
                "source_sha256": observed_sha,
                "provenance_status": "immutable_legacy_baseline",
            }
        )
    output_membership = _stable(
        pd.DataFrame(membership_records),
        MODULE_MEMBERSHIP_COLUMNS,
        ["feature_id", "gene_symbol"],
    )

    member_counts = output_membership.groupby("feature_id").size().to_dict()
    dictionary_records: list[dict[str, object]] = []
    for row in dictionary.to_dict("records"):
        feature_id = _text(row["frozen_module_id"])
        dictionary_records.append(
            {
                "feature_id": feature_id,
                "feature_family": "legacy_fm",
                "feature_level": "component",
                "parent_axis_id": "",
                "feature_name": feature_id,
                "display_name": feature_id,
                "component_role": "legacy_comparability_baseline",
                "scoring_semantics": "weighted_topic_mass_per_million_log1p",
                "signed": False,
                "scoreable": True,
                "n_genes": int(member_counts[feature_id]),
                "source_resource": "v6_frozen_lda_topic_membership",
                "source_set": _text(row["primary_run"]),
                "source_sha256": observed_sha,
                "provenance_status": "immutable_legacy_baseline_remeasure_only",
                "legacy_primary_method": _text(row["primary_method"]),
                "legacy_primary_run": _text(row["primary_run"]),
                "legacy_primary_module_id": _text(row["primary_module_id"]),
                "support_tier": _text(row["support_tier"]),
                "n_supporting_methods": int(row["n_supporting_methods"]),
                "bootstrap_median_top40_jaccard": float(
                    row["bootstrap_median_top40_jaccard"]
                ),
                "leakage_detected": _boolean(row["leakage_detected"]),
                "default_for_phase7": _boolean(row["default_for_phase7"]),
                "notes": "Legacy annotations, if added downstream, remain provisional.",
            }
        )
    output_dictionary = _stable(
        pd.DataFrame(dictionary_records),
        MODULE_DICTIONARY_COLUMNS,
        ["feature_id"],
    )
    return output_dictionary, output_membership


def _read_gmt(path: Path) -> Mapping[str, tuple[str, ...]]:
    sets: dict[str, tuple[str, ...]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 3:
                raise Stage2Error(
                    f"STAGE2_BLOCKED_GMT: {path}:{line_number} has fewer than three fields"
                )
            name = fields[0].strip()
            if name in sets:
                raise Stage2Error(f"STAGE2_BLOCKED_GMT: duplicate set {name} in {path}")
            genes = tuple(sorted({_gene(gene) for gene in fields[2:] if _text(gene)}))
            if not genes:
                raise Stage2Error(f"STAGE2_BLOCKED_GMT: empty set {name} in {path}")
            sets[name] = genes
    return sets


def _curated_components(path: Path) -> Mapping[str, tuple[tuple[str, int], ...]]:
    frame = pd.read_csv(path)
    _require_columns(
        frame,
        ["feature_name", "gene_symbol", "direction", "source_resource"],
        "v6 curated signature registry",
    )
    # The registry contains response-derived MVP records.  A row is eligible
    # only when the curated core is an explicit source; mixed rows are retained
    # solely through their curated provenance.
    curated = frame[
        frame["source_resource"].astype(str).str.contains(
            "curated_step2_9_core", regex=False, na=False
        )
    ].copy()
    if curated.empty:
        raise Stage2Error("STAGE2_BLOCKED_CURATED_SIGNATURES: no curated core rows")
    curated["direction"] = pd.to_numeric(curated["direction"], errors="raise").astype(int)
    if not curated["direction"].isin([-1, 1]).all():
        raise Stage2Error("STAGE2_BLOCKED_CURATED_SIGNATURES: direction must be -1 or 1")

    output: dict[str, tuple[tuple[str, int], ...]] = {}
    for name, group in curated.groupby("feature_name", sort=True):
        records = sorted({(_gene(row.gene_symbol), int(row.direction)) for row in group.itertuples()})
        by_gene: dict[str, set[int]] = {}
        for gene, direction in records:
            by_gene.setdefault(gene, set()).add(direction)
        conflicts = sorted(gene for gene, directions in by_gene.items() if len(directions) > 1)
        if conflicts:
            raise Stage2Error(
                f"STAGE2_BLOCKED_CURATED_SIGNATURES: opposing directions for {name}: "
                + ",".join(conflicts)
            )
        output[str(name)] = tuple(records)
    return output


def compile_mechanism_axes(
    signature_registry_path: str | Path,
    hallmark_gmt_path: str | Path,
    reactome_gmt_path: str | Path,
    sentinel_config_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compile ten mechanism parents and their explicit score components."""

    signature_file = Path(signature_registry_path)
    hallmark_file = Path(hallmark_gmt_path)
    reactome_file = Path(reactome_gmt_path)
    sentinel_file = Path(sentinel_config_path)
    curated = _curated_components(signature_file)
    hallmark = _read_gmt(hallmark_file)
    reactome = _read_gmt(reactome_file)
    sentinel_payload = yaml.safe_load(sentinel_file.read_text(encoding="utf-8")) or {}
    sentinels = sentinel_payload.get("sentinel_programs") or {}

    source_hashes = {
        "curated": sha256(signature_file),
        "hallmark": sha256(hallmark_file),
        "reactome": sha256(reactome_file),
        "sentinel": sha256(sentinel_file),
    }
    source_resources = {
        "curated": "v6_curated_step2_9_core",
        "hallmark": "MSigDB_Hallmark_v2025.1",
        "reactome": "MSigDB_Reactome_v2025.1",
        "sentinel": "v6_2_response_blind_sentinel_config",
    }

    axis_rows: list[dict[str, object]] = []
    for slug, display_name, notes in _AXES:
        axis_rows.append(
            {
                "feature_id": f"axis__{slug}",
                "feature_family": "mechanism_axis",
                "feature_level": "parent_axis",
                "parent_axis_id": "",
                "feature_name": slug,
                "display_name": display_name,
                "component_role": "conceptual_parent",
                "scoring_semantics": "vector_of_explicit_components",
                "signed": False,
                "scoreable": False,
                "n_genes": 0,
                "source_resource": "v7_stage2_explicit_axis_definition",
                "source_set": slug,
                "source_sha256": "",
                "provenance_status": "response_blind_conceptual_parent",
                "legacy_primary_method": "",
                "legacy_primary_run": "",
                "legacy_primary_module_id": "",
                "support_tier": "",
                "n_supporting_methods": 0,
                "bootstrap_median_top40_jaccard": "",
                "leakage_detected": False,
                "default_for_phase7": False,
                "notes": notes,
            }
        )

    dictionary_rows: list[dict[str, object]] = axis_rows
    membership_rows: list[dict[str, object]] = []
    for spec in _COMPONENTS:
        if spec.source_kind == "curated":
            if spec.source_set not in curated:
                raise Stage2Error(
                    f"STAGE2_BLOCKED_MECHANISM_COMPONENT: missing curated set {spec.source_set}"
                )
            members = curated[spec.source_set]
        elif spec.source_kind in {"hallmark", "reactome"}:
            collection = hallmark if spec.source_kind == "hallmark" else reactome
            if spec.source_set not in collection:
                raise Stage2Error(
                    f"STAGE2_BLOCKED_MECHANISM_COMPONENT: missing GMT set {spec.source_set}"
                )
            members = tuple((gene, 1) for gene in collection[spec.source_set])
        elif spec.source_kind == "sentinel":
            payload = sentinels.get(spec.source_set) or {}
            genes = payload.get("genes") if isinstance(payload, Mapping) else None
            if not genes:
                raise Stage2Error(
                    f"STAGE2_BLOCKED_MECHANISM_COMPONENT: missing sentinel {spec.source_set}"
                )
            members = tuple((_gene(gene), 1) for gene in sorted(set(genes)))
        else:  # pragma: no cover - module-owned specification invariant
            raise AssertionError(spec.source_kind)

        members = tuple(sorted(set(members), key=lambda item: (item[0], item[1])))
        directions_by_gene: dict[str, set[int]] = {}
        for gene, direction in members:
            directions_by_gene.setdefault(gene, set()).add(int(direction))
        conflicts = sorted(gene for gene, directions in directions_by_gene.items() if len(directions) > 1)
        if conflicts:
            raise Stage2Error(
                f"STAGE2_BLOCKED_MECHANISM_COMPONENT: opposing directions in {spec.source_set}: "
                + ",".join(conflicts)
            )

        feature_id = f"mechanism__{spec.axis_slug}__{spec.component_slug}"
        parent_id = f"axis__{spec.axis_slug}"
        is_signed = any(direction < 0 for _, direction in members)
        dictionary_rows.append(
            {
                "feature_id": feature_id,
                "feature_family": "mechanism_component",
                "feature_level": "component",
                "parent_axis_id": parent_id,
                "feature_name": spec.component_slug,
                "display_name": spec.source_set,
                "component_role": spec.role,
                "scoring_semantics": (
                    "signed_mean_log_expression"
                    if is_signed
                    else "unsigned_mean_log_expression"
                ),
                "signed": is_signed,
                "scoreable": True,
                "n_genes": len(members),
                "source_resource": source_resources[spec.source_kind],
                "source_set": spec.source_set,
                "source_sha256": source_hashes[spec.source_kind],
                "provenance_status": "response_blind_curated_or_public_pathway",
                "legacy_primary_method": "",
                "legacy_primary_run": "",
                "legacy_primary_module_id": "",
                "support_tier": "",
                "n_supporting_methods": 0,
                "bootstrap_median_top40_jaccard": "",
                "leakage_detected": False,
                "default_for_phase7": False,
                "notes": "Components remain separate; the parent axis is not a collapsed scalar.",
            }
        )
        relative_weight = 1.0 / len(members)
        for gene, direction in members:
            membership_rows.append(
                {
                    "feature_id": feature_id,
                    "feature_family": "mechanism_component",
                    "parent_axis_id": parent_id,
                    "gene_symbol": gene,
                    "direction": int(direction),
                    "weight": 1.0,
                    "relative_weight": relative_weight,
                    "is_top_gene": False,
                    "source_resource": source_resources[spec.source_kind],
                    "source_set": spec.source_set,
                    "source_record": spec.source_set,
                    "source_sha256": source_hashes[spec.source_kind],
                    "provenance_status": "response_blind_curated_or_public_pathway",
                }
            )

    dictionary = _stable(
        pd.DataFrame(dictionary_rows),
        MODULE_DICTIONARY_COLUMNS,
        ["feature_family", "parent_axis_id", "feature_id"],
    )
    membership = _stable(
        pd.DataFrame(membership_rows),
        MODULE_MEMBERSHIP_COLUMNS,
        ["parent_axis_id", "feature_id", "gene_symbol"],
    )
    if dictionary.loc[dictionary["feature_level"] == "parent_axis", "feature_id"].nunique() != 10:
        raise Stage2Error("STAGE2_BLOCKED_MECHANISM_AXES: expected ten parent axes")
    return dictionary, membership


_ANNOTATION_SENTINELS = ("Unknown", "Mixed", "Low_quality_or_ambient")


def _ontology_id(level: str, label: str, node_kind: str) -> str:
    return f"{node_kind}::{level}::{label}"


def _compile_ontology_table(path: Path) -> tuple[pd.DataFrame, Mapping[tuple[str, str], str]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    levels = payload.get("levels") or {}
    if not all(level in levels and isinstance(levels[level], Mapping) for level in ("coarse", "mid", "fine")):
        raise Stage2Error("STAGE2_BLOCKED_ONTOLOGY: coarse/mid/fine mappings are required")
    source_version = _text(payload.get("ontology_version")) or "unknown"
    source_hash = sha256(path)
    hooks = set(payload.get("b2_annotation_hooks") or [])

    material: dict[str, dict[str, object]] = {
        level: {str(label): dict(details or {}) for label, details in levels[level].items()}
        for level in ("coarse", "mid", "fine")
    }
    for level in ("coarse", "mid", "fine"):
        for sentinel in _ANNOTATION_SENTINELS:
            if sentinel not in material[level]:
                material[level][sentinel] = {
                    "parent": sentinel if level != "coarse" else "",
                    "definition": f"Synthetic explicit {sentinel} mapping state for v7 interoperability.",
                    "allowed_downstream_use": "excluded",
                    "_synthetic": True,
                }

    lookup: dict[tuple[str, str], str] = {}
    for level in ("coarse", "mid", "fine"):
        for label in material[level]:
            node_kind = "annotation_status" if label in _ANNOTATION_SENTINELS else (
                "state" if level == "fine" else "identity"
            )
            lookup[(level, label)] = _ontology_id(level, label, node_kind)

    records: list[dict[str, object]] = []
    for level in ("coarse", "mid", "fine"):
        for label, details in material[level].items():
            node_kind = "annotation_status" if label in _ANNOTATION_SENTINELS else (
                "state" if level == "fine" else "identity"
            )
            parent_label = _text(details.get("parent"))
            parent_id = ""
            if parent_label:
                candidate_levels = ("coarse",) if level == "mid" else (("mid", "coarse") if level == "fine" else ())
                for candidate_level in candidate_levels:
                    if (candidate_level, parent_label) in lookup:
                        parent_id = lookup[(candidate_level, parent_label)]
                        break
                if not parent_id:
                    raise Stage2Error(
                        f"STAGE2_BLOCKED_ONTOLOGY_PARENT: {level}:{label} -> {parent_label}"
                    )
            allowed = details.get("allowed_downstream_use", "")
            if isinstance(allowed, list):
                allowed = ";".join(map(str, allowed))
            records.append(
                {
                    "ontology_id": lookup[(level, label)],
                    "level": level,
                    "node_kind": node_kind,
                    "label": label,
                    "parent_label": parent_label,
                    "parent_ontology_id": parent_id,
                    "definition": _text(details.get("definition")),
                    "reliability_rule": _text(details.get("reliability_rule")),
                    "allowed_downstream_use": _text(allowed),
                    "is_exclusive": node_kind == "identity",
                    "is_synthetic": bool(details.get("_synthetic", False)),
                    "legacy_annotation_role": "annotation_radar" if label in hooks else "",
                    "source_version": source_version,
                    "source_sha256": source_hash,
                    "provenance_status": (
                        "v7_explicit_mapping_status"
                        if details.get("_synthetic", False)
                        else "v6_response_blind_ontology_inherited"
                    ),
                }
            )
    ontology = _stable(
        pd.DataFrame(records),
        ONTOLOGY_COLUMNS,
        ["level", "node_kind", "label"],
    )
    if ontology["ontology_id"].duplicated().any():
        raise Stage2Error("STAGE2_BLOCKED_ONTOLOGY: duplicate ontology id")
    return ontology, lookup


def _resolution_status(coarse: str, mid: str, fine: str) -> str:
    values = {coarse, mid, fine}
    if "Low_quality_or_ambient" in values:
        return "low_quality"
    if "Mixed" in values:
        return "mixed"
    if coarse == "Unknown" or mid == "Unknown":
        return "unknown"
    return "resolved"


def _canonical_label_mapping(
    path: Path,
    lookup: Mapping[tuple[str, str], str],
) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = [
        "cohort_id",
        "annotation_field",
        "original_label",
        "harmonized_coarse",
        "harmonized_mid",
        "harmonized_fine",
        "confidence",
    ]
    _require_columns(frame, required, "original label mapping")
    if frame.duplicated(["cohort_id", "annotation_field", "original_label"]).any():
        raise Stage2Error("STAGE2_BLOCKED_LABEL_MAPPING: duplicate canonical mapping key")
    source_hash = sha256(path)
    records: list[dict[str, object]] = []
    for row in frame.to_dict("records"):
        coarse = _text(row["harmonized_coarse"]) or "Unknown"
        mid = _text(row["harmonized_mid"]) or "Unknown"
        fine = _text(row["harmonized_fine"]) or "Unknown"
        targets = {
            "coarse": lookup.get(("coarse", coarse), ""),
            "mid": lookup.get(("mid", mid), ""),
            "fine": lookup.get(("fine", fine), ""),
        }
        missing = [level for level, target in targets.items() if not target]
        if missing:
            raise Stage2Error(
                "STAGE2_BLOCKED_LABEL_MAPPING_TARGET: "
                f"{_text(row['cohort_id'])}/{_text(row['original_label'])} missing {','.join(missing)}"
            )
        key = [row[column] for column in required]
        records.append(
            {
                "mapping_id": "canonical__" + _mapping_id(key),
                "cohort_id": _text(row["cohort_id"]),
                "annotation_field": _text(row["annotation_field"]),
                "original_label": _text(row["original_label"]),
                "harmonized_coarse": coarse,
                "harmonized_mid": mid,
                "harmonized_fine": fine,
                "coarse_ontology_id": targets["coarse"],
                "mid_ontology_id": targets["mid"],
                "fine_ontology_id": targets["fine"],
                "mapping_confidence": _text(row["confidence"]) or "unknown",
                "mapping_status": "canonical",
                "mapping_resolution_status": _resolution_status(coarse, mid, fine),
                "shadow_candidate_level": "",
                "shadow_candidate_label": "",
                "shadow_candidate_ontology_id": "",
                "n_shadow_observations": 0,
                "mapping_source": "v6_original_to_harmonized_label_map",
                "source_sha256": source_hash,
            }
        )
    return pd.DataFrame(records)


def _shadow_level_and_column(columns: Iterable[str]) -> tuple[str, str]:
    candidates: list[tuple[str, str]] = []
    for column in columns:
        lowered = str(column).lower()
        if not lowered.startswith("marker_based_"):
            continue
        for level in ("coarse", "mid", "fine"):
            if level in lowered:
                candidates.append((level, str(column)))
                break
    if len(candidates) != 1:
        raise Stage2Error(
            "STAGE2_BLOCKED_SHADOW_MAPPING: expected one marker_based coarse/mid/fine column"
        )
    return candidates[0]


def _shadow_label_mapping(
    paths: Sequence[Path],
    lookup: Mapping[tuple[str, str], str],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for path in sorted(paths, key=lambda item: str(item)):
        header = pd.read_csv(path, nrows=0)
        _require_columns(header, ["cohort_id"], f"shadow mapping {path.name}")
        level, candidate_column = _shadow_level_and_column(header.columns)
        counts: Counter[tuple[str, str]] = Counter()
        for chunk in pd.read_csv(
            path,
            usecols=["cohort_id", candidate_column],
            chunksize=250_000,
        ):
            for row in chunk.itertuples(index=False, name=None):
                cohort, candidate = _text(row[0]), _text(row[1])
                if cohort and candidate:
                    counts[(cohort, candidate)] += 1
        source_hash = sha256(path)
        for (cohort, candidate), count in sorted(counts.items()):
            records.append(
                {
                    "mapping_id": "shadow__"
                    + _mapping_id([path.name, cohort, level, candidate]),
                    "cohort_id": cohort,
                    "annotation_field": f"marker_based_{level}",
                    "original_label": candidate,
                    "harmonized_coarse": "Unknown",
                    "harmonized_mid": "Unknown",
                    "harmonized_fine": "Unknown",
                    "coarse_ontology_id": lookup[("coarse", "Unknown")],
                    "mid_ontology_id": lookup[("mid", "Unknown")],
                    "fine_ontology_id": lookup[("fine", "Unknown")],
                    "mapping_confidence": "candidate_only",
                    "mapping_status": "shadow_candidate",
                    "mapping_resolution_status": "shadow_candidate",
                    "shadow_candidate_level": level,
                    "shadow_candidate_label": candidate,
                    "shadow_candidate_ontology_id": lookup.get((level, candidate), ""),
                    "n_shadow_observations": int(count),
                    "mapping_source": path.name,
                    "source_sha256": source_hash,
                }
            )
    return pd.DataFrame(records)


def compile_cell_state_ontology(
    ontology_path: str | Path,
    original_mapping_path: str | Path,
    shadow_mapping_paths: Sequence[str | Path] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compile identity/state ontology and canonical plus optional shadow mappings."""

    ontology, lookup = _compile_ontology_table(Path(ontology_path))
    canonical = _canonical_label_mapping(Path(original_mapping_path), lookup)
    shadow = _shadow_label_mapping([Path(path) for path in shadow_mapping_paths], lookup)
    mappings = pd.concat([canonical, shadow], ignore_index=True, sort=False)
    mappings = _stable(
        mappings,
        LABEL_MAPPING_COLUMNS,
        ["mapping_status", "cohort_id", "annotation_field", "original_label", "mapping_id"],
    )
    if mappings["mapping_id"].duplicated().any():
        raise Stage2Error("STAGE2_BLOCKED_LABEL_MAPPING: duplicate mapping id")
    return ontology, mappings


def compile_stage2_vocabulary(
    *,
    legacy_membership_path: str | Path,
    legacy_dictionary_path: str | Path,
    signature_registry_path: str | Path,
    hallmark_gmt_path: str | Path,
    reactome_gmt_path: str | Path,
    sentinel_config_path: str | Path,
    ontology_path: str | Path,
    original_mapping_path: str | Path,
    shadow_mapping_paths: Sequence[str | Path] = (),
    expected_legacy_membership_sha256: str = FROZEN_FM_MEMBERSHIP_SHA256,
) -> Stage2Vocabulary:
    """Compile the complete in-memory Stage 2 vocabulary bundle."""

    legacy_dictionary, legacy_membership = compile_legacy_modules(
        legacy_membership_path,
        legacy_dictionary_path,
        expected_sha256=expected_legacy_membership_sha256,
    )
    mechanism_dictionary, mechanism_membership = compile_mechanism_axes(
        signature_registry_path,
        hallmark_gmt_path,
        reactome_gmt_path,
        sentinel_config_path,
    )
    ontology, label_mapping = compile_cell_state_ontology(
        ontology_path,
        original_mapping_path,
        shadow_mapping_paths=shadow_mapping_paths,
    )
    module_dictionary = _stable(
        pd.concat([legacy_dictionary, mechanism_dictionary], ignore_index=True),
        MODULE_DICTIONARY_COLUMNS,
        ["feature_family", "parent_axis_id", "feature_id"],
    )
    module_membership = _stable(
        pd.concat([legacy_membership, mechanism_membership], ignore_index=True),
        MODULE_MEMBERSHIP_COLUMNS,
        ["feature_family", "parent_axis_id", "feature_id", "gene_symbol"],
    )
    if module_dictionary["feature_id"].duplicated().any():
        raise Stage2Error("STAGE2_BLOCKED_VOCABULARY: duplicate feature id")
    return Stage2Vocabulary(
        module_dictionary=module_dictionary,
        module_membership=module_membership,
        cell_state_ontology=ontology,
        label_mapping=label_mapping,
    )


__all__ = [
    "FROZEN_FM_MEMBERSHIP_SHA256",
    "LABEL_MAPPING_COLUMNS",
    "MODULE_DICTIONARY_COLUMNS",
    "MODULE_MEMBERSHIP_COLUMNS",
    "ONTOLOGY_COLUMNS",
    "Stage2Vocabulary",
    "compile_cell_state_ontology",
    "compile_legacy_modules",
    "compile_mechanism_axes",
    "compile_stage2_vocabulary",
]
