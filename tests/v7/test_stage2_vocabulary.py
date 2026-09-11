from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.v7.ontology.contracts import Stage2Error
from src.v7.ontology.vocabulary import (
    FROZEN_FM_IDS,
    compile_cell_state_ontology,
    compile_legacy_modules,
    compile_mechanism_axes,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _legacy_inputs(tmp_path: Path) -> tuple[Path, Path, str]:
    membership = tmp_path / "membership.csv"
    dictionary = tmp_path / "dictionary.csv"
    pd.DataFrame(
        [
            {
                "frozen_module_id": feature_id,
                "primary_method": "lda_topic",
                "primary_run": "frozen_run",
                "primary_module_id": f"source_{feature_id}",
                "gene": f"GENE{index}",
                "membership_weight": float(index),
                "relative_weight": float(index) / 10,
                "is_top_gene": "yes",
            }
            for index, feature_id in enumerate(FROZEN_FM_IDS, start=1)
        ]
    ).to_csv(membership, index=False)
    pd.DataFrame(
        [
            {
                "frozen_module_id": feature_id,
                "primary_method": "lda_topic",
                "primary_run": "frozen_run",
                "primary_module_id": f"source_{feature_id}",
                "support_tier": "multi_method_supported",
                "n_supporting_methods": 2,
                "bootstrap_median_top40_jaccard": 0.7,
                "leakage_detected": False,
                "default_for_phase7": True,
            }
            for feature_id in FROZEN_FM_IDS
        ]
    ).to_csv(dictionary, index=False)
    return membership, dictionary, _sha256(membership)


def test_legacy_fm_is_hash_locked_and_long_form(tmp_path: Path) -> None:
    membership, dictionary, digest = _legacy_inputs(tmp_path)
    compiled_dictionary, compiled_membership = compile_legacy_modules(
        membership, dictionary, expected_sha256=digest
    )

    assert tuple(compiled_dictionary["feature_id"]) == FROZEN_FM_IDS
    assert set(compiled_dictionary["feature_family"]) == {"legacy_fm"}
    assert compiled_dictionary["scoreable"].all()
    assert compiled_dictionary["n_genes"].eq(1).all()
    assert set(compiled_membership["direction"]) == {1}
    assert compiled_membership["is_top_gene"].all()
    assert set(compiled_membership["source_sha256"]) == {digest}

    with pytest.raises(Stage2Error, match="STAGE2_BLOCKED_FROZEN_FM_HASH"):
        compile_legacy_modules(membership, dictionary, expected_sha256="0" * 64)


def _mechanism_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    registry = tmp_path / "signatures.csv"
    hallmark = tmp_path / "hallmark.gmt"
    reactome = tmp_path / "reactome.gmt"
    sentinel = tmp_path / "sentinel.yaml"

    curated_names = [
        "angiogenesis",
        "stromal_TGFb",
        "hypoxia",
        "myeloid_suppression",
        "APC_DC_activation",
        "IFN_response",
        "antigen_presentation",
        "cytotoxicity",
        "NK_cytotoxicity",
        "exhaustion",
    ]
    rows = []
    for index, name in enumerate(curated_names):
        rows.append(
            {
                "feature_name": name,
                "gene_symbol": f"CURATED{index}",
                "direction": -1 if name == "IFN_response" else 1,
                "source_resource": "curated_step2_9_core",
            }
        )
    rows.extend(
        [
            {
                "feature_name": "antigen_presentation",
                "gene_symbol": "MIXED_CURATED",
                "direction": 1,
                "source_resource": "/mvp/program_signatures.json;curated_step2_9_core",
            },
            {
                "feature_name": "hcc_responder_like",
                "gene_symbol": "MVP_ONLY",
                "direction": 1,
                "source_resource": "/mvp/program_signatures.json",
            },
        ]
    )
    pd.DataFrame(rows).to_csv(registry, index=False)

    hallmark_sets = [
        "HALLMARK_ALLOGRAFT_REJECTION",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_HYPOXIA",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    ]
    hallmark.write_text(
        "".join(f"{name}\tdescription\t{name}_GENE\n" for name in hallmark_sets),
        encoding="utf-8",
    )
    reactome_sets = [
        "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
        "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
        "REACTOME_INTERFERON_GAMMA_SIGNALING",
        "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
        "REACTOME_NEUTROPHIL_DEGRANULATION",
        "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR",
        "REACTOME_SIGNALING_BY_VEGF",
        "REACTOME_SIGNALING_BY_WNT",
    ]
    reactome.write_text(
        "".join(f"{name}\tdescription\t{name}_GENE\n" for name in reactome_sets),
        encoding="utf-8",
    )
    yaml.safe_dump(
        {
            "sentinel_programs": {
                "TLS_B": {"genes": ["TLS1", "TLS2"]},
                "NEUTROPHIL_NET": {"genes": ["NET1", "NET2"]},
                "TUMOR_WNT_EXCLUSION": {"genes": ["WNT1", "WNT2"]},
                "CD8_PROGENITOR_EXHAUSTION": {"genes": ["PROG1", "PROG2"]},
                "CD8_TERMINAL_EXHAUSTION": {"genes": ["TERM1", "TERM2"]},
            }
        },
        sentinel.open("w", encoding="utf-8"),
        sort_keys=False,
    )
    return registry, hallmark, reactome, sentinel


def test_ten_mechanism_axes_are_explicit_and_response_blind(tmp_path: Path) -> None:
    inputs = _mechanism_inputs(tmp_path)
    dictionary, membership = compile_mechanism_axes(*inputs)

    parents = dictionary[dictionary["feature_level"] == "parent_axis"]
    components = dictionary[dictionary["feature_level"] == "component"]
    assert len(parents) == 10
    assert len(components) == 31
    assert not parents["scoreable"].any()
    assert components["scoreable"].all()
    assert set(components["parent_axis_id"]) == set(parents["feature_id"])
    assert "MVP_ONLY" not in set(membership["gene_symbol"])
    assert "MIXED_CURATED" in set(membership["gene_symbol"])
    assert "/mvp/" not in " ".join(membership["source_resource"].astype(str))

    unsigned = membership[membership["source_resource"].str.startswith("MSigDB_")]
    assert set(unsigned["direction"]) == {1}
    ifn = membership[
        membership["feature_id"].eq("mechanism__apc_ifn_apm__curated_ifn_response")
    ]
    assert set(ifn["direction"]) == {-1}
    pd.testing.assert_frame_equal(dictionary, compile_mechanism_axes(*inputs)[0])
    pd.testing.assert_frame_equal(membership, compile_mechanism_axes(*inputs)[1])


def _ontology_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    ontology = tmp_path / "ontology.yaml"
    mapping = tmp_path / "mapping.csv"
    shadow = tmp_path / "marker_based_mid_annotation.csv"
    ontology.write_text(
        yaml.safe_dump(
            {
                "ontology_version": "toy-v1",
                "levels": {
                    "coarse": {
                        "T_NK": {"definition": "T and NK"},
                        "Unknown": {"definition": "unknown", "allowed_downstream_use": "excluded"},
                        "Low_quality_or_ambient": {
                            "definition": "low quality",
                            "allowed_downstream_use": "excluded",
                        },
                    },
                    "mid": {
                        "CD8_T": {"parent": "T_NK", "definition": "CD8 T"},
                        "Unknown": {"parent": "Unknown", "definition": "unknown"},
                    },
                    "fine": {
                        "CD8_exhausted_like": {
                            "parent": "CD8_T",
                            "definition": "exhausted-like state",
                        }
                    },
                },
                "b2_annotation_hooks": ["CD8_exhausted_like"],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "cohort_id": "COHORT_A",
                "annotation_field": "cell_type",
                "original_label": "old_cd8",
                "harmonized_coarse": "T_NK",
                "harmonized_mid": "CD8_T",
                "harmonized_fine": "Unknown",
                "confidence": "high",
            }
        ]
    ).to_csv(mapping, index=False)
    pd.DataFrame(
        [
            {
                "cell_key": "a",
                "cohort_id": "COHORT_A",
                "marker_based_mid_label": "CD8_T",
            },
            {
                "cell_key": "b",
                "cohort_id": "COHORT_A",
                "marker_based_mid_label": "CD8_T",
            },
        ]
    ).to_csv(shadow, index=False)
    return ontology, mapping, shadow


def test_identity_state_separation_and_shadow_candidates(tmp_path: Path) -> None:
    ontology_path, mapping_path, shadow_path = _ontology_inputs(tmp_path)
    ontology, mappings = compile_cell_state_ontology(
        ontology_path, mapping_path, shadow_mapping_paths=[shadow_path]
    )

    assert set(ontology.loc[ontology["level"].isin(["coarse", "mid"]), "node_kind"]) == {
        "identity",
        "annotation_status",
    }
    fine_state = ontology[ontology["label"].eq("CD8_exhausted_like")].iloc[0]
    assert fine_state["node_kind"] == "state"
    assert fine_state["is_exclusive"] == False  # noqa: E712
    assert fine_state["legacy_annotation_role"] == "annotation_radar"

    sentinels = ontology[ontology["node_kind"].eq("annotation_status")]
    assert set(sentinels["label"]) == {"Unknown", "Mixed", "Low_quality_or_ambient"}
    assert sentinels.groupby("label").size().eq(3).all()

    canonical = mappings[mappings["mapping_status"].eq("canonical")].iloc[0]
    assert canonical["harmonized_mid"] == "CD8_T"
    assert canonical["shadow_candidate_label"] == ""
    shadow = mappings[mappings["mapping_status"].eq("shadow_candidate")].iloc[0]
    assert shadow["shadow_candidate_level"] == "mid"
    assert shadow["shadow_candidate_label"] == "CD8_T"
    assert shadow["n_shadow_observations"] == 2
    assert shadow["harmonized_mid"] == "Unknown"
    assert shadow["shadow_candidate_ontology_id"].startswith("identity::mid::")
