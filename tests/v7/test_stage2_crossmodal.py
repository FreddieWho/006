from __future__ import annotations

import numpy as np
import pandas as pd

from src.v7.ontology.diagnostics import (
    cross_modal_mapping,
    gse193736_repeatability_diagnostics,
    gse207422_pair_vector_diagnostics,
    gse207422_paired_descriptive_diagnostics,
)


def test_gse207422_exact_pairs_are_primary_and_p07_is_support_only() -> None:
    rows = []
    for feature_id in ("F1", "F2"):
        for patient_id, match_status, bridge_role, sc, bulk in (
            ("P05", "exact_matched", "paired", 1.0, 1.5),
            ("P08", "exact_matched", "paired", 2.0, 2.5),
            ("P07", "timepoint_mismatch", "support_only", 9.0, -9.0),
        ):
            rows.extend(
                [
                    {
                        "feature_id": feature_id,
                        "patient_id": patient_id,
                        "match_status": match_status,
                        "bridge_role": bridge_role,
                        "modality": "scRNA",
                        "score_native": sc,
                    },
                    {
                        "feature_id": feature_id,
                        "patient_id": patient_id,
                        "match_status": match_status,
                        "bridge_role": bridge_role,
                        "modality": "bulkRNA",
                        "score_native": bulk,
                    },
                ]
            )

    result = gse207422_paired_descriptive_diagnostics(pd.DataFrame(rows))

    assert result.n_exact_pairs.eq(2).all()
    assert result.n_support_only_pairs.eq(1).all()
    assert result.primary_pair_ids.eq("P05|P08").all()
    assert result.support_only_pair_ids.eq("P07").all()
    assert result.paired_diagnostic_status.eq("paired_descriptive_insufficient_n").all()
    assert result.paired_spearman.isna().all()

    route = gse207422_pair_vector_diagnostics(pd.DataFrame(rows))
    assert set(route.comparison_id) == {"P05", "P07", "P08", "P08-minus-P05"}
    assert (
        route.set_index("comparison_id").loc["P07", "comparison_status"]
        == "SUPPORT_ONLY"
    )
    assert (
        route.set_index("comparison_id").loc[
            "P08-minus-P05", "comparison_status"
        ]
        == "PAIRED_DESCRIPTIVE_INSUFFICIENT_PATIENT_N"
    )


def test_gse193736_repeatability_uses_replicate_condition_profiles_only() -> None:
    rows = []
    design_rows = []
    for lineage in ("CD4", "CD8"):
        for perturbation in ("LTBR", "tNGFR"):
            for rest_stim in ("rest", "stim"):
                base = float(len(rows) + 1)
                for replicate in (1, 2, 3):
                    sample = f"{lineage}_{perturbation}_{rest_stim}_{replicate}"
                    design_rows.append(
                        {
                            "sample_id": sample,
                            "lineage": lineage,
                            "perturbation": perturbation,
                            "rest_stim": rest_stim,
                            "replicate": replicate,
                        }
                    )
                    rows.append(
                        {
                            "feature_id": "F1",
                            "sample_key": sample,
                            "score_native": base + replicate * 0.01,
                        }
                    )

    result = gse193736_repeatability_diagnostics(
        pd.DataFrame(rows), pd.DataFrame(design_rows)
    ).iloc[0]

    assert result["n_conditions"] == 8
    assert result["n_replicate_pairs"] == 3
    assert np.isclose(result["replicate_profile_median_spearman"], 1.0)
    assert result["repeatability_status"] == "REPEATABILITY_COMPUTED_NO_DIRECTIONAL_CLAIM"


def test_cross_modal_mapping_never_copies_scrna_measurable_to_other_modalities() -> None:
    scores = pd.DataFrame(
        {
            "feature_family": ["mechanism_component"] * 3,
            "feature_id": ["F"] * 3,
            "modality": ["scRNA", "bulkRNA", "bulk_perturbation"],
            "aggregation_level": ["expression_unit", "bulk_sample", "perturbation_sample"],
            "score_native": [1.0, 2.0, 3.0],
            "cohort_id": ["SC", "BULK", "PERT"],
            "patient_key": ["P", "P", ""],
            "coverage": [1.0, 1.0, 1.0],
        }
    )
    reliability = pd.DataFrame(
        {
            "feature_id": ["F"],
            "technical_status": ["valid"],
            "D2_class": ["measurable"],
            "D2_reason": ["scRNA_stability"],
        }
    )
    paired = pd.DataFrame(
        {
            "feature_id": ["F"],
            "n_exact_pairs": [2],
            "n_support_only_pairs": [1],
            "paired_spearman": [np.nan],
            "paired_diagnostic_status": ["paired_descriptive_insufficient_n"],
        }
    )
    repeatability = pd.DataFrame(
        {
            "feature_id": ["F"],
            "replicate_profile_median_spearman": [0.8],
            "n_replicate_pairs": [3],
            "repeatability_status": ["REPEATABILITY_COMPUTED_NO_DIRECTIONAL_CLAIM"],
        }
    )

    mapped = cross_modal_mapping(
        scores,
        reliability,
        paired_bulk=paired,
        perturb_repeatability=repeatability,
    ).set_index("modality")

    assert mapped.loc["scRNA", "D2_class"] == "measurable"
    assert mapped.loc["bulkRNA", "D2_class"] == "paired_descriptive_insufficient_n"
    assert mapped.loc["bulkRNA", "n_paired_bulk_evidence"] == 2
    assert (
        mapped.loc["bulk_perturbation", "D2_class"]
        == "interface_only_no_directional_validation"
    )
    assert mapped.loc["bulk_perturbation", "perturbation_repeatability"] == 0.8
