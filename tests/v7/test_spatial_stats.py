from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.v7.spatial_stats import (
    build_section_graph,
    edge_association,
    moran_permutation_test,
    moran_statistic,
    patient_nested_summary,
    permute_within_sections,
    sparse_variogram,
    validate_section_local_graph,
)


def test_knn_graph_is_sparse_symmetric_and_never_crosses_sections() -> None:
    coordinates = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 0.0], [1.0, 0.0]]
    )
    sections = np.array(["A", "A", "A", "B", "B"])

    graph = build_section_graph(coordinates, sections, method="knn", k=1)

    assert sparse.isspmatrix_csr(graph)
    assert (graph != graph.T).nnz == 0
    assert graph.diagonal().sum() == 0
    assert validate_section_local_graph(graph, sections)["status"] == "PASS"
    rows, columns = graph.nonzero()
    assert np.all(sections[rows] == sections[columns])
    assert graph.nnz <= 2 * len(coordinates)


def test_radius_graph_uses_only_pairs_inside_radius_and_section() -> None:
    coordinates = np.array([[0.0, 0.0], [0.4, 0.0], [2.0, 0.0], [0.1, 0.0]])
    sections = np.array(["A", "A", "A", "B"])

    graph = build_section_graph(
        coordinates, sections, method="radius", radius=0.5
    )

    assert graph.nnz == 2
    assert graph[0, 1] == graph[1, 0] == 1.0
    assert graph[0, 3] == 0.0


def test_section_validator_rejects_cross_section_edges() -> None:
    graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))

    with pytest.raises(ValueError, match="cross-section"):
        validate_section_local_graph(graph, np.array(["A", "B"]))


def test_moran_statistic_detects_clustered_signal() -> None:
    coordinates = np.column_stack([np.arange(8, dtype=float), np.zeros(8)])
    sections = np.array(["S"] * 8)
    graph = build_section_graph(coordinates, sections, method="knn", k=2)
    values = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 5.0, 5.0, 5.0])

    result = moran_statistic(values, graph, section_ids=sections)

    assert result["status"] == "ESTIMABLE"
    assert result["estimate"] > 0
    assert result["n_observations"] == 8


def test_moran_with_sections_does_not_turn_between_section_offset_into_signal() -> None:
    coordinates = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]
    )
    sections = np.array(["A"] * 3 + ["B"] * 3)
    graph = build_section_graph(coordinates, sections, method="knn", k=1)

    result = moran_statistic(
        np.array([0.0, 0.0, 0.0, 100.0, 100.0, 100.0]),
        graph,
        section_ids=sections,
    )

    assert result["status"] == "NOT_ESTIMABLE_ZERO_VARIANCE"


def test_sparse_variogram_uses_graph_edges_and_increases_for_linear_field() -> None:
    coordinates = np.column_stack([np.arange(7, dtype=float), np.zeros(7)])
    sections = np.array(["S"] * 7)
    graph = build_section_graph(
        coordinates, sections, method="radius", radius=3.1
    )

    result = sparse_variogram(
        np.arange(7, dtype=float),
        graph,
        coordinates,
        section_ids=sections,
        distance_bins=np.array([0.0, 1.5, 2.5, 3.5]),
    )
    estimated = result.loc[result.status.eq("ESTIMABLE")]

    assert len(estimated) == 3
    assert estimated.semivariance.is_monotonic_increasing
    assert estimated.n_edges.sum() == sparse.triu(graph, k=1).nnz


def test_edge_association_requires_composition_and_uses_sparse_edges() -> None:
    coordinates = np.column_stack([np.arange(6, dtype=float), np.zeros(6)])
    sections = np.array(["S"] * 6)
    graph = build_section_graph(coordinates, sections, method="knn", k=1)
    composition = pd.DataFrame(
        {
            "T": [0.0, 0.1, 0.2, 0.8, 0.9, 1.0],
            "Myeloid": [0.0, 0.2, 0.3, 0.7, 0.8, 1.0],
        }
    )

    missing = edge_association(graph, None, source="T", target="Myeloid")
    observed = edge_association(
        graph,
        composition,
        source="T",
        target="Myeloid",
        section_ids=sections,
    )

    assert missing["status"] == "NOT_ESTIMABLE_MISSING_COMPOSITION"
    assert np.isnan(missing["estimate"])
    assert observed["status"] == "ESTIMABLE"
    assert observed["estimate"] > 0


def test_spatial_permutation_is_within_section_and_reproducible() -> None:
    values = np.array([1.0, 2.0, 3.0, 10.0, 20.0, np.nan])
    sections = np.array(["A", "A", "A", "B", "B", "B"])
    first = permute_within_sections(values, sections, np.random.default_rng(4))
    assert sorted(first[:3]) == [1.0, 2.0, 3.0]
    assert sorted(first[3:5]) == [10.0, 20.0]
    assert np.isnan(first[5])

    coordinates = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]
    )
    graph = build_section_graph(coordinates, sections, method="knn", k=1)
    one = moran_permutation_test(values, graph, sections, permutations=29, seed=7)
    two = moran_permutation_test(values, graph, sections, permutations=29, seed=7)
    assert one["status"] == "ESTIMABLE"
    assert one["p_value"] == two["p_value"]
    assert np.array_equal(one["permuted_statistics"], two["permuted_statistics"])


def test_patient_nested_summary_weights_patients_not_sections() -> None:
    sections = pd.DataFrame(
        {
            "cohort_id": ["C"] * 4,
            "patient_id": ["P1", "P1", "P1", "P2"],
            "section_id": ["S1", "S2", "S3", "S4"],
            "feature": ["F"] * 4,
            "estimate": [0.0, 0.0, 0.0, 10.0],
            "status": ["ESTIMABLE"] * 4,
        }
    )

    patient, cohort = patient_nested_summary(
        sections, group_columns=("feature",)
    )

    assert patient.set_index("patient_id").loc["P1", "estimate"] == 0.0
    assert patient.set_index("patient_id").loc["P2", "estimate"] == 10.0
    assert cohort.loc[0, "estimate"] == 5.0
    assert cohort.loc[0, "n_patients"] == 2


def test_patient_nested_summary_keeps_not_estimable_explicit() -> None:
    sections = pd.DataFrame(
        {
            "cohort_id": ["C"],
            "patient_id": ["P"],
            "section_id": ["S"],
            "estimate": [np.nan],
            "status": ["NOT_ESTIMABLE_MISSING_COMPOSITION"],
        }
    )

    patient, cohort = patient_nested_summary(sections)

    assert patient.loc[0, "status"] == "NOT_ESTIMABLE_NO_VALID_SECTIONS"
    assert cohort.loc[0, "status"] == "NOT_ESTIMABLE_NO_VALID_PATIENTS"
