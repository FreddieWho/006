"""Sparse spatial statistics and section-preserving null permutations."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from .graph import validate_section_local_graph


def _graph_and_values(
    values: np.ndarray | Sequence[float], adjacency: sparse.spmatrix
) -> tuple[np.ndarray, sparse.csr_matrix]:
    array = np.asarray(values, dtype=float)
    graph = sparse.csr_matrix(adjacency, dtype=float)
    if array.ndim != 1 or graph.shape != (len(array), len(array)):
        raise ValueError("values and adjacency have incompatible shapes")
    if graph.data.size and (not np.isfinite(graph.data).all() or (graph.data < 0).any()):
        raise ValueError("adjacency weights must be finite and non-negative")
    graph = graph.copy()
    graph.setdiag(0)
    graph.eliminate_zeros()
    return array, graph


def moran_statistic(
    values: np.ndarray | Sequence[float],
    adjacency: sparse.spmatrix,
    *,
    section_ids: np.ndarray | Sequence[object] | None = None,
) -> dict[str, float | int | str]:
    """Compute a global Moran-type statistic over finite nodes and sparse edges."""

    array, graph = _graph_and_values(values, adjacency)
    if section_ids is not None:
        validate_section_local_graph(graph, section_ids)
    valid = np.isfinite(array)
    if valid.sum() < 3:
        return _not_estimable("NOT_ESTIMABLE_TOO_FEW_OBSERVATIONS", int(valid.sum()))
    retained = np.flatnonzero(valid)
    local_graph = graph[retained][:, retained]
    weight_sum = float(local_graph.sum())
    retained_values = array[retained]
    if section_ids is None:
        centered = retained_values - float(np.mean(retained_values))
    else:
        retained_sections = np.asarray(section_ids)[retained]
        centered = retained_values.copy()
        for section in np.unique(retained_sections):
            section_mask = retained_sections == section
            centered[section_mask] -= float(np.mean(centered[section_mask]))
    denominator = float(centered @ centered)
    if local_graph.nnz == 0 or weight_sum <= 0:
        return _not_estimable("NOT_ESTIMABLE_NO_EDGES", len(retained))
    if denominator <= 0:
        return _not_estimable("NOT_ESTIMABLE_ZERO_VARIANCE", len(retained))
    estimate = len(retained) / weight_sum * float(
        centered @ (local_graph @ centered)
    ) / denominator
    return {
        "status": "ESTIMABLE",
        "estimate": float(estimate),
        "n_observations": int(len(retained)),
        "n_directed_edges": int(local_graph.nnz),
        "weight_sum": weight_sum,
        "centering": "within_section" if section_ids is not None else "global",
    }


def _not_estimable(status: str, n: int = 0) -> dict[str, float | int | str]:
    return {
        "status": status,
        "estimate": float("nan"),
        "n_observations": int(n),
        "n_directed_edges": 0,
        "weight_sum": 0.0,
    }


def sparse_variogram(
    values: np.ndarray | Sequence[float],
    adjacency: sparse.spmatrix,
    coordinates: np.ndarray | Sequence[Sequence[float]],
    *,
    section_ids: np.ndarray | Sequence[object] | None = None,
    distance_bins: int | np.ndarray | Sequence[float] = 10,
) -> pd.DataFrame:
    """Estimate a semivariogram from graph edges only, never all spot pairs."""

    array, graph = _graph_and_values(values, adjacency)
    points = np.asarray(coordinates, dtype=float)
    if points.ndim != 2 or len(points) != len(array) or not np.isfinite(points).all():
        raise ValueError("coordinates must be a finite n-by-d array")
    if section_ids is not None:
        validate_section_local_graph(graph, section_ids)
    upper = sparse.triu(graph, k=1).tocoo()
    valid = np.isfinite(array[upper.row]) & np.isfinite(array[upper.col])
    rows = upper.row[valid]
    columns = upper.col[valid]
    if not len(rows):
        return pd.DataFrame(
            [
                {
                    "bin_left": np.nan,
                    "bin_right": np.nan,
                    "distance_mean": np.nan,
                    "semivariance": np.nan,
                    "n_edges": 0,
                    "status": "NOT_ESTIMABLE_NO_EDGES",
                }
            ]
        )
    distances = np.linalg.norm(points[rows] - points[columns], axis=1)
    semivariance = 0.5 * (array[rows] - array[columns]) ** 2
    if isinstance(distance_bins, (int, np.integer)):
        if int(distance_bins) < 1:
            raise ValueError("distance_bins must be positive")
        maximum = float(distances.max())
        if maximum == 0:
            edges = np.array([0.0, np.finfo(float).eps])
        else:
            edges = np.linspace(0.0, maximum, int(distance_bins) + 1)
    else:
        edges = np.asarray(distance_bins, dtype=float)
        if edges.ndim != 1 or len(edges) < 2 or not np.all(np.diff(edges) > 0):
            raise ValueError("distance_bins edges must be strictly increasing")
    assignments = np.searchsorted(edges, distances, side="right") - 1
    assignments[distances == edges[-1]] = len(edges) - 2
    output: list[dict[str, float | int | str]] = []
    for index in range(len(edges) - 1):
        selected = assignments == index
        n_edges = int(selected.sum())
        output.append(
            {
                "bin_left": float(edges[index]),
                "bin_right": float(edges[index + 1]),
                "distance_mean": float(np.mean(distances[selected])) if n_edges else np.nan,
                "semivariance": float(np.mean(semivariance[selected])) if n_edges else np.nan,
                "n_edges": n_edges,
                "status": "ESTIMABLE" if n_edges else "NOT_ESTIMABLE_EMPTY_BIN",
            }
        )
    return pd.DataFrame(output)


def edge_association(
    adjacency: sparse.spmatrix,
    composition: pd.DataFrame | np.ndarray | None,
    *,
    source: str | int,
    target: str | int,
    section_ids: np.ndarray | Sequence[object] | None = None,
) -> dict[str, float | int | str]:
    """Compute a standardized source-to-neighbor-target edge association."""

    if composition is None:
        return {
            "status": "NOT_ESTIMABLE_MISSING_COMPOSITION",
            "estimate": np.nan,
            "n_observations": 0,
            "n_directed_edges": 0,
        }
    if isinstance(composition, pd.DataFrame):
        if source not in composition.columns or target not in composition.columns:
            return {
                "status": "NOT_ESTIMABLE_MISSING_COMPOSITION_FEATURE",
                "estimate": np.nan,
                "n_observations": int(len(composition)),
                "n_directed_edges": 0,
            }
        source_values = pd.to_numeric(composition[source], errors="coerce").to_numpy(float)
        target_values = pd.to_numeric(composition[target], errors="coerce").to_numpy(float)
    else:
        matrix = np.asarray(composition, dtype=float)
        if matrix.ndim != 2 or not isinstance(source, (int, np.integer)) or not isinstance(target, (int, np.integer)):
            raise ValueError("array composition requires integer source and target columns")
        source_values = matrix[:, int(source)]
        target_values = matrix[:, int(target)]
    _, graph = _graph_and_values(source_values, adjacency)
    if len(target_values) != len(source_values):
        raise ValueError("composition and adjacency have incompatible shapes")
    if section_ids is not None:
        validate_section_local_graph(graph, section_ids)

    def standardized(values: np.ndarray) -> np.ndarray:
        output = np.full(len(values), np.nan, dtype=float)
        groups = (
            np.asarray(section_ids)
            if section_ids is not None
            else np.repeat("ALL", len(values))
        )
        for group in np.unique(groups):
            selected = (groups == group) & np.isfinite(values)
            if selected.sum() < 2:
                continue
            scale = float(np.std(values[selected]))
            if scale > 0:
                output[selected] = (values[selected] - float(np.mean(values[selected]))) / scale
        return output

    source_z = standardized(source_values)
    target_z = standardized(target_values)
    rows, columns = graph.nonzero()
    edge_valid = np.isfinite(source_z[rows]) & np.isfinite(target_z[columns])
    rows = rows[edge_valid]
    columns = columns[edge_valid]
    weights = graph[rows, columns].A1
    node_valid = np.isfinite(source_z) & np.isfinite(target_z)
    raw_node_valid = np.isfinite(source_values) & np.isfinite(target_values)
    if len(rows) == 0 or node_valid.sum() < 3:
        return {
            "status": (
                "NOT_ESTIMABLE_ZERO_WITHIN_SECTION_VARIANCE"
                if raw_node_valid.sum() >= 3 and node_valid.sum() == 0
                else "NOT_ESTIMABLE_NO_VALID_EDGES"
            ),
            "estimate": np.nan,
            "n_observations": int(node_valid.sum()),
            "n_directed_edges": int(len(rows)),
        }
    products = source_z[rows] * target_z[columns]
    return {
        "status": "ESTIMABLE",
        "estimate": float(np.average(products, weights=weights)),
        "n_observations": int(node_valid.sum()),
        "n_directed_edges": int(len(rows)),
    }


def permute_within_sections(
    values: np.ndarray | Sequence[float],
    section_ids: np.ndarray | Sequence[object],
    rng: np.random.Generator,
) -> np.ndarray:
    """Shuffle finite values within sections while leaving missing positions fixed."""

    array = np.asarray(values, dtype=float)
    sections = np.asarray(section_ids)
    if array.ndim != 1 or sections.ndim != 1 or len(array) != len(sections):
        raise ValueError("values and section_ids must be aligned one-dimensional arrays")
    output = array.copy()
    for section in np.unique(sections):
        indices = np.flatnonzero((sections == section) & np.isfinite(array))
        output[indices] = rng.permutation(array[indices])
    return output


def moran_permutation_test(
    values: np.ndarray | Sequence[float],
    adjacency: sparse.spmatrix,
    section_ids: np.ndarray | Sequence[object],
    *,
    permutations: int = 999,
    seed: int = 20260904,
    alternative: str = "two-sided",
) -> dict[str, object]:
    """Moran test using exchangeability restricted to each section."""

    if permutations < 1:
        raise ValueError("permutations must be positive")
    if alternative not in {"two-sided", "greater", "less"}:
        raise ValueError("alternative must be 'two-sided', 'greater', or 'less'")
    array, graph = _graph_and_values(values, adjacency)
    sections = np.asarray(section_ids)
    validate_section_local_graph(graph, sections)
    observed = moran_statistic(array, graph, section_ids=sections)
    if observed["status"] != "ESTIMABLE":
        return {
            "status": observed["status"],
            "observed": np.nan,
            "p_value": np.nan,
            "n_permutations": 0,
            "permuted_statistics": np.array([], dtype=float),
        }
    rng = np.random.default_rng(seed)
    # Reuse the valid-node graph and the within-section centering contract for
    # every draw.  Calling ``moran_statistic`` recursively would rebuild the
    # same sparse subgraph and validation 99--999 times per feature, which is
    # unnecessarily expensive for large Xenium captures.
    valid = np.isfinite(array)
    retained = np.flatnonzero(valid)
    local_graph = graph[retained][:, retained]
    retained_sections = sections[retained]
    weight_sum = float(local_graph.sum())
    def centered_within(values: np.ndarray) -> np.ndarray:
        retained_values = values[retained].astype(float, copy=True)
        for section in np.unique(retained_sections):
            selected = retained_sections == section
            retained_values[selected] -= float(np.mean(retained_values[selected]))
        return retained_values

    observed_centered = centered_within(array)
    denominator = float(observed_centered @ observed_centered)

    # Within-section permutations preserve the section means. Cache the
    # centering and group membership once instead of sorting string section
    # identifiers and recentering for every draw. Bound dense working memory.
    groups = [np.flatnonzero(retained_sections == s) for s in np.unique(retained_sections)]
    batches = []
    for start in range(0, permutations, 32):
        count = min(32, permutations-start)
        null_centered = np.empty((len(retained), count), dtype=float)
        for column in range(count):
            for indices in groups:
                null_centered[indices,column] = rng.permutation(observed_centered[indices])
        numerators = np.asarray(null_centered * (local_graph @ null_centered)).sum(axis=0)
        batches.append(len(retained) / weight_sum * numerators / denominator)
    null = np.concatenate(batches)
    null = null[np.isfinite(null)]
    if not len(null):
        return {
            "status": "NOT_ESTIMABLE_PERMUTATION_NULL",
            "observed": observed["estimate"],
            "p_value": np.nan,
            "n_permutations": 0,
            "permuted_statistics": null,
        }
    observed_value = float(observed["estimate"])
    if alternative == "greater":
        exceedances = int(np.count_nonzero(null >= observed_value))
    elif alternative == "less":
        exceedances = int(np.count_nonzero(null <= observed_value))
    else:
        center = float(np.median(null))
        exceedances = int(np.count_nonzero(np.abs(null - center) >= abs(observed_value - center)))
    return {
        "status": "ESTIMABLE",
        "observed": observed_value,
        "p_value": float((exceedances + 1) / (len(null) + 1)),
        "n_permutations": int(len(null)),
        "alternative": alternative,
        "permuted_statistics": null,
    }
