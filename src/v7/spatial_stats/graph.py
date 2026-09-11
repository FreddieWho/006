"""Sparse, section-local spatial graph construction."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree


def _validated_inputs(
    coordinates: np.ndarray | Sequence[Sequence[float]],
    section_ids: np.ndarray | Sequence[object],
) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(coordinates, dtype=float)
    sections = np.asarray(section_ids)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError("coordinates must be an n-by-d array with d >= 2")
    if sections.ndim != 1 or len(sections) != len(points):
        raise ValueError("section_ids must contain one value per coordinate row")
    if not np.isfinite(points).all():
        raise ValueError("coordinates contain non-finite values")
    if any(value is None or (isinstance(value, float) and np.isnan(value)) for value in sections):
        raise ValueError("section_ids contain missing values")
    return points, sections


def _edge_weights(distances: np.ndarray, mode: str) -> np.ndarray:
    if mode == "binary":
        return np.ones(len(distances), dtype=float)
    if mode == "inverse_distance":
        positive = distances[distances > 0]
        floor = float(positive.min()) * 1e-6 if positive.size else 1.0
        return 1.0 / np.maximum(distances, floor)
    raise ValueError("weight_mode must be 'binary' or 'inverse_distance'")


def build_section_graph(
    coordinates: np.ndarray | Sequence[Sequence[float]],
    section_ids: np.ndarray | Sequence[object],
    *,
    method: str = "knn",
    k: int = 6,
    radius: float | None = None,
    weight_mode: str = "binary",
) -> sparse.csr_matrix:
    """Build a symmetric graph with a separate ``cKDTree`` for each section.

    The implementation never allocates an all-by-all distance matrix.  kNN
    edges are symmetrized by union; radius edges are emitted as undirected
    pairs.  Cross-section edges cannot be generated.
    """

    points, sections = _validated_inputs(coordinates, section_ids)
    if method not in {"knn", "radius"}:
        raise ValueError("method must be 'knn' or 'radius'")
    if method == "knn" and (not isinstance(k, (int, np.integer)) or k < 1):
        raise ValueError("k must be a positive integer")
    if method == "radius" and (radius is None or not np.isfinite(radius) or radius <= 0):
        raise ValueError("radius must be finite and positive for a radius graph")

    edge_rows: list[np.ndarray] = []
    edge_columns: list[np.ndarray] = []
    edge_distances: list[np.ndarray] = []
    for section in np.unique(sections):
        indices = np.flatnonzero(sections == section)
        if len(indices) < 2:
            continue
        local = points[indices]
        tree = cKDTree(local)
        if method == "knn":
            n_neighbors = min(int(k), len(indices) - 1)
            distances, neighbors = tree.query(local, k=n_neighbors + 1)
            distances = np.atleast_2d(distances)
            neighbors = np.atleast_2d(neighbors)
            if distances.shape[0] != len(indices):
                distances = distances.T
                neighbors = neighbors.T
            local_rows = np.repeat(np.arange(len(indices)), n_neighbors + 1)
            local_columns = neighbors.reshape(-1)
            local_distances = distances.reshape(-1)
            keep = (
                (local_rows != local_columns)
                & (local_columns < len(indices))
                & np.isfinite(local_distances)
            )
            rows = indices[local_rows[keep]]
            columns = indices[local_columns[keep]]
            distances_kept = local_distances[keep]
        else:
            pairs = tree.query_pairs(float(radius), output_type="ndarray")
            if pairs.size == 0:
                continue
            pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
            first = indices[pairs[:, 0]]
            second = indices[pairs[:, 1]]
            pair_distances = np.linalg.norm(
                points[first] - points[second], axis=1
            )
            rows = np.concatenate([first, second])
            columns = np.concatenate([second, first])
            distances_kept = np.concatenate([pair_distances, pair_distances])
        edge_rows.append(np.asarray(rows, dtype=np.int64))
        edge_columns.append(np.asarray(columns, dtype=np.int64))
        edge_distances.append(np.asarray(distances_kept, dtype=float))

    n = len(points)
    if not edge_rows:
        return sparse.csr_matrix((n, n), dtype=float)
    rows = np.concatenate(edge_rows)
    columns = np.concatenate(edge_columns)
    distances = np.concatenate(edge_distances)
    graph = sparse.coo_matrix(
        (_edge_weights(distances, weight_mode), (rows, columns)), shape=(n, n)
    ).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    graph = graph.maximum(graph.T).tocsr()
    graph.sort_indices()
    validate_section_local_graph(graph, sections)
    return graph


def validate_section_local_graph(
    adjacency: sparse.spmatrix,
    section_ids: np.ndarray | Sequence[object],
) -> dict[str, int | str]:
    """Raise on any cross-section edge and return a compact audit record."""

    graph = sparse.csr_matrix(adjacency)
    sections = np.asarray(section_ids)
    if graph.shape[0] != graph.shape[1] or graph.shape[0] != len(sections):
        raise ValueError("adjacency and section_ids have incompatible shapes")
    rows, columns = graph.nonzero()
    cross = sections[rows] != sections[columns]
    n_cross = int(np.count_nonzero(cross))
    if n_cross:
        raise ValueError(f"graph contains {n_cross} cross-section directed edges")
    return {
        "status": "PASS",
        "n_nodes": int(graph.shape[0]),
        "n_directed_edges": int(graph.nnz),
        "n_cross_section_edges": 0,
    }
