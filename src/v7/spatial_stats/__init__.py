"""Sparse, patient-aware spatial statistics for v7 Stage 3."""

from .aggregation import patient_nested_summary
from .graph import build_section_graph, validate_section_local_graph
from .resampling import patient_block_folds, patient_block_split
from .statistics import (
    edge_association,
    moran_permutation_test,
    moran_statistic,
    permute_within_sections,
    sparse_variogram,
)

__all__ = [
    "build_section_graph",
    "edge_association",
    "moran_permutation_test",
    "moran_statistic",
    "patient_block_folds",
    "patient_block_split",
    "patient_nested_summary",
    "permute_within_sections",
    "sparse_variogram",
    "validate_section_local_graph",
]
