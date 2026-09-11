"""Platform adapters returning the shared :class:`SpatialUnitData` contract."""

from .h5ad_counts import load_h5ad_counts
from .visium_10x import load_visium_10x
from .xenium import load_xenium_h5ad

__all__ = ["load_h5ad_counts", "load_visium_10x", "load_xenium_h5ad"]
