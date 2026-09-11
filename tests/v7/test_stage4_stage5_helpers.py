import numpy as np
from scipy import sparse

from v7.clinical_baseline import _pre_flag, _match_dataset
from v7.spatial_discovery import _marker_proxy, _nmf_summary
from v7.spatial_stats.graph import build_section_graph


def test_pre_timepoint_classification_is_conservative():
    assert _pre_flag("T0") == "primary_pre"
    assert _pre_flag("baseline") == "primary_pre"
    assert _pre_flag("on_treatment|pre") == "ambiguous_pre"
    assert _pre_flag("post_treatment") == "not_pre"


def test_dataset_aliases_are_case_normalized():
    assert _match_dataset("lambrecht_hcc") == "LAMBRECHT_HCC"
    assert _match_dataset("GSE301741") == "GSE301741"
    assert _match_dataset("unregistered") is None


def test_marker_proxy_keeps_missing_markers_explicit():
    counts = sparse.csr_matrix(np.asarray([[2, 0], [0, 4]], dtype=int))
    proxies, audit = _marker_proxy(counts, ["CD3D", "LST1"])
    assert proxies.shape == (2, 6)
    assert audit.loc[audit.proxy_id.eq("T_NK_proxy"), "status"].item() == "proxy_available"
    assert audit.loc[audit.proxy_id.eq("T_NK_proxy"), "coverage"].item() > 0
    assert audit.loc[audit.proxy_id.eq("B_Plasma_proxy"), "status"].item() == "NOT_ESTIMABLE_NO_MARKERS"
    assert np.isnan(proxies["B_Plasma_proxy"].to_numpy()).all()


def test_nmf_offset_is_derived_before_subsampling():
    matrix = np.ones((30, 4), dtype=float)
    matrix[:, 1:] += np.linspace(0.0, 2.0, 30)[:, None]
    matrix[0, 0] = -100.0
    coordinates = np.column_stack([np.arange(30, dtype=float), np.zeros(30)])
    sections = np.repeat("section-1", 30)
    graph = build_section_graph(coordinates, sections, method="knn", k=2)
    result = _nmf_summary(
        matrix,
        graph,
        section_ids=sections,
        unit_id="offset-regression",
        seeds=(20260911,),
        ks=(2,),
        sample_limit=5,
    )
    assert result.loc[result.k.eq(2), "status"].eq("ESTIMABLE").all()
