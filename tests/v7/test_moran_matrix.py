import numpy as np
from scipy import sparse
from v7.moran_matrix import matrix_moran
from v7.spatial_stats.statistics import moran_permutation_test


def test_matrix_kernel_single_feature_matches_seeded_scalar_null():
    graph=sparse.diags([np.ones(5),np.ones(5)],[-1,1],shape=(6,6),format="csr")
    values=np.array([1.,2.,np.nan,9.,2.,6.])
    sections=np.array(["a"]*6)
    a=matrix_moran(values[:,None],graph,sections,permutations=99,seed=2)[0]
    b=moran_permutation_test(values,graph,sections,permutations=99,seed=2)
    np.testing.assert_allclose(a["observed"],b["observed"])
    assert a["p_value"]==b["p_value"]


def test_matrix_kernel_preserves_distinct_missingness_and_constant_status():
    graph=sparse.diags([np.ones(5),np.ones(5)],[-1,1],shape=(6,6),format="csr")
    x=np.array([[1.,np.nan,2.],[2.,1.,2.],[np.nan,3.,2.],[4.,2.,2.],[5.,4.,2.],[2.,6.,2.]])
    results=matrix_moran(x,graph,np.array(["a"]*6),permutations=19)
    assert [r["status"] for r in results]==["ESTIMABLE","ESTIMABLE","NOT_ESTIMABLE_ZERO_VARIANCE"]
