"""Evaluate marginal Moran nulls jointly for features with the same observed mask."""
import numpy as np
from scipy import sparse
from .spatial_stats.graph import validate_section_local_graph


def matrix_moran(values, adjacency, sections, *, permutations=999, seed=20260911):
    values=np.asarray(values,float)
    if values.ndim!=2 or len(sections)!=len(values):
        raise ValueError("unaligned matrix and sections")
    graph=sparse.csr_matrix(adjacency,dtype=float).copy()
    graph.setdiag(0);graph.eliminate_zeros()
    validate_section_local_graph(graph,sections)
    masks={}
    for j in range(values.shape[1]):
        key=np.packbits(np.isfinite(values[:,j])).tobytes()
        masks.setdefault(key,[]).append(j)
    results=[None]*values.shape[1]
    for ordinal,columns in enumerate(masks.values()):
        valid=np.isfinite(values[:,columns[0]])
        index=np.flatnonzero(valid)
        if len(index)<3:
            for j in columns:
                results[j]=dict(status="NOT_ESTIMABLE_TOO_FEW_OBSERVATIONS",observed=np.nan,p_value=np.nan,n_permutations=0)
            continue
        local=graph[index][:,index]
        weight=float(local.sum())
        if weight<=0:
            for j in columns:
                results[j]=dict(status="NOT_ESTIMABLE_NO_EDGES",observed=np.nan,p_value=np.nan,n_permutations=0)
            continue
        section=np.asarray(sections)[index]
        groups=[np.flatnonzero(section==s) for s in np.unique(section)]
        centered=values[np.ix_(index,columns)].copy()
        for group in groups:
            centered[group]-=centered[group].mean(axis=0)
        denominator=(centered**2).sum(axis=0)
        usable=denominator>0
        scale=np.divide(len(index)/weight,denominator,out=np.zeros_like(denominator),where=usable)
        observed=(centered*(local@centered)).sum(axis=0)*scale
        rng=np.random.default_rng(seed+ordinal)
        null=np.empty((permutations,len(columns)))
        for draw in range(permutations):
            order=np.arange(len(index))
            for group in groups:
                order[group]=rng.permutation(group)
            permuted=centered[order]
            null[draw]=(permuted*(local@permuted)).sum(axis=0)*scale
        for pos,j in enumerate(columns):
            if not usable[pos]:
                results[j]=dict(status="NOT_ESTIMABLE_ZERO_VARIANCE",observed=np.nan,p_value=np.nan,n_permutations=0)
            else:
                center=np.median(null[:,pos])
                exceed=np.sum(np.abs(null[:,pos]-center)>=abs(observed[pos]-center))
                results[j]=dict(status="ESTIMABLE",observed=float(observed[pos]),
                                p_value=float((1+exceed)/(1+permutations)),n_permutations=permutations)
    return results
