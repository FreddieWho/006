import sys, numpy as np, pandas as pd
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

cid='GSE115978'
ad=sc.read_h5ad(f'data/processed/srt/raw/{cid.lower()}.h5ad', backed='r')
syms=r.resolve_var_symbols(ad)
print('total genes', len(syms))
# some known markers
for g in ['CD3E','CD3D','CD79A','CD79B','LYZ','CD68','EPCAM','KRT5','COL1A1','PECAM1']:
    if g in syms:
        idx=syms.index(g)
        print(g, 'idx', idx)
        vals=ad.layers['counts'][:, idx].toarray().flatten()
        print(' nonzero', np.count_nonzero(vals), 'max', vals.max())
    else:
        print(g, 'not found')
# W for T_NK
mods = r.MODULES
print('T_NK genes', mods.get('T_NK', [])[:10])
# overlap with symbols
for mod, genes in list(mods.items())[:3]:
    present=[g for g in genes if g in syms]
    print(mod, 'present', len(present), present[:10])
