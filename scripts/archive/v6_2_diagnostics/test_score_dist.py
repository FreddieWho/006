import sys, numpy as np, pandas as pd
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

cid = 'GSE115978'
path = f'data/processed/srt/raw/{cid.lower()}.h5ad'
print('===', cid, flush=True)
ad = sc.read_h5ad(path, backed='r')
syms = r.resolve_var_symbols(ad)
print('overlap', len(set(syms)&r.ALL_MARKER_GENES), 'example symbols', syms[:10], flush=True)
print('layers', list(ad.layers.keys()), flush=True)
scores = r.compute_module_scores_for_object(ad, 'counts', cid)
print('scores shape', scores.shape, flush=True)
print('max per module:', scores.max().sort_values(ascending=False).head(10), flush=True)
print('nonzero cells per module:', (scores>0).sum().sort_values(ascending=False).head(10), flush=True)
labels = r.assign_labels_from_scores(scores, cid)
print(labels.value_counts().head(10), flush=True)
ad.file.close()
