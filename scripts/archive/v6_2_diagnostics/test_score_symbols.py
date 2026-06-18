import sys
import time
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

for cid, layer in [('GSE225063', None), ('GSE272993', None), ('GSE273718', None), ('krishna_2021_rcc', None)]:
    path = f'data/processed/srt/raw/{cid.lower()}.h5ad'
    print(f'\n=== {cid} ===', flush=True)
    ad = sc.read_h5ad(path, backed='r')
    syms = r.resolve_var_symbols(ad)
    overlap = len(set(syms) & r.ALL_MARKER_GENES)
    print('resolved overlap', overlap, 'var head', ad.var_names[:3].tolist(), flush=True)
    start = time.time()
    scores = r.compute_module_scores_for_object(ad, layer, cid)
    print('scores shape', scores.shape, 'time', time.time()-start, flush=True)
    print(scores.iloc[:2, :5], flush=True)
    ad.file.close()
