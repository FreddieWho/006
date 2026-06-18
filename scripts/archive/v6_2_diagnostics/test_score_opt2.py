import sys
import time
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

for cid, layer in [('GSE115978', 'counts'), ('GSE243013', None)]:
    path = f'data/processed/srt/raw/{cid.lower()}.h5ad'
    print(f'\n=== {cid} layer={layer} ===', flush=True)
    ad = sc.read_h5ad(path, backed='r')
    print('opened', ad.shape, 'layers', list(ad.layers.keys()), flush=True)
    start = time.time()
    scores = r.compute_module_scores_for_object(ad, layer, cid)
    elapsed = time.time() - start
    print('scores shape', scores.shape, 'time', elapsed, flush=True)
    print(scores.iloc[:2, :5], flush=True)
    ad.file.close()
