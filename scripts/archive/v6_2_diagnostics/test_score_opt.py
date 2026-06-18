import sys
import time
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

path = 'data/processed/srt/raw/gse243013.h5ad'
print('opening', path, flush=True)
ad = sc.read_h5ad(path, backed='r')
print('opened', ad.shape, flush=True)
start = time.time()
scores = r.compute_module_scores_for_object(ad, None, 'GSE243013')
elapsed = time.time() - start
print('scores shape', scores.shape, 'time', elapsed, flush=True)
print(scores.iloc[:3, :5], flush=True)
ad.file.close()
