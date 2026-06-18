import sys
sys.path.insert(0,'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc, pandas as pd, numpy as np
for cid in ['GSE272734','GSE225063','GSE272993','GSE288199','GSE169246','TASK02']:
    print('===',cid, flush=True)
    try:
        ad=sc.read_h5ad(f'data/processed/srt/raw/{cid.lower()}.h5ad', backed='r')
    except Exception as e:
        print('open error', e); continue
    layer = r.parse_layer('counts') if 'counts' in ad.layers else (r.parse_layer('data') if 'data' in ad.layers else None)
    syms=r.resolve_var_symbols(ad)
    print(' overlap', len(set(syms)&r.ALL_MARKER_GENES), 'layer', layer, flush=True)
    scores = r.compute_module_scores_for_object(ad, layer, cid)
    print(' max zscore coarse', scores[['T_NK','B_Plasma','Myeloid','DC_APC','CAF_Stromal','Endothelial','Tumor_like','Cycling','Low_quality_or_ambient']].max().sort_values(ascending=False).head(5).to_dict(), flush=True)
    labels = r.assign_labels_from_scores(scores, cid)
    print(labels['marker_based_coarse_label'].value_counts().to_dict(), flush=True)
    ad.file.close()
