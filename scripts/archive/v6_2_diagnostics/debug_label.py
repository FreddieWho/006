import sys, numpy as np, pandas as pd
sys.path.insert(0, 'results/v6_2/phase4a_cell_state_harmonization/scripts')
import run_phase4a as r
import scanpy as sc

cid='GSE115978'
ad=sc.read_h5ad(f'data/processed/srt/raw/{cid.lower()}.h5ad', backed='r')
scores = r.compute_module_scores_for_object(ad, 'counts', cid)
# emulate loop for one coarse/mid
mask = (scores[['T_NK','B_Plasma','Myeloid','DC_APC','CAF_Stromal','Endothelial','Tumor_like','Cycling','Low_quality_or_ambient']].idxmax(axis=1) == 'T_NK').values
print('T_NK mask sum', mask.sum())
cand_mid = [m for m,p in r.MID_PARENT.items() if p=='T_NK']
sub = scores.loc[mask, cand_mid]
best_mid = sub.idxmax(axis=1)
best_mid_score = sub.max(axis=1)
print('best_mid_score type', type(best_mid_score), best_mid_score.shape)
mid = pd.Series(['Unknown']*len(scores), index=scores.index)
mid.loc[mask] = np.where(best_mid_score > -0.5, best_mid, 'Unknown')
m='CD8_T'
mmask = mask & (mid == m).values
print('mmask sum', mmask.sum())
cand_fine = [f for f,p in r.FINE_PARENT.items() if p==m]
fsub = scores.loc[mmask, cand_fine]
best_fine = fsub.idxmax(axis=1)
best_fine_score = fsub.max(axis=1)
print('best_fine_score type', type(best_fine_score), best_fine_score.shape)
parent_score = scores.loc[mmask, m]
print('parent_score type', type(parent_score), parent_score.shape)
print('comparison', (best_fine_score > parent_score).head())
fine = pd.Series(['Unknown']*len(scores), index=scores.index)
fine.loc[mmask] = np.where(best_fine_score > parent_score, best_fine, 'Unknown')
print('fine assigned', fine.value_counts().head())
