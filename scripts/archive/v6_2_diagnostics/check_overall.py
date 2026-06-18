import pandas as pd
m='results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz'
total=0
unknown=0
coarse_known=0
mid_known=0
fine_known=0
by_tier={}
for chunk in pd.read_csv(m, chunksize=50000, compression='gzip', low_memory=False):
    total+=len(chunk)
    unknown+=(chunk['harmonized_fine_label']=='Unknown').sum()
    coarse_known+=(chunk['harmonized_coarse_label']!='Unknown').sum()
    mid_known+=(chunk['harmonized_mid_label']!='Unknown').sum()
    fine_known+=(chunk['harmonized_fine_label']!='Unknown').sum()
    for tier, g in chunk.groupby('integration_tier'):
        by_tier[tier]=by_tier.get(tier,0)+len(g)
print('total', total)
print('unknown_fine', unknown, 'pct', unknown/total*100)
print('coarse_known', coarse_known, coarse_known/total*100)
print('mid_known', mid_known, mid_known/total*100)
print('fine_known', fine_known, fine_known/total*100)
print('by_tier', by_tier)
