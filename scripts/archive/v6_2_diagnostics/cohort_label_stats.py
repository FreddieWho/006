import pandas as pd
m='results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz'
agg=[]
for chunk in pd.read_csv(m, chunksize=50000, compression='gzip', low_memory=False):
    agg.append(chunk.groupby(['cohort_id','integration_tier']).agg(
        total=('cell_key','size'),
        unknown_fine=('harmonized_fine_label', lambda s: (s=='Unknown').sum()),
        coarse_known=('harmonized_coarse_label', lambda s: (s!='Unknown').sum()),
        mid_known=('harmonized_mid_label', lambda s: (s!='Unknown').sum()),
        fine_known=('harmonized_fine_label', lambda s: (s!='Unknown').sum()),
    ).reset_index())
s=pd.concat(agg).groupby(['cohort_id','integration_tier']).sum().reset_index()
s['fine_known_pct']=s['fine_known']/s['total']*100
s['mid_known_pct']=s['mid_known']/s['total']*100
s['coarse_known_pct']=s['coarse_known']/s['total']*100
s=s.sort_values('fine_known_pct')
print(s[['cohort_id','integration_tier','total','coarse_known_pct','mid_known_pct','fine_known_pct']].to_string(index=False))
