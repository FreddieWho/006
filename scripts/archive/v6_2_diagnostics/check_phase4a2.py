import pandas as pd
m='results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz'
agg=[]
for chunk in pd.read_csv(m, chunksize=50000, compression='gzip', low_memory=False):
    agg.append(chunk.groupby(['cohort_id','harmonized_fine_label']).size())
s = pd.concat(agg).groupby(level=[0,1]).sum().reset_index(name='n_cells')
by_cohort = s.groupby('cohort_id')['n_cells'].sum().reset_index(name='total')
unknown = s[s['harmonized_fine_label']=='Unknown'].groupby('cohort_id')['n_cells'].sum().reset_index(name='unknown')
by_cohort = by_cohort.merge(unknown, on='cohort_id', how='left').fillna(0)
by_cohort['known_frac'] = 1 - by_cohort['unknown']/by_cohort['total']
print(by_cohort.sort_values('known_frac').head(10).to_string(index=False))
print('--- top known fine states ---')
print(s[s['harmonized_fine_label']!='Unknown'].sort_values('n_cells', ascending=False).head(20).to_string(index=False))
