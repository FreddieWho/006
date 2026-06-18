import pandas as pd
m='results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz'
total=0
cohorts=set()
for chunk in pd.read_csv(m, chunksize=50000, compression='gzip', low_memory=False):
    total+=len(chunk)
    cohorts.update(chunk['cohort_id'].unique())
print('total cells', total)
print('cohorts', len(cohorts))
rel=pd.read_csv('results/v6_2/phase4a_cell_state_harmonization/qc/cell_state_reliability_scores.csv')
print('reliability rows', len(rel))
print(rel['reliability_level'].value_counts())
print('top fine states:')
print(rel.sort_values('n_cells', ascending=False).head(10)[['cohort_id','harmonized_fine_label','n_cells','reliability_level']].to_string(index=False))
