import pandas as pd
cohorts=set()
for chunk in pd.read_csv('results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz', chunksize=50000, compression='gzip', low_memory=False):
    cohorts.update(chunk['cohort_id'].unique())
print('n', len(cohorts))
print(sorted(cohorts)[:20])
print('GSE272734' in cohorts)
