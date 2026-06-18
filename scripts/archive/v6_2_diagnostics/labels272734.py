import pandas as pd
for chunk in pd.read_csv('results/v6_2/phase4a_cell_state_harmonization/handoff/cell_state_annotation_master.csv.gz', chunksize=50000, compression='gzip', low_memory=False):
    sub=chunk[chunk['cohort_id']=='GSE272734']
    if not sub.empty:
        print(sub[['marker_based_coarse_label','harmonized_coarse_label','harmonized_fine_label','final_label_level']].head(5).to_string(index=False))
        print('harmonized coarse counts:', sub['harmonized_coarse_label'].value_counts().to_dict())
        break
