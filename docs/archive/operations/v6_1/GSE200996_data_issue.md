# GSE200996 Data Issue Note

## Summary
GSE200996 cannot be joined on raw 10x barcode alone.

## Why
- The scRNA matrices use standard 10x barcodes like `AAAAAAAAAAAAAA-1`.
- The metadata rows are keyed with patient/timepoint context, e.g. `barcode_patient_stage`.
- The dataset includes pooled PBMC libraries with hashtag demultiplexing, so the same 16nt barcode can occur in different patient pools.

## Consequence
- Stripping patient/stage context causes barcode collisions.
- Some cells are silently dropped.
- Some cells can be assigned to the wrong patient metadata.

## Handling rule
- Keep `patient + stage + barcode` as the join key.
- For pooled PBMC, use the hashtag mapping / pool context as part of the assignment logic.
- Do not use bare 16nt barcode as a global key.

## Notes
- This is a dataset-structure issue, not a corruption issue.
- The raw barcode whitelist is expected for 10x h5 outputs.

## Implementation update (2026-04-05)
- `scripts/preprocessing/script_sc/gse200996.R` now builds a context-aware join key as `matrix_context::barcode`.
- PBMC matrix context is resolved from `GSE200996_hashtag_info_sc_GEX_PBMC.csv` (pool column -> patient/stage sample IDs).
- Tumor matrix context is resolved from the matrix filename sample label.
- Global bare 16nt barcode matching is removed; joins are now context-scoped.
- Companion TCR export remains separate in `data/processed/tcr/raw/gse200996_tcr_cellid_barcode.tsv.gz`.
