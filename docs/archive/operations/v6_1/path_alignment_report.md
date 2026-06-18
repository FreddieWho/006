# path_alignment_report (v6.1 Task4)

## 结论
- 当前 `results/v6_1/step1/*` 主要承担的是合同 `S0 Cohort/Label/Context Lock` 的语义，而不是合同原文里的 `S1 Immune-state Measurement`。
- 因此下一个真实可执行步骤应读取 `results/v6_1/step1/*` 作为治理输入，并在 `results/v6_1/step2/*` 生成免疫状态测量产物。
- 合同中声明的 `results/v6_1/step1/*.parquet/csv` 免疫状态产物当前均未生成，不能用 symlink 伪造存在。

## Contract Path References Found
- `results/v6_1/step0/cohort_qc_audit.csv`
- `results/v6_1/step0/cohort_registry_v6_1.csv`
- `results/v6_1/step0/response_label_dictionary.csv`
- `results/v6_1/step0/treatment_context_flags.csv`
- `results/v6_1/step1/cell_state_fraction_by_sample.parquet`
- `results/v6_1/step1/immune_state_features.parquet`
- `results/v6_1/step1/pseudobulk_by_sample.parquet`
- `results/v6_1/step1/sample_index_master_v6_1.csv`
- `results/v6_1/step1/state_qc_report.csv`

## Contract-to-Actual Mapping
- contract: `results/v6_1/step0/cohort_registry_v6_1.csv`
  actual: `results/v6_1/step1/cohort_registry_v6_1.csv`
  status: `mapped_exact_semantics`
  symlink_candidate: `True`
  note: Actual step1 governance output semantically serves contract S0 cohort registry.
- contract: `results/v6_1/step0/treatment_context_flags.csv`
  actual: `results/v6_1/step1/treatment_context_flags.csv`
  status: `mapped_exact_semantics`
  symlink_candidate: `True`
  note: Actual step1 governance output semantically serves contract S0 treatment context lock.
- contract: `results/v6_1/step0/response_label_dictionary.csv`
  actual: `results/v6_1/step1/response_mapping_dictionary_v6_1.csv`
  status: `mapped_renamed_semantics_close`
  symlink_candidate: `True`
  note: Renamed in implementation; file includes binary mapping, rule type, evidence, and notes.
- contract: `results/v6_1/step0/cohort_qc_audit.csv`
  actual: `results/v6_1/step1/step1_summary.md;results/v6_1/step1/data_leakage_risk_log.md;results/v6_1/step1/cohort_role_gap_report.md`
  status: `no_single_equivalent`
  symlink_candidate: `False`
  note: Current implementation splits cohort audit across summary/risk/gap reports; do not symlink to a fake CSV.
- contract: `results/v6_1/step1/sample_index_master_v6_1.csv`
  actual: `results/v6_1/step1/sample_metadata_master_v6_1.csv`
  status: `proxy_only`
  symlink_candidate: `False`
  note: Sample metadata can seed a later sample index, but it is broader than a dedicated measurement-step sample index.
- contract: `results/v6_1/step1/cell_state_fraction_by_sample.parquet`
  actual: `none`
  status: `missing_not_generated`
  symlink_candidate: `False`
  note: Not yet generated; must be produced by the next executable measurement step.
- contract: `results/v6_1/step1/pseudobulk_by_sample.parquet`
  actual: `none`
  status: `missing_not_generated`
  symlink_candidate: `False`
  note: Not yet generated; MVP pseudobulk remains proxy-only by contract.
- contract: `results/v6_1/step1/immune_state_features.parquet`
  actual: `none`
  status: `missing_not_generated`
  symlink_candidate: `False`
  note: Not yet generated; must come from the next measurement step.
- contract: `results/v6_1/step1/state_qc_report.csv`
  actual: `none`
  status: `missing_not_generated`
  symlink_candidate: `False`
  note: Not yet generated; do not backfill from governance reports.

## Step2 Final Read Paths
- governance manifest: `docs/v6_1/step2_input_manifest_v6_1.yaml`
- Step2 should read governance tables from `results/v6_1/step1/`, sidecars from `mvp/sidecars/*/`, raw h5ad from `data/processed/srt/raw/`, and pathway resources from `data/pathway/`.

## Symlink Proposal (Do Not Execute Automatically)
- Safe symlink candidates:
  - `mkdir -p results/v6_1/step0`
  - `ln -s ../step1/cohort_registry_v6_1.csv results/v6_1/step0/cohort_registry_v6_1.csv`
  - `ln -s ../step1/treatment_context_flags.csv results/v6_1/step0/treatment_context_flags.csv`
  - `ln -s ../step1/response_mapping_dictionary_v6_1.csv results/v6_1/step0/response_label_dictionary.csv`
- Do not symlink:
  - `cohort_qc_audit.csv` because there is no single CSV equivalent.
  - any contract `S1` output (`cell_state_fraction_by_sample.parquet`, `pseudobulk_by_sample.parquet`, `immune_state_features.parquet`, `state_qc_report.csv`) because these files do not exist yet.
  - `sample_index_master_v6_1.csv` because `sample_metadata_master_v6_1.csv` is only a proxy, not an exact semantic match.

## Recommended Interpretation for Next Work
- Treat current implementation numbering as: actual `step1` = contract `S0`.
- Reserve actual `step2` for contract-like immune-state measurement outputs.
- Only after actual `step2` generates measurement matrices should strong baseline modeling begin.
