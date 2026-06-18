# Step2.6 Myeloid Adjudication Audit Findings

**Date**: 2026-05-19
**Auditor**: Third-party review
**Scope**: `scripts/v6_1/step2_6_myeloid_adjudication_v1.py` and outputs under `results/v6_1/step2/step2_v6_1_0505_0319/05_myeloid_qc/myeloid_adjudication_v1/`
**Run ID**: `step2_v6_1_0505_0319`

---

## Executive Summary

Core adjudication logic is correct and all row/key integrity validations pass. However, **one source (`krishna_2021_rcc.h5ad`) suffered total zero marker coverage due to a gene-naming mismatch**, causing 2,337 candidate cells to be silently downgraded to `Myeloid_unspecified`. Additional issues include panel-order override bias, a redundant TAM-like threshold, and an under-informative QC report.

**Verdict**: Conditional pass — fix the zero-coverage source before Step2.7.

---

## Findings

### 🔴 Critical: `krishna_2021_rcc.h5ad` Zero Marker Coverage

- **Manifest entry**: `source_h5ad: krishna_2021_rcc.h5ad`, `n_cells: 2337`, all `panel_gene_coverage: 0.0`
- **Root cause**: The h5ad stores gene identifiers in `var/gene_id` as Ensembl IDs (e.g. `ENSG00000237683`). `build_var_lookup()` only searches `var_names`, `gene_symbol`, `geneSymbol`, `gene_symbols`, `features`, `feature_name`, and `gene_name`. It does not look at `gene_id`, nor does it map Ensembl IDs to HGNC symbols.
- **Impact**: All 2,337 candidate cells from this source were adjudicated as `Myeloid_unspecified` (`downgrade_to_myeloid_unspecified`) because no marker genes could be resolved.
- **Location**: `scripts/v6_1/step2_6_myeloid_adjudication_v1.py:156-172`

### 🟡 Medium: `gse272735.h5ad` Low Panel Coverage

- **Coverage**: `cdc1` 33.3%, `macro` 60%, `tam_like` 66.7%, `mast_minor_immune` 20%, `nonimmune_exclusion` 44.4%
- **Root cause**: The dataset lacks many marker genes (verified: missing `CLEC9A`, `XCR1`, `C1QB`, `C1QC`, `SPP1`, `TREM2`, `TPSAB1`, etc.). This is a platform limitation, not a naming issue.
- **Impact**: Subtyping accuracy for this source is degraded. DC and macrophage calls rely on incomplete marker sets.

### 🟡 Medium: Panel Override Ordering Bias

- **Issue**: `specific_panels` is processed as a fixed list (`mono_inflammatory` → `macro` → `tam_like` → `cdc1` → `cdc2` → `pdc` → `mast_minor_immune`). When a cell satisfies multiple panels, the later panel always overwrites the earlier one.
- **Impact**: `mast_minor_immune` has the highest effective priority and `mono_inflammatory` the lowest, without biological justification.
- **Location**: `scripts/v6_1/step2_6_myeloid_adjudication_v1.py:248-265`

### 🟢 Minor: Redundant `tam_like` Threshold

- **Issue**: `tam_like` is processed with `strong=True`, which already enforces `z >= strong_z_min (1.0)` and `detected >= strong_detected_min (0.20)`. The additional `tam_like_z_min=0.75` check is logically masked by the stronger threshold and never affects the outcome.
- **Location**: `scripts/v6_1/step2_6_myeloid_adjudication_v1.py:256-259`

### 🟢 Minor: Candidate Selection Broader Than Plan

- **Issue**: The plan states review candidates should be `final_major_lineage == "Immune"` **with** myeloid/DC-compatible evidence. The code includes all `major == "Immune"` cells unconditionally.
- **Mitigation**: All 483,912 such cells happen to have `final_confidence == "low"`, so they are captured by `low_immune_review` and safely marked `sensitivity_only`. No incorrect main labels were produced.

### 🟢 Minor: Test Environment Failure

- **Issue**: `mvp/tests/test_step2_6_myeloid_adjudication_v1.py` fails because `anndata` is not installed in the current environment.
- **Note**: This is an environment issue, not a code bug.

### 🟢 Minor: QC Report Lacks Coverage Warnings

- **Issue**: `myeloid_qc_report.md` does not flag zero-coverage or low-coverage sources. A reader has no way to know that `krishna_2021_rcc` was effectively un-annotated.
- **Location**: `scripts/v6_1/step2_6_myeloid_adjudication_v1.py:484-511`

---

## Validation Results (Passed)

| Check | Result |
|---|---|
| Row count integrity (candidate = score = adjudication) | ✅ 1,522,971 |
| Candidate count matches locked plan | ✅ 1,522,971 |
| Unique composite key | ✅ True |
| Missing keys | ✅ False |
| Unknown sample_id | ✅ 0 |
| Allowed statuses only | ✅ True |
| Allowed labels only | ✅ True |
| Forbidden response columns | ✅ [] |
| No global clustering / no DEA / no response | ✅ True |

---

## Recommended Actions

1. **Critical**: Provide Ensembl→HGNC mapping for `krishna_2021_rcc.h5ad` and re-run scoring for that source, or explicitly exclude it from myeloid subtyping with a documented justification.
2. **High**: Add a zero/low-coverage source warning table to `myeloid_qc_report.md`.
3. **Medium**: Refactor `adjudicate_myeloid()` to use highest-supported score (or `top_panel`) as the tie-breaker instead of fixed sequential override.
4. **Low**: Remove redundant `tam_like` extra-threshold logic or adjust the `strong` flag so the extra threshold actually matters.
5. **Low**: Install `anndata` in the test environment or refactor tests to avoid importing the full script module.
