# Step2.7 Cell Fraction Feature Factory Third-Party Audit Report

- **audit_date**: 2026-05-21
- **auditor**: third-party review
- **target_run_id**: step2_v6_1_0505_0319
- **target_branch**: fraction_feature_v1
- **script**: scripts/v6_1/step2_7_cell_fraction_feature_factory_v1.py

## Executive Summary

Validation passed (2292 samples, 71 features, fraction sums ~1.0). No response leakage detected. No immediate data corruption. Three issues identified: one unused input path (minor), one label-set inconsistency w/ upstream (moderate), one directory naming mismatch (minor). Test suite could not be executed due to environment dependency conflict.

---

## Findings

### 1. FLAGS_PATH declared in manifest but never read in load_inputs() — Minor

**Location**: `step2_7_cell_fraction_feature_factory_v1.py` L26, L141-168, L505-507

**Description**: `FLAGS_PATH` points to `03_qc/final_cell_inclusion_flags.qc0513_balanced.parquet`. It is hashed and recorded in the decision manifest as an input, but `load_inputs()` reads `included_in_main_annotation` directly from the consensus parquet instead of loading the standalone QC flags file. This means the manifest claims a dependency that is not actually exercised at runtime. If the consensus parquet and the flags file ever diverge, provenance will be misleading.

**Recommended action**: Either load and cross-check `FLAGS_PATH` against consensus `included_in_main_annotation` to assert they match, or remove it from `input_paths` / `input_hashes` if it is not used.

---

### 2. Macro_inflammatory label inconsistency between PARENT_CHILDREN and MYELOID_MAIN_LABELS — Moderate

**Location**: `step2_7_cell_fraction_feature_factory_v1.py` L36-47, L62-79

**Description**: `PARENT_CHILDREN["Myeloid_DC"]` lists `"Macro_inflammatory"` as a valid child label, but `MYELOID_MAIN_LABELS` does **not** include it. Consequently, any cell adjudicated by Step2.6 as `Macro_inflammatory` would fail the `included_in_myeloid_main` gate (L195) and be downgraded to `Myeloid_unspecified` (L199-204). As a result, `Macro_inflammatory` never appears in `fraction_state_label`, and no `frac_parent_Myeloid_DC__Macro_inflammatory` feature is ever created.

**Data check**: Step2.6 output contains **zero** `Macro_inflammatory` cells (`final_myeloid_label` count = 0). However, the Step2.5 consensus `final_immune_mid_state` contains **4,779** `Macro_inflammatory` cells. These 4,779 cells were re-labeled during Step2.6 adjudication and therefore never propagate to Step2.7 under their original mid-state name.

**Impact**: Currently latent because Step2.6 does not emit `Macro_inflammatory`. If Step2.6 is ever updated to emit this label, Step2.7 will silently downgrade it to `Myeloid_unspecified`, and the intended parent feature will remain absent.

**Recommended action**:
- Add `"Macro_inflammatory"` to `MYELOID_MAIN_LABELS` if it is intended to be a main fraction label, **or**
- Remove `"Macro_inflammatory"` from `PARENT_CHILDREN["Myeloid_DC"]` if it is intentionally excluded from main fractions, **or**
- Document explicitly that `Macro_inflammatory` is excluded by design and is summarized only as QC/sensitivity.

---

### 3. Step2.8 plan document stored under `07_pseudobulk/` directory — Minor

**Location**: `results/v6_1/step2/plans/07_pseudobulk/STEP2_8_pseudobulk_feature_factory_plan.md`

**Description**: The plan for Step2.8 resides in a directory named `07_pseudobulk/`. There is no `08_pseudobulk/` directory. This breaks the convention that plan directory names match their step numbers (Step2.7 → `06_fraction/`, Step2.6 → `05_myeloid_qc/`). It is confusing for downstream navigation.

**Recommended action**: Rename directory to `08_pseudobulk/` and move the Step2.8 plan into it. Update any internal cross-references if they hard-code the old path.

---

### 4. Test suite cannot execute due to missing `requests_toolbelt` — Minor (Process)

**Location**: `mvp/tests/test_step2_7_cell_fraction_feature_factory_v1.py`

**Description**: `pytest` crashes during plugin collection because `langsmith` (pulled in via test dependency chain) requires `requests_toolbelt`, which is not installed in `proj006`. Passing `-p no:langsmith` does not prevent the import error because it occurs during plugin discovery before argument parsing.

**Recommended action**: Pin or install `requests_toolbelt` in the test environment, or isolate project tests from langsmith by running `pytest mvp/tests/test_step2_7_cell_fraction_feature_factory_v1.py --ignore-glob='*langsmith*'` / using a dedicated test requirements file.

---

## Validation Results (Reproduced from Output)

| Check | Result |
|-------|--------|
| Sample count matches expected | True (2292) |
| All features have dictionary row | True (71 features) |
| All-cell fraction sum max abs error | 5.55e-16 |
| Immune fraction sum max abs error | 5.55e-16 |
| Forbidden response columns | [] |
| Provenance columns present | True |
| Overall validation passed | True |

## Compliance Checklist

- [x] No response fields used or propagated.
- [x] No responder/non-responder comparison performed.
- [x] No differential analysis performed.
- [x] No model training performed.

## Recommended Priority Order

1. **Resolve Finding #2** (`Macro_inflammatory` label alignment) before any Step2.6 logic change that could start emitting this label.
2. **Resolve Finding #1** (unused `FLAGS_PATH`) to keep manifest accuracy.
3. **Resolve Finding #3** (directory rename) during next documentation pass.
4. **Resolve Finding #4** (test environment) to restore automated test coverage.
