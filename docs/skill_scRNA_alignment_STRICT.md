# SKILL: strict-scrna-raw-alignment

## 1) Purpose
Align a new single-cell dataset into project raw SRT outputs for downstream merge:
- `data/processed/srt/raw/{dataset_id}.qs`
- `data/processed/srt/raw/{dataset_id}.h5ad`

Core goals are **non-negotiable**:
1. No data loss
2. No alignment drift
3. Uniform quality across datasets

---

## 2) Environment Constraints (MANDATORY)
1. Python/tooling must use conda env: `proj006`
   - Use only for python-side checks/tools (e.g., `anndata`, `h5py`, `ripgrep`).
2. **R must NOT run inside conda env**.
   - Always use external/system R environment.
3. If a required R package is missing in external R, **STOP and ask user** before installing or changing runtime.

---

## 3) Stop-and-Ask Policy (MANDATORY)
You must pause and ask user before deciding anything that can affect merge consistency, including but not limited to:
- metadata standard fields (required/optional set)
- metadata value normalization strategy
- filtering thresholds (cells/genes/biotypes/response filtering)
- duplicate gene conflict policy (if alternative exists)
- missing/ambiguous metadata mapping
- controlled vocabulary mapping (disease/organ/tissue/lib/response/timepoint/target)

If uncertain, ask with explicit options and impact:
1. Option A (conservative)
2. Option B (aggressive)
3. Recommendation and why

---

## 4) Canonical Pipeline (Invariant Order)
1. Read expression matrix (adapter-specific)
2. Gene symbol alignment (`gene_symbol_align`)
3. Keep gene types: `protein_coding`, `IG_C_gene`, `TR_C_gene`
4. Resolve duplicated aligned symbols by max `rowMeans`
5. Build/clean metadata to canonical schema
6. Enforce ID prefixing with `{dataset_id}_`
7. Build Seurat object (`min.cells=0`, `min.features=0` unless user changes)
8. Export `.qs` and `.h5ad`

---

## 5) Metadata Contract

### 5.1 Core fields (default baseline, verify with user)
- `cell_id`
- `sample_id`
- `donor_id`
- `disease`
- `organ`
- `tissue`
- `lib`

### 5.2 Common conditional fields
- `cell_sorting`, `sex`, `age`, `target`, `timepoint`
- `response`, `resp_standard`, `anno_orig`
- `surv`, `time_type`, `surv_type`, `surv_status`, `stage`

### 5.3 Alignment invariants
1. `colnames(counts) == rownames(meta_final)` in same order
2. `rownames(meta_final) == meta_final$cell_id`
3. `cell_id`, `sample_id`, `donor_id` carry dataset prefix
4. No post-subset mutation that breaks rowname alignment

---

## 6) Hard Validation Gates (must pass all)

### G1 Input integrity
- expression matrix readable
- dimensions > 0
- gene IDs and cell IDs non-empty

### G2 Gene alignment integrity
- aligned symbols non-empty rate acceptable (report exact numbers)
- duplicate symbol resolution complete
- final `rownames(counts_aligned)` unique

### G3 Metadata integrity
- all user-confirmed required fields exist
- no unintended NA spikes in required fields
- `nrow(meta_final) == ncol(counts_aligned)`

### G4 ID integrity
- `cell_id` unique
- prefix uniformity check for `cell_id/sample_id/donor_id`

### G5 Output integrity
- both output files exist and non-empty
- h5ad basic read check (python in `proj006`)

---

## 7) Known Drift/Bug Traps (must proactively check)
1. Mutating `meta` after creating `meta_final` (changes not propagated)
2. Mutating IDs after assigning `rownames(meta_final)`
3. Typo fields (e.g., `target` vs misspelled source column)
4. Duplicate field names in `selInfo`
5. Silent missing-field checks that do not stop execution

On any trap hit: **STOP and ask user**.

---

## 8) Adapter Design for New Dataset
Only these are dataset-specific:
1. `read_expression()`
2. `map_metadata()`
3. optional controlled-vocab mapper (must be user-confirmed)

Everything else must reuse invariant template.

---

## 9) Execution Checklist (copy/paste)
1. Confirm `dataset_id`, input paths, output paths
2. Confirm required metadata field set with user
3. Confirm filtering/threshold policy with user
4. Run invariant pipeline
5. Run G1-G5 gates
6. If any gate fails or ambiguous mapping appears: STOP + ask user
7. Export `.qs` + `.h5ad`
8. Emit final audit summary (counts, fields, missingness, ID checks)

---

## 10) Operational Commands

### Python/tooling (allowed in conda)
```bash
/home/huyudi/.local/share/r-miniconda/bin/conda run -n proj006 python -c "import anndata,h5py,pandas,numpy; print('ok')"
/home/huyudi/.local/share/r-miniconda/bin/conda run -n proj006 rg --version
```

### R runtime policy
- Use external R only (not `conda run -n proj006 Rscript`).

---

## 11) Fail-Safe Output Rule
If a choice can influence downstream merge consistency and user has not approved that policy yet:
- Do not guess
- Do not auto-fix
- Ask first, then continue
