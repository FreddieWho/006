from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = ROOT / "docs" / "v6_1"
RESULTS_STEP1_DIR = ROOT / "results" / "v6_1" / "step1"
ANALYSIS_CONTRACT = DOCS_DIR / "analysis_contract.md"
SCOPE_LOCK = DOCS_DIR / "scope_lock_v6_1.md"

OUTPUT_MANIFEST = DOCS_DIR / "step2_input_manifest_v6_1.yaml"
OUTPUT_REPORT = DOCS_DIR / "path_alignment_report.md"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def quote_yaml(value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def list_files(base: Path, pattern: str) -> List[Path]:
    return sorted(base.glob(pattern))


def extract_contract_paths(text: str) -> List[str]:
    return sorted(set(re.findall(r"results/v6_1/step[01]/[A-Za-z0-9_./-]+", text)))


def build_contract_mapping() -> List[Dict[str, object]]:
    return [
        {
            "contract_path": "results/v6_1/step0/cohort_registry_v6_1.csv",
            "resolved_path": "results/v6_1/step1/cohort_registry_v6_1.csv",
            "exists": True,
            "alignment_status": "mapped_exact_semantics",
            "symlink_candidate": True,
            "notes": "Actual step1 governance output semantically serves contract S0 cohort registry.",
        },
        {
            "contract_path": "results/v6_1/step0/treatment_context_flags.csv",
            "resolved_path": "results/v6_1/step1/treatment_context_flags.csv",
            "exists": True,
            "alignment_status": "mapped_exact_semantics",
            "symlink_candidate": True,
            "notes": "Actual step1 governance output semantically serves contract S0 treatment context lock.",
        },
        {
            "contract_path": "results/v6_1/step0/response_label_dictionary.csv",
            "resolved_path": "results/v6_1/step1/response_mapping_dictionary_v6_1.csv",
            "exists": True,
            "alignment_status": "mapped_renamed_semantics_close",
            "symlink_candidate": True,
            "notes": "Renamed in implementation; file includes binary mapping, rule type, evidence, and notes.",
        },
        {
            "contract_path": "results/v6_1/step0/cohort_qc_audit.csv",
            "resolved_path": "results/v6_1/step1/step1_summary.md;results/v6_1/step1/data_leakage_risk_log.md;results/v6_1/step1/cohort_role_gap_report.md",
            "exists": True,
            "alignment_status": "no_single_equivalent",
            "symlink_candidate": False,
            "notes": "Current implementation splits cohort audit across summary/risk/gap reports; do not symlink to a fake CSV.",
        },
        {
            "contract_path": "results/v6_1/step1/sample_index_master_v6_1.csv",
            "resolved_path": "results/v6_1/step1/sample_metadata_master_v6_1.csv",
            "exists": True,
            "alignment_status": "proxy_only",
            "symlink_candidate": False,
            "notes": "Sample metadata can seed a later sample index, but it is broader than a dedicated measurement-step sample index.",
        },
        {
            "contract_path": "results/v6_1/step1/cell_state_fraction_by_sample.parquet",
            "resolved_path": "",
            "exists": False,
            "alignment_status": "missing_not_generated",
            "symlink_candidate": False,
            "notes": "Not yet generated; must be produced by the next executable measurement step.",
        },
        {
            "contract_path": "results/v6_1/step1/pseudobulk_by_sample.parquet",
            "resolved_path": "",
            "exists": False,
            "alignment_status": "missing_not_generated",
            "symlink_candidate": False,
            "notes": "Not yet generated; MVP pseudobulk remains proxy-only by contract.",
        },
        {
            "contract_path": "results/v6_1/step1/immune_state_features.parquet",
            "resolved_path": "",
            "exists": False,
            "alignment_status": "missing_not_generated",
            "symlink_candidate": False,
            "notes": "Not yet generated; must come from the next measurement step.",
        },
        {
            "contract_path": "results/v6_1/step1/state_qc_report.csv",
            "resolved_path": "",
            "exists": False,
            "alignment_status": "missing_not_generated",
            "symlink_candidate": False,
            "notes": "Not yet generated; do not backfill from governance reports.",
        },
    ]


def write_step2_manifest() -> None:
    sidecar_sample_summary = list_files(ROOT / "mvp" / "sidecars", "*/sample_summary.csv")
    sidecar_obs = list_files(ROOT / "mvp" / "sidecars", "*/obs.parquet")
    sidecar_var = list_files(ROOT / "mvp" / "sidecars", "*/var.parquet")
    raw_h5ad = list_files(ROOT / "data" / "processed" / "srt" / "raw", "*.h5ad")

    governance_inputs = [
        ROOT / "docs" / "data_collection.csv",
        RESULTS_STEP1_DIR / "cohort_registry_v6_1.csv",
        RESULTS_STEP1_DIR / "sample_metadata_master_v6_1.csv",
        RESULTS_STEP1_DIR / "patient_metadata_master_v6_1.csv",
        RESULTS_STEP1_DIR / "treatment_context_flags.csv",
        RESULTS_STEP1_DIR / "response_mapping_dictionary_v6_1.csv",
        RESULTS_STEP1_DIR / "frozen_cohort_inclusion_v6_1.csv",
        RESULTS_STEP1_DIR / "frozen_patient_split_v6_1.csv",
        RESULTS_STEP1_DIR / "source_records_local.tsv",
        RESULTS_STEP1_DIR / "source_records_web_verified.tsv",
    ]

    pathway_resources = [
        ROOT / "data" / "pathway" / "h.all.v2025.1.Hs.symbols.gmt",
        ROOT / "data" / "pathway" / "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
        ROOT / "data" / "pathway" / "c5.go.v2025.1.Hs.symbols.gmt",
        ROOT / "data" / "pathway" / "c7.immunesigdb.v2025.1.Hs.symbols.gmt",
        ROOT / "data" / "pathway" / "immport.gmt",
        ROOT / "data" / "pathway" / "immport_gene.gmt",
    ]

    lines: List[str] = []
    lines.append('manifest_version: "v6.1"')
    lines.append('step_id: "v6_1_step2_input_manifest"')
    lines.append('step_name: "immune_state_measurement_next_executable_step"')
    lines.append('contract_interpretation: "actual results/v6_1/step1 corresponds to contract S0 governance outputs"')
    lines.append('output_target_dir: "results/v6_1/step2"')
    lines.append("governance_inputs:")
    for path in governance_inputs:
        lines.append(f"  - path: {quote_yaml(rel(path))}")
        lines.append(f"    exists: {str(path.exists()).lower()}")
        lines.append('    role: "governance_table"')
    lines.append("sidecar_inputs:")
    lines.append("  sample_summary_csv:")
    for path in sidecar_sample_summary:
        lines.append(f"    - {quote_yaml(rel(path))}")
    lines.append("  obs_parquet:")
    for path in sidecar_obs:
        lines.append(f"    - {quote_yaml(rel(path))}")
    lines.append("  var_parquet:")
    for path in sidecar_var:
        lines.append(f"    - {quote_yaml(rel(path))}")
    lines.append("raw_single_cell_inputs:")
    for path in raw_h5ad:
        lines.append(f"  - {quote_yaml(rel(path))}")
    lines.append("feature_resource_inputs:")
    for path in pathway_resources:
        lines.append(f"  - path: {quote_yaml(rel(path))}")
        lines.append(f"    exists: {str(path.exists()).lower()}")
    lines.append("pending_contract_s1_outputs:")
    for item in [
        "results/v6_1/step1/cell_state_fraction_by_sample.parquet",
        "results/v6_1/step1/pseudobulk_by_sample.parquet",
        "results/v6_1/step1/immune_state_features.parquet",
        "results/v6_1/step1/state_qc_report.csv",
    ]:
        lines.append(f"  - contract_path: {quote_yaml(item)}")
        lines.append('    exists: false')
        lines.append('    note: "Must be generated by the next executable measurement step before strict contract-S2 baseline modeling."')

    OUTPUT_MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report() -> None:
    contract_text = ANALYSIS_CONTRACT.read_text(encoding="utf-8")
    scope_text = SCOPE_LOCK.read_text(encoding="utf-8")
    contract_paths = extract_contract_paths(contract_text + "\n" + scope_text)
    mappings = build_contract_mapping()

    lines: List[str] = []
    lines.append("# path_alignment_report (v6.1 Task4)")
    lines.append("")
    lines.append("## 结论")
    lines.append("- 当前 `results/v6_1/step1/*` 主要承担的是合同 `S0 Cohort/Label/Context Lock` 的语义，而不是合同原文里的 `S1 Immune-state Measurement`。")
    lines.append("- 因此下一个真实可执行步骤应读取 `results/v6_1/step1/*` 作为治理输入，并在 `results/v6_1/step2/*` 生成免疫状态测量产物。")
    lines.append("- 合同中声明的 `results/v6_1/step1/*.parquet/csv` 免疫状态产物当前均未生成，不能用 symlink 伪造存在。")
    lines.append("")
    lines.append("## Contract Path References Found")
    for item in contract_paths:
        lines.append(f"- `{item}`")
    lines.append("")
    lines.append("## Contract-to-Actual Mapping")
    for item in mappings:
        lines.append(f"- contract: `{item['contract_path']}`")
        lines.append(f"  actual: `{item['resolved_path'] or 'none'}`")
        lines.append(f"  status: `{item['alignment_status']}`")
        lines.append(f"  symlink_candidate: `{item['symlink_candidate']}`")
        lines.append(f"  note: {item['notes']}")
    lines.append("")
    lines.append("## Step2 Final Read Paths")
    lines.append(f"- governance manifest: `{rel(OUTPUT_MANIFEST)}`")
    lines.append("- Step2 should read governance tables from `results/v6_1/step1/`, sidecars from `mvp/sidecars/*/`, raw h5ad from `data/processed/srt/raw/`, and pathway resources from `data/pathway/`.")
    lines.append("")
    lines.append("## Symlink Proposal (Do Not Execute Automatically)")
    lines.append("- Safe symlink candidates:")
    lines.append("  - `mkdir -p results/v6_1/step0`")
    lines.append("  - `ln -s ../step1/cohort_registry_v6_1.csv results/v6_1/step0/cohort_registry_v6_1.csv`")
    lines.append("  - `ln -s ../step1/treatment_context_flags.csv results/v6_1/step0/treatment_context_flags.csv`")
    lines.append("  - `ln -s ../step1/response_mapping_dictionary_v6_1.csv results/v6_1/step0/response_label_dictionary.csv`")
    lines.append("- Do not symlink:")
    lines.append("  - `cohort_qc_audit.csv` because there is no single CSV equivalent.")
    lines.append("  - any contract `S1` output (`cell_state_fraction_by_sample.parquet`, `pseudobulk_by_sample.parquet`, `immune_state_features.parquet`, `state_qc_report.csv`) because these files do not exist yet.")
    lines.append("  - `sample_index_master_v6_1.csv` because `sample_metadata_master_v6_1.csv` is only a proxy, not an exact semantic match.")
    lines.append("")
    lines.append("## Recommended Interpretation for Next Work")
    lines.append("- Treat current implementation numbering as: actual `step1` = contract `S0`.")
    lines.append("- Reserve actual `step2` for contract-like immune-state measurement outputs.")
    lines.append("- Only after actual `step2` generates measurement matrices should strong baseline modeling begin.")

    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    write_step2_manifest()
    write_report()


if __name__ == "__main__":
    main()
