from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Set

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
STEP1_DIR = ROOT / "results" / "v6_1" / "step1"

REGISTRY_PATH = STEP1_DIR / "cohort_registry_v6_1.csv"
ROLE_PATH = STEP1_DIR / "dataset_role_assignment.csv"
SAMPLE_PATH = STEP1_DIR / "sample_metadata_master_v6_1.csv"
SOURCE_LOCAL_PATH = STEP1_DIR / "source_records_local.tsv"
SOURCE_WEB_PATH = STEP1_DIR / "source_records_web_verified.tsv"

OUTPUT_CSV = STEP1_DIR / "frozen_cohort_inclusion_v6_1.csv"
OUTPUT_MD = STEP1_DIR / "cohort_role_gap_report.md"


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "n/a", "none", "null"}:
        return ""
    return text


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def parse_roles(value: object) -> List[str]:
    text = clean_str(value)
    if not text:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def parse_paths(value: object) -> List[Path]:
    text = clean_str(value)
    if not text:
        return []
    return [Path(x.strip()) for x in text.split(";") if x.strip()]


def has_existing_path(values: Iterable[object]) -> bool:
    for value in values:
        for path in parse_paths(value):
            if path.exists():
                return True
    return False


def build_sample_summary(sample: pd.DataFrame) -> pd.DataFrame:
    def count_true(series: pd.Series) -> int:
        return int(series.astype(str).str.lower().isin(["true", "1", "yes", "y", "t"]).sum())

    return (
        sample.groupby("cohort_id")
        .agg(
            sample_count=("sample_id", "size"),
            patient_count=("patient_id", "nunique"),
            usable_main_samples=("usable_in_main_analysis", count_true),
            usable_sensitivity_samples=("usable_in_sensitivity_analysis", count_true),
            scRNA_samples=("scRNA_available", count_true),
            spatial_samples=("spatial_available", count_true),
            bulk_samples=("bulk_available", count_true),
        )
        .reset_index()
    )


def choose_executable_roles(
    declared_roles: List[str],
    treatment_context: str,
    sample_count: int,
    usable_main_samples: int,
    usable_sensitivity_samples: int,
    spatial_samples: int,
    has_bulk_input: bool,
    has_perturb_input: bool,
) -> List[str]:
    roles: List[str] = []

    if "perturb_prior" in declared_roles and has_perturb_input:
        roles.append("perturb_prior")

    if "bulk_external_anchor" in declared_roles and has_bulk_input:
        roles.append("bulk_external_anchor")

    if "spatial_adjudication" in declared_roles and spatial_samples > 0:
        roles.append("spatial_adjudication")

    if sample_count > 0:
        if "main_scRNA_training" in declared_roles and usable_main_samples > 0:
            roles.append("main_scRNA_training")
        if "PD1_anchor_training" in declared_roles and treatment_context == "PD1_ICI_anchor" and usable_main_samples > 0:
            roles.append("PD1_anchor_training")
        if "PD1X_extension_support" in declared_roles and treatment_context == "PD1X_extension" and usable_main_samples > 0:
            roles.append("PD1X_extension_support")
        if "HCC_specific_deep_dive" in declared_roles and (
            usable_main_samples > 0 or usable_sensitivity_samples > 0 or spatial_samples > 0
        ):
            roles.append("HCC_specific_deep_dive")
        if "sensitivity_only" in declared_roles and usable_sensitivity_samples > 0:
            roles.append("sensitivity_only")

    return roles


def decide_inclusion_status(
    declared_roles: List[str],
    executable_roles: List[str],
    sample_count: int,
    usable_main_samples: int,
    usable_sensitivity_samples: int,
    spatial_samples: int,
    has_any_local_input: bool,
    has_local_source_record: bool,
    has_web_verification: bool,
) -> str:
    if "excluded" in declared_roles:
        return "excluded_no_usable_input"

    if "perturb_prior" in executable_roles:
        return "perturb_prior_only"

    if "spatial_adjudication" in executable_roles and usable_main_samples == 0:
        return "spatial_only"

    if "bulk_external_anchor" in executable_roles and usable_main_samples == 0 and usable_sensitivity_samples == 0:
        return "external"

    if any(role in executable_roles for role in ["main_scRNA_training", "PD1_anchor_training", "PD1X_extension_support"]):
        return "main"

    if usable_sensitivity_samples > 0 or "sensitivity_only" in executable_roles:
        return "sensitivity"

    if sample_count > 0 or has_any_local_input or has_local_source_record or has_web_verification:
        return "excluded_pending_metadata"

    return "excluded_no_usable_input"


def build_status_reason(
    inclusion_status: str,
    declared_roles: List[str],
    executable_roles: List[str],
    sample_count: int,
    usable_main_samples: int,
    usable_sensitivity_samples: int,
    spatial_samples: int,
    has_any_local_input: bool,
    local_source_count: int,
    web_source_count: int,
) -> str:
    if inclusion_status == "main":
        return (
            f"sample_ingested={sample_count};usable_main_samples={usable_main_samples};"
            f"executable_roles={';'.join(executable_roles)}"
        )
    if inclusion_status == "sensitivity":
        return (
            f"sample_ingested={sample_count};usable_sensitivity_samples={usable_sensitivity_samples};"
            f"executable_roles={';'.join(executable_roles) or 'sensitivity_only'}"
        )
    if inclusion_status == "external":
        return (
            f"bulk_or_external_input_available={has_any_local_input};"
            f"local_source_records={local_source_count};web_records={web_source_count}"
        )
    if inclusion_status == "spatial_only":
        return (
            f"spatial_samples={spatial_samples};usable_main_samples={usable_main_samples};"
            f"declared_roles={';'.join(declared_roles)}"
        )
    if inclusion_status == "perturb_prior_only":
        return (
            f"perturb_input_available={has_any_local_input};"
            f"declared_roles={';'.join(declared_roles)}"
        )
    if inclusion_status == "excluded_pending_metadata":
        return (
            f"declared_but_not_executable;sample_ingested={sample_count};usable_main_samples={usable_main_samples};"
            f"usable_sensitivity_samples={usable_sensitivity_samples};local_input={has_any_local_input};"
            f"local_source_records={local_source_count};web_records={web_source_count}"
        )
    return f"no_executable_input_or_explicitly_excluded;declared_roles={';'.join(declared_roles)}"


def build_gap_report(df: pd.DataFrame) -> str:
    lines: List[str] = []
    lines.append("# cohort_role_gap_report (v6.1 Task2)")
    lines.append("")
    lines.append("## 冻结原则")
    lines.append("- `inclusion_status` 代表当前工作区下的可执行状态，而不是概念上可能获得的数据。")
    lines.append("- `main` 仅授予已入库且存在 `usable_main_samples > 0` 的主 scRNA cohort。")
    lines.append("- `sensitivity` 仅授予已入库且存在 `usable_sensitivity_samples > 0` 的高混杂/支持性 cohort。")
    lines.append("- `external` 仅授予已有本地 bulk/clinical 输入的外部 anchor。")
    lines.append("- `spatial_only` 与 `perturb_prior_only` 单独冻结，避免误入主训练。")
    lines.append("- 缺少样本主表、缺少可执行输入、或 metadata 仍不足者，冻结为 `excluded_pending_metadata` 或 `excluded_no_usable_input`。")
    lines.append("")

    lines.append("## inclusion_status 统计")
    for status, count in df["inclusion_status"].value_counts().sort_index().items():
        lines.append(f"- `{status}`: {int(count)} cohorts")
    lines.append("")

    role_specs = [
        ("HCC deep-dive", "HCC_specific_deep_dive"),
        ("shared module", "main_scRNA_training"),
        ("bulk external anchor", "bulk_external_anchor"),
        ("spatial adjudication", "spatial_adjudication"),
        ("perturb prior", "perturb_prior"),
    ]

    for label, role_key in role_specs:
        declared = df[df["declared_roles"].str.contains(role_key, na=False)].copy()
        if role_key == "HCC_specific_deep_dive":
            ingested = declared[declared["sample_ingested"]]
        elif role_key == "main_scRNA_training":
            ingested = declared[declared["sample_ingested"]]
        elif role_key == "bulk_external_anchor":
            ingested = declared[declared["inclusion_status"] == "external"]
        elif role_key == "spatial_adjudication":
            ingested = declared[declared["inclusion_status"] == "spatial_only"]
        else:
            ingested = declared[declared["inclusion_status"] == "perturb_prior_only"]

        gap = declared[~declared["cohort_id"].isin(ingested["cohort_id"])]

        lines.append(f"## {label}")
        lines.append(f"- 声明 cohort: `{len(declared)}`")
        lines.append(f"- 已入库 cohort: `{len(ingested)}`")
        lines.append(f"- 缺口 cohort: `{len(gap)}`")
        lines.append(f"- declared: `{';'.join(declared['cohort_id'].tolist()) or 'none'}`")
        lines.append(f"- ingested: `{';'.join(ingested['cohort_id'].tolist()) or 'none'}`")
        lines.append(f"- gap: `{';'.join(gap['cohort_id'].tolist()) or 'none'}`")
        lines.append("")

    lines.append("## 主要阻断")
    blocked = df[df["inclusion_status"] == "excluded_pending_metadata"].copy()
    blocked = blocked.sort_values(
        ["primary_role", "sample_count", "has_any_local_input"],
        ascending=[True, False, False],
    )
    for _, row in blocked.head(15).iterrows():
        lines.append(
            f"- `{row['cohort_id']}`: primary_role=`{row['primary_role']}`, "
            f"declared_roles=`{row['declared_roles']}`, reason=`{row['status_reason']}`"
        )
    lines.append("")

    lines.append("## 结论")
    lines.append("- `shared module` 与 `HCC deep-dive` 已有可执行主干，但 cohort 声明数明显高于当前已入库数。")
    lines.append("- `bulk external anchor` 目前只有已落地本地输入的 cohort 应进入后续 external anchor 步骤，其余保留为 gap。")
    lines.append("- `spatial adjudication` 当前只有一个已冻结 cohort，足够作为 biological adjudicator 原型，但不构成广覆盖验证面。")
    lines.append("- `perturb prior` 角色在 registry 层定义完整，但需要严格与主训练隔离。")
    return "\n".join(lines) + "\n"


def main() -> None:
    registry = pd.read_csv(REGISTRY_PATH)
    roles = pd.read_csv(ROLE_PATH)
    sample = pd.read_csv(SAMPLE_PATH)
    source_local = pd.read_csv(SOURCE_LOCAL_PATH, sep="\t")
    source_web = pd.read_csv(SOURCE_WEB_PATH, sep="\t")

    sample_summary = build_sample_summary(sample)
    df = registry.merge(roles, on="cohort_id", how="left").merge(sample_summary, on="cohort_id", how="left")
    df = df.fillna(
        {
            "assigned_roles": "",
            "primary_role": "unknown",
            "rationale": "",
            "sample_count": 0,
            "patient_count": 0,
            "usable_main_samples": 0,
            "usable_sensitivity_samples": 0,
            "scRNA_samples": 0,
            "spatial_samples": 0,
            "bulk_samples": 0,
        }
    )

    local_counts = source_local["cohort_id"].value_counts().to_dict()
    web_counts = source_web["cohort_id"].value_counts().to_dict()

    rows: List[Dict[str, object]] = []
    for _, row in df.iterrows():
        cohort_id = row["cohort_id"]
        declared_roles = parse_roles(row["assigned_roles"])
        sample_count = int(row["sample_count"])
        usable_main_samples = int(row["usable_main_samples"])
        usable_sensitivity_samples = int(row["usable_sensitivity_samples"])
        spatial_samples = int(row["spatial_samples"])

        has_local_source_record = cohort_id in local_counts
        has_web_verification = cohort_id in web_counts
        has_any_local_input = has_existing_path(
            [
                row["raw_data_path"],
                row["processed_data_path"],
                row["clinical_metadata_path"],
                row["sidecar_path"],
            ]
        )
        has_bulk_input = truthy(row["bulk_available"]) and (
            has_any_local_input or int(row["bulk_samples"]) > 0 or has_local_source_record
        )
        has_perturb_input = truthy(row["perturb_available"]) and has_any_local_input

        executable_roles = choose_executable_roles(
            declared_roles=declared_roles,
            treatment_context=clean_str(row["treatment_context"]),
            sample_count=sample_count,
            usable_main_samples=usable_main_samples,
            usable_sensitivity_samples=usable_sensitivity_samples,
            spatial_samples=spatial_samples,
            has_bulk_input=has_bulk_input,
            has_perturb_input=has_perturb_input,
        )
        inclusion_status = decide_inclusion_status(
            declared_roles=declared_roles,
            executable_roles=executable_roles,
            sample_count=sample_count,
            usable_main_samples=usable_main_samples,
            usable_sensitivity_samples=usable_sensitivity_samples,
            spatial_samples=spatial_samples,
            has_any_local_input=has_any_local_input,
            has_local_source_record=has_local_source_record,
            has_web_verification=has_web_verification,
        )
        status_reason = build_status_reason(
            inclusion_status=inclusion_status,
            declared_roles=declared_roles,
            executable_roles=executable_roles,
            sample_count=sample_count,
            usable_main_samples=usable_main_samples,
            usable_sensitivity_samples=usable_sensitivity_samples,
            spatial_samples=spatial_samples,
            has_any_local_input=has_any_local_input,
            local_source_count=int(local_counts.get(cohort_id, 0)),
            web_source_count=int(web_counts.get(cohort_id, 0)),
        )

        rows.append(
            {
                "cohort_id": cohort_id,
                "cohort_name": row["cohort_name"],
                "treatment_context": row["treatment_context"],
                "declared_roles": ";".join(declared_roles),
                "primary_role": row["primary_role"],
                "executable_role": ";".join(executable_roles) if executable_roles else "none",
                "inclusion_status": inclusion_status,
                "sample_ingested": sample_count > 0,
                "sample_count": sample_count,
                "patient_count": int(row["patient_count"]),
                "usable_main_samples": usable_main_samples,
                "usable_sensitivity_samples": usable_sensitivity_samples,
                "has_any_local_input": has_any_local_input,
                "has_local_source_record": has_local_source_record,
                "has_web_verification": has_web_verification,
                "status_reason": status_reason,
                "rationale": row["rationale"],
                "notes": clean_str(row["notes"]),
            }
        )

    out = pd.DataFrame(rows).sort_values(["inclusion_status", "primary_role", "cohort_id"]).reset_index(drop=True)
    out.to_csv(OUTPUT_CSV, index=False)
    OUTPUT_MD.write_text(build_gap_report(out), encoding="utf-8")


if __name__ == "__main__":
    main()
