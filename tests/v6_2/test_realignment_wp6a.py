from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import yaml
import pytest


SCRIPT = Path(__file__).parents[2] / "scripts" / "v6_2" / "run_realignment_wp6a.py"
spec = importlib.util.spec_from_file_location("run_realignment_wp6a", SCRIPT)
wp6a = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(wp6a)


def _metadata_frames():
    registry = pd.DataFrame(
        {
            "cohort_id": ["C1", "C2"],
            "disease": ["d1", "d2"],
            "cancer_group": ["hcc_liver", "non_hcc"],
            "data_modality": ["scRNA", "bulkRNA"],
            "response_label_type": ["RECIST", "RECIST"],
            "response_label_quality": ["high", "medium"],
            "paired_available": [True, False],
            "spatial_available": [False, False],
            "tcr_raw_available": [False, False],
            "perturb_available": [False, False],
            "bulk_available": [False, True],
            "etiology_metadata": ["known", "unknown"],
            "b1_etiology_applicable": [False, False],
            "reaudit_status": ["RE-AUDIT", "RE-AUDIT"],
        }
    )
    environment = pd.DataFrame(
        {
            "cohort_id": ["C1", "C2"],
            "cancer_group": ["hcc_liver", "non_hcc"],
            "treatment_context": ["PD1_ICI_anchor", "PD1X_extension"],
            "response_label_type": ["RECIST", "RECIST"],
            "response_label_quality": ["high", "medium"],
            "timepoint_schema": ["pre/post", "unknown"],
            "paired_available": [True, False],
            "environment_role": ["other_or_external", "stratification_or_excluded"],
            "etiology_env_B1": ["not_applicable", "not_applicable"],
        }
    )
    return registry, environment


def test_wp6a_rejects_value_bearing_columns(tmp_path):
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"cohort_id": ["C1"], "response_raw": ["R"]}).to_csv(bad, index=False)
    with pytest.raises(ValueError, match="forbidden value-bearing"):
        wp6a.require_columns(bad, ["cohort_id"])


def test_wp6a_run_is_metadata_only(tmp_path):
    registry, environment = _metadata_frames()
    registry_path = tmp_path / "registry.csv"
    environment_path = tmp_path / "environment.csv"
    registry.to_csv(registry_path, index=False)
    environment.to_csv(environment_path, index=False)
    out = tmp_path / "out"
    wp6a.run(registry_path, environment_path, out, "test_wp6a")
    summary = yaml.safe_load((out / "anchor_metadata_summary.yaml").read_text())
    assert summary["status"] == "COMPLETE_METADATA_ONLY"
    assert summary["response_values_read"] is False
    assert summary["numeric_contrast_run"] is False
    assert summary["estimand_status"] == "NOT_ESTIMABLE_UNDER_RESPONSE_LOCK"
    audit = pd.read_csv(out / "anchor_metadata_audit.csv")
    assert len(audit) == 2
    assert audit["metadata_only"].all()
    assert not audit["response_values_read"].any()
