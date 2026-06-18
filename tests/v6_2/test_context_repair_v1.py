from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def load_module(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


def test_response_mapping_preserves_ordinal_and_ne():
    rules = load_module("context_repair_rules", "scripts/v6_2/context_repair/cohort_rules.py")
    assert rules.response_binary("High", "GSE200996")[0] == "unknown"
    assert rules.response_binary("Medium", "GSE200996")[0] == "unknown"
    assert rules.response_binary("NE", "GSE266919")[0] == "unknown"
    assert rules.response_binary("PD", "GSE266919")[0] == "non_responder"
    assert rules.response_binary("PR", "TASK01")[0] == "responder"
    assert rules.response_binary("SD", "TASK01")[0] == "non_responder"


def test_task01_timepoints_use_source_proven_visit_order():
    common = load_module("phase7_common_task01", "scripts/v6_2/phase7_common.py")
    assert common.normalize_timepoint("T0", "TASK01") == "baseline"
    assert common.normalize_timepoint("T1", "TASK01") == "on_treatment"
    assert common.normalize_timepoint("T2", "TASK01") == "post_treatment"
    assert common.normalize_timepoint("T0", "other_cohort") == "t0"


def test_local_h5ad_modalities_are_not_legacy_bulk_or_spatial():
    rules = load_module("context_repair_rules_modality", "scripts/v6_2/context_repair/cohort_rules.py")
    assert rules.COHORT_CORRECTIONS["GSE200996"]["modality"] == "scRNA"
    assert rules.COHORT_CORRECTIONS["GSE179994"]["modality"] == "scRNA"
    assert rules.COHORT_CORRECTIONS["mendeley_skrx2fz79n"]["modality"] == "scRNA"


def test_su010_has_one_study_subject_key():
    rules = load_module("context_repair_rules_subject", "scripts/v6_2/context_repair/cohort_rules.py")
    bcc = {"cohort_id": "GSE123813_bcc", "patient_id": "su010", "patient_key": "GSE123813_bcc::su010"}
    scc = {"cohort_id": "GSE123813_scc", "patient_id": "su010", "patient_key": "GSE123813_scc::su010"}
    assert rules.default_study_subject_key(bcc) == "GSE123813::su010"
    assert rules.default_study_subject_key(scc) == "GSE123813::su010"


def test_author_label_mapping_is_response_blind_and_conservative():
    repair = load_module("rebuild_phase4a4b", "scripts/v6_2/context_repair/rebuild_phase4a4b.py")
    mapped = repair.author_annotation_labels(pd.Series(["CD8_tumor:1", "NK cell", "Myeloid", "mystery"]))
    assert mapped.coarse.tolist() == ["T_NK", "T_NK", "Myeloid", "Unknown"]
    assert mapped.mid.tolist() == ["CD8_T", "NK", "Unknown", "Unknown"]


def test_identity_overlay_never_changes_cell_count():
    repair = load_module("rebuild_phase4a4b_overlay", "scripts/v6_2/context_repair/rebuild_phase4a4b.py")
    master = pd.DataFrame({
        "cohort_id": ["X", "X"], "source_object_id": ["X::object_1"] * 2,
        "original_cell_barcode": ["a", "b"], "sample_key": ["old", "old"],
        "patient_key": ["X::p", "X::p"],
    })
    identity = pd.DataFrame({
        "cohort_id": ["X", "X"], "source_object_id": ["X::object_1"] * 2,
        "source_cell_id": ["a", "b"], "study_subject_key": ["X::p", "X::p"],
        "analysis_unit_key": ["X::p::baseline::tumor"] * 2,
        "normalized_timepoint": ["baseline"] * 2, "tissue_context": ["tumor"] * 2,
        "assignment_status": ["resolved", "resolved"],
    })
    fixed, audit = repair.overlay_group(master, repair.normalize_identity(identity, "X"))
    assert len(fixed) == len(master)
    assert audit["status"] == "PASS"
    assert fixed.sample_key.nunique() == 1


def test_identity_overlay_falls_back_to_frozen_patient_sample_key():
    repair = load_module("rebuild_phase4a4b_fallback", "scripts/v6_2/context_repair/rebuild_phase4a4b.py")
    master = pd.DataFrame({
        "cohort_id": ["X"], "source_object_id": ["X::object_1"],
        "original_cell_barcode": ["a"], "sample_key": ["X::sample1"],
        "patient_key": ["X::patient1"],
    })
    identity = pd.DataFrame({
        "source_cell_id": ["a"], "source_object_id": ["X::object_1"],
        "study_subject_key": ["unknown"], "patient_key": ["unknown"],
        "analysis_unit_key": ["analysis::unknown"], "normalized_timepoint": ["unknown"],
        "tissue_context": ["unknown"], "assignment_status": ["quarantined"],
    })
    fixed, audit = repair.overlay_group(master, repair.normalize_identity(identity, "X"))
    assert fixed.loc[0, "patient_key"] == "X::patient1"
    assert fixed.loc[0, "sample_key"] == "X::sample1"
    assert fixed.loc[0, "assignment_status"] == "resolved_legacy_validated_fallback"
    assert audit["n_legacy_fallback"] == 1


def test_phase7_config_environment_override(tmp_path, monkeypatch):
    cfg = {
        "output_dir": "results/v6_2/test_override",
        "inputs": {},
    }
    path = tmp_path / "phase7.yaml"
    path.write_text(__import__("yaml").safe_dump(cfg))
    monkeypatch.setenv("PHASE7_CONFIG_PATH", str(path))
    module = load_module("phase7_common_override", "scripts/v6_2/phase7_common.py")
    loaded = module.load_config()
    assert loaded["out"] == ROOT / "results/v6_2/test_override"


def test_treatment_arm_is_part_of_analysis_unit_identity():
    identity = load_module("context_repair_identity_arm", "scripts/v6_2/context_repair/identity.py")
    obs = pd.DataFrame({
        "patient_id": ["14-1189", "14-1189"],
        "timepoint": ["Follow-up 4", "Follow-up 4"],
        "tissue_source": ["tumor", "tumor"],
        "treatment": ["aPD1", "aCTLA-4"],
        "sample_id": ["pd1", "ctla"],
    }, index=["cell-a", "cell-b"])
    fixed = identity.build_cell_identity(obs, "GSE272993", "GSE272993::object_1")
    assert fixed.study_subject_key.nunique() == 1
    assert fixed.analysis_unit_key.nunique() == 2
    assert fixed.biological_sample_key.nunique() == 2
    assert set(fixed.treatment_arm) == {"aPD1", "aCTLA-4"}


def test_local_filename_cohort_aliases_are_stable():
    identity = load_module("context_repair_identity_alias", "scripts/v6_2/context_repair/identity.py")
    empty = pd.DataFrame(index=["cell"])
    assert identity.cohort_from_obs("mendeley_skrx2fz79n", empty) == "mendeley_skrx2fz79n"
    assert identity.cohort_from_obs("bi2021rcc", empty) == "bi_2021_rcc"
    assert identity.cohort_from_obs("gse200996", empty) == "GSE200996"


def test_generic_identity_accepts_source_proven_donor_and_tissue_aliases():
    identity = load_module("context_repair_identity_donor", "scripts/v6_2/context_repair/identity.py")
    obs = pd.DataFrame({
        "donor_id": ["P1", "P1"],
        "sample_id": ["S1", "S1"],
        "lib": ["L1", "L1"],
        "timepoint": ["pre", "pre"],
        "tissue": ["tumor", "tumor"],
    }, index=["c1", "c2"])
    fixed = identity.build_cell_identity(obs, "TASK01", "TASK01::object_1")
    assert fixed.assignment_status.eq("resolved").all()
    assert fixed.study_subject_key.eq("TASK01::P1").all()
    assert fixed.normalized_timepoint.eq("baseline").all()
    assert fixed.tissue_context.eq("tumor").all()


def test_frozen_projection_supports_raw_counts_in_x_and_preserves_all_modules(tmp_path):
    import anndata as ad
    import numpy as np

    module_dir = ROOT / "scripts/v6_2/context_repair"
    sys.path.insert(0, str(module_dir))
    try:
        import frozen_projection as projection
    finally:
        sys.path.pop(0)
    h5ad = tmp_path / "mini.h5ad"
    a = ad.AnnData(
        X=np.array([[1, 0, 2], [0, 3, 0], [2, 1, 0]], dtype=np.int32),
        obs=pd.DataFrame(index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["A", "B", "C"]),
    )
    a.write_h5ad(h5ad)
    identity = pd.DataFrame({
        "source_object_id": ["X::object_1"] * 3,
        "source_cell_id": ["c1", "c2", "c3"],
        "study_subject_key": ["X::p1"] * 3,
        "analysis_unit_key": ["analysis::1"] * 3,
        "normalized_timepoint": ["baseline"] * 3,
        "tissue_context": ["tumor"] * 3,
        "assignment_status": ["resolved"] * 3,
    })
    identity_path = tmp_path / "identity.parquet"
    identity.to_parquet(identity_path)
    annotation = pd.DataFrame({
        "source_object_id": ["X::object_1"] * 3,
        "source_cell_id": ["c1", "c2", "c3"],
        "harmonized_coarse_label": ["T_NK"] * 3,
        "harmonized_mid_label": ["CD8_T"] * 3,
    })
    annotation_path = tmp_path / "annotation.parquet"
    annotation.to_parquet(annotation_path)
    membership = pd.DataFrame({
        "frozen_module_id": ["FM01", "FM02"],
        "gene": ["A", "MISSING"],
        "relative_weight": [1.0, 1.0],
    })
    membership_path = tmp_path / "membership.csv"
    membership.to_csv(membership_path, index=False)
    reference = pd.DataFrame({"expression_unit_id": ["old"], "A": [0.0], "B": [0.0]})
    reference_path = tmp_path / "reference.parquet"
    reference.to_parquet(reference_path)
    config_path = tmp_path / "projection.yaml"
    config_path.write_text(__import__("yaml").safe_dump({"inputs": {
        "objects": [{"cohort_id": "X", "object_id": "X::object_1", "h5ad": str(h5ad),
                     "count_layer": "X", "identity_sidecar": str(identity_path)}],
        "phase4a_annotations": str(annotation_path), "membership": str(membership_path),
        "expression_gene_reference": str(reference_path),
        "coarse_state_column": "harmonized_coarse_label", "mid_state_column": "harmonized_mid_label",
    }}))
    config, inputs = projection.load_config(config_path)
    out = tmp_path / "out"
    projection.run_projection(config, inputs, out)
    score = pd.read_parquet(out / "frozen_module_score_matrix.parquet")
    registry = pd.read_csv(out / "unified_expression_unit_registry.csv")
    assert {"FM01", "FM02"}.issubset(score.columns)
    assert score.FM02.isna().all()
    assert set(score.patient_key) == {"X::p1"}
    assert {"n_genes", "inclusion_status", "layer_family"}.issubset(registry.columns)


def test_phase4b_exact_legacy_sample_binding_preserves_cross_arm_labels(tmp_path):
    repair = load_module("rebuild_phase4b_exact_binding", "scripts/v6_2/context_repair/rebuild_phase4a4b.py")
    out4a, metadata, out4b = tmp_path / "p4a", tmp_path / "meta", tmp_path / "p4b"
    (out4a / "handoff/cell_state_annotation_master.parquet").mkdir(parents=True)
    metadata.mkdir()
    cells = pd.DataFrame({
        "analysis_unit_key": ["unit_pd1", "unit_pd1", "unit_ctla", "unit_ctla"],
        "cohort_id": ["GSE272993"] * 4,
        "study_subject_key": ["GSE272993::14-1189"] * 4,
        "normalized_timepoint": ["baseline"] * 4,
        "tissue_context": ["tumor"] * 4,
        "lesion_context": ["unknown"] * 4,
        "biological_sample_key": ["bio_pd1"] * 2 + ["bio_ctla"] * 2,
        "legacy_sample_key": ["legacy_pd1"] * 2 + ["legacy_ctla"] * 2,
        "library_key": ["lib_pd1"] * 2 + ["lib_ctla"] * 2,
        "treatment_arm": ["aPD1"] * 2 + ["aCTLA-4"] * 2,
        "assignment_status": ["resolved"] * 4,
        "allowed_phase4b_fraction": ["yes"] * 4,
        "harmonized_coarse_label": ["T_NK"] * 4,
        "harmonized_mid_label": ["CD8_T"] * 4,
        "harmonized_fine_label": ["Unknown"] * 4,
        "reliability_level": ["mid_reliable"] * 4,
    })
    cells.to_parquet(out4a / "handoff/cell_state_annotation_master.parquet/cohort_id=GSE272993.parquet")
    pd.DataFrame({
        "sample_key": ["legacy_pd1", "legacy_ctla"], "cohort_id": ["GSE272993"] * 2,
        "study_subject_key": ["GSE272993::14-1189"] * 2,
        "timepoint": ["baseline", "baseline"], "tissue_source": ["tumor", "tumor"],
        "treatment_context": ["aPD1", "aCTLA-4"],
        "response_raw": ["R", "NR"], "response_endpoint_type": ["RECIST", "RECIST"],
        "response_binary_harmonized": ["responder", "non_responder"],
        "response_harmonization_confidence": ["high", "high"],
        "supervised_use_allowed": ["yes", "yes"], "support_use_allowed": ["yes", "yes"],
        "dataset_role": ["conditional", "support"], "cancer_type": ["melanoma", "melanoma"],
    }).to_csv(metadata / "sample_metadata_master.effective_v1.csv", index=False)
    pd.DataFrame({
        "patient_key": ["legacy_patient"], "study_subject_key": ["GSE272993::14-1189"],
        "cancer_type": ["melanoma"], "split": ["test"],
        "response_binary_harmonized": ["unknown_or_mixed"],
    }).to_csv(metadata / "patient_metadata_master.effective_v1.csv", index=False)
    repair.build_phase4b(out4a, metadata, out4b)
    bound = pd.read_csv(out4b / "response_environment/sample_metadata_analysis_unit_v1.csv")
    assert bound.set_index("sample_key").loc["unit_pd1", "response_binary_harmonized"] == "responder"
    assert bound.set_index("sample_key").loc["unit_ctla", "response_binary_harmonized"] == "non_responder"
