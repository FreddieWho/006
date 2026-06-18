#!/usr/bin/env python3
"""Finish Phase4B rescue-expression activity patch from existing parquet."""
from __future__ import annotations

import pandas as pd

from run_phase4b_rescue_expression_patch import (
    ACTIVITY,
    PSEUDOBULK,
    append_audit_tables,
    feature_row,
    refresh_output_index,
    score_signatures,
    update_feature_dictionaries,
    update_manifests,
    update_report,
    update_support_matrix,
    write_csv,
)


def main() -> None:
    pb = pd.read_parquet(PSEUDOBULK / "pseudobulk_matrix_normalized_expression.parquet")
    source_audit = pd.read_csv(PSEUDOBULK / "pseudobulk_rescue_expression_source_audit.csv", dtype=str, keep_default_na=False)
    pb_dict = pd.DataFrame(
        [
            feature_row(
                feature_id="pseudobulk_normalized_expression__rescue_addendum_gene_matrix",
                family="pseudobulk_normalized_expression",
                source_file="pseudobulk/pseudobulk_matrix_normalized_expression.parquet",
                source_layer="normalized_expression_rescue_addendum",
                method="local_addendum_pseudobulk_standardization",
                gene_set="",
                allowed="support_or_sensitivity_expression_only",
                notes="normalized/support pseudobulk only; not raw count and not cell-state-specific",
            )
        ]
    )
    sig_matrix, sig_cov, sig_dict = score_signatures(pb)
    tf_matrix = pb[["sample_key", "cohort_id"]].copy()
    tf_cov = pd.DataFrame(
        [
            {
                "feature_id": "tf_activity",
                "resource": "saezlab_tf_regulon_consensus_v1",
                "n_genes_total": 0,
                "n_genes_present": 0,
                "coverage": 0.0,
                "status": "skipped_lightweight_patch",
            }
        ]
    )
    tf_dict = pd.DataFrame()

    write_csv(sig_matrix, ACTIVITY / "signature_activity_matrix.sample_level.csv")
    write_csv(tf_matrix, ACTIVITY / "tf_activity_matrix.sample_level.csv")
    write_csv(pd.DataFrame(columns=["sample_key", "cohort_id"]), ACTIVITY / "pathway_activity_matrix.sample_level.csv")
    write_csv(pd.concat([sig_cov, tf_cov], ignore_index=True, sort=False), ACTIVITY / "gene_set_coverage_report.csv")

    activity = sig_matrix.merge(tf_matrix, on=["sample_key", "cohort_id"], how="outer")
    update_support_matrix(activity)
    update_feature_dictionaries(pb_dict, sig_dict, tf_dict)
    append_audit_tables(activity)
    update_manifests(pb, sig_matrix, tf_matrix)
    update_report(pb, sig_matrix, tf_matrix, source_audit)
    refresh_output_index()
    print(f"activity patch complete: samples={len(pb)} signatures={max(sig_matrix.shape[1] - 2, 0)}")


if __name__ == "__main__":
    main()
