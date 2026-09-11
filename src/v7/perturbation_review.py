"""Measured LTBR-labelled versus tNGFR-labelled contrasts, without repair claims."""
from itertools import combinations
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .spatial_control import sha256_file


def condition_contrast(treatment, control):
    """Exact two-sided sample-label null; replicate numbers do not imply pairs."""
    a, b = np.asarray(treatment,float), np.asarray(control,float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a)<2 or len(b)<2:
        return dict(status="NOT_TESTABLE_TOO_FEW_SAMPLES")
    x = np.concatenate([a,b])
    observed = float(a.mean()-b.mean())
    null = []
    for selected in combinations(range(len(x)),len(a)):
        mask = np.zeros(len(x),bool)
        mask[list(selected)] = True
        null.append(float(x[mask].mean()-x[~mask].mean()))
    p = float(np.mean(np.abs(null)>=abs(observed)-1e-12))
    return dict(status="MEASURED_CONDITION_CONTRAST",difference=observed,
                treatment_mean=float(a.mean()),control_mean=float(b.mean()),
                n_treatment=len(a),n_control=len(b),exact_sample_label_p=p,
                n_label_assignments=len(null),
                independence="sample_replicates_donor_pairing_not_established",
                claim="measured_condition_direction_not_PD1X_repair_or_drug_effect")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    output=Path(args.output)
    output.mkdir(parents=True,exist_ok=False)
    start=time.time()
    score_path=Path("results/v7/ontology/scores/expression_unit_scores.parquet")
    design_path=Path("results/v7/ontology/gse193736_perturbation_design.tsv")
    raw_path=Path("data/perturb_seq/GSE193736/GSE193736_bulkRNAseq_counts.csv.gz")
    scores=pd.read_parquet(score_path,filters=[("cohort_id","==","GSE193736")])
    design=pd.read_csv(design_path,sep="\t")
    if design.sample_id.duplicated().any() or scores.duplicated(["sample_key","feature_id"]).any():
        raise ValueError("duplicate sample/feature or design identity")
    if set(scores.sample_key)!=set(design.sample_id):
        raise ValueError("sample identity does not close")
    merged=scores.merge(design,left_on="sample_key",right_on="sample_id",validate="many_to_one")
    rows=[]
    for (lineage,stim,feature),frame in merged.groupby(["lineage","rest_stim","feature_id"],sort=True):
        result=condition_contrast(frame.loc[frame.perturbation.eq("LTBR"),"score_native"],
                                  frame.loc[frame.perturbation.eq("tNGFR"),"score_native"])
        rows.append(dict(lineage=lineage,rest_stim=stim,feature_id=feature,
                         min_coverage=frame.coverage.min(),**result))
    result=pd.DataFrame(rows)
    result.to_csv(output/"measured_perturbation_direction.tsv",sep="\t",index=False)
    transfer=result.pivot(index="feature_id",columns=["lineage","rest_stim"],values="difference")
    transfer.columns=["__".join(c) for c in transfer.columns]
    transfer.to_csv(output/"context_direction_matrix.tsv",sep="\t")
    receipt=dict(status="COMPUTATION_COMPLETE",n_samples=len(design),n_contrasts=len(result),
                 stage_closure="OPEN_PD1X_TARGET_LINK_EXTERNAL_VALIDATION",
                 multiple_testing="no_feature_significance_claim_or_selection_from_nominal_p",
                 raw_input_sha256=sha256_file(raw_path),score_input_sha256=sha256_file(score_path),
                 design_sha256=sha256_file(design_path),code_sha256=sha256_file(__file__),
                 elapsed_seconds=time.time()-start,
                 outputs={p.name:sha256_file(p) for p in output.iterdir() if p.is_file()})
    (output/"RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":
    main()
