"""Cross-fitted molecular direction sensitivity; does not certify a repair barrier."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from .spatial_control import sha256_file


def fit_direction(pre, labels, train):
    """Training baseline only. No target-timepoint observations enter fitting."""
    x=pre[train]
    mean=np.mean(x,axis=0)
    scale=np.std(x,axis=0)
    scale[scale<1e-8]=1.
    z=(x-mean)/scale
    d=z[labels[train]==1].mean(axis=0)-z[labels[train]==0].mean(axis=0)
    norm=np.linalg.norm(d)
    return (d/norm if norm>1e-8 else np.zeros_like(d)),scale


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    output=Path(args.output);output.mkdir(parents=True,exist_ok=False)
    source=Path("results/v7/repair/review_rerun_2026-09-11/paired_displacements.tsv")
    reliability=Path("results/v7/ontology/measurement_reliability.tsv")
    paired=pd.read_csv(source,sep="\t")
    rel=pd.read_csv(reliability,sep="\t")
    features=sorted(rel.loc[rel.D2_class.eq("measurable"),"feature_id"].astype(str))
    features=[f for f in features if f in set(paired.feature_id)]
    records,summary,directions=[],[],[]
    for env,frame in paired[paired.feature_id.isin(features)].groupby("environment",sort=True):
        baseline=frame.drop_duplicates(["patient_key","feature_id"]).pivot(index="patient_key",columns="feature_id",values="score_standardized_pre")[features]
        if not np.isfinite(baseline.to_numpy()).all():
            raise ValueError("direction requires complete primary measurements; do not silently impute")
        y=frame.drop_duplicates("patient_key").set_index("patient_key").response_binary.reindex(baseline.index).to_numpy(int)
        pre=baseline.to_numpy(float)
        folds=list(StratifiedKFold(min(5,int(y.sum()),int((1-y).sum())),shuffle=True,random_state=20260911).split(pre,y))
        for target,sub in frame.groupby("target_timepoint",sort=True):
            post=sub.pivot(index="patient_key",columns="feature_id",values="score_standardized_post").reindex(index=baseline.index,columns=features).to_numpy(float)
            valid=np.isfinite(post).all(axis=1)
            def project(labels,record=False):
                projection=np.full(len(y),np.nan)
                for fold,(train,test) in enumerate(folds):
                    direction,scale=fit_direction(pre,labels,train)
                    projection[test]=((post[test]-pre[test])/scale) @ direction
                    if record:
                        for feature,value in zip(features,direction):
                            directions.append(dict(environment=env,target_timepoint=target,fold=fold,feature_id=feature,weight=value,n_training_patients=len(train)))
                        for i in test:
                            if valid[i]:
                                records.append(dict(environment=env,target_timepoint=target,patient_key=baseline.index[i],fold=fold,response_binary=int(y[i]),alignment=projection[i]))
                return projection
            observed=project(y,True)
            contrast=float(np.median(observed[valid & (y==1)])-np.median(observed[valid & (y==0)]))
            rng=np.random.default_rng(20260911)
            null=[]
            for _ in range(199):
                labels=y.copy()
                # Keep fold label counts and refit directions on every shuffle.
                for _,test in folds:
                    labels[test]=rng.permutation(labels[test])
                projected=project(labels)
                null.append(float(np.median(projected[valid & (labels==1)])-np.median(projected[valid & (labels==0)])))
            summary.append(dict(environment=env,target_timepoint=target,n_patients=int(valid.sum()),n_features=len(features),
                                responder_minus_nonresponder_alignment=contrast,
                                permutation_p=float((1+np.sum(np.abs(null)>=abs(contrast)))/200),
                                null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)),
                                result="EXPLORATORY_CROSSFITTED_DIRECTION_NOT_VALIDATED_REPAIR",
                                claim_boundary="baseline_discrimination_weak_no_external_or_spatial_validation"))
    pd.DataFrame(records).to_csv(output/"patient_direction_alignment.tsv",sep="\t",index=False)
    pd.DataFrame(directions).to_csv(output/"training_only_direction_weights.tsv",sep="\t",index=False)
    pd.DataFrame(summary).to_csv(output/"direction_permutation_results.tsv",sep="\t",index=False)
    receipt=dict(status="COMPUTATION_COMPLETE",stage_closure="OPEN_BARRIER_QUALIFICATION_AND_EXTERNAL_VALIDATION",features=features,
                 selection="all_frozen_D2_measurable_features_no_post_or_validation_label_feature_selection",
                 independent_review="NOT_RUN_NEW_IMPLEMENTATION_SELF_TESTED",
                 inputs={str(p):sha256_file(p) for p in [source,reliability]},code_sha256=sha256_file(__file__),
                 outputs={p.name:sha256_file(p) for p in output.iterdir() if p.is_file()})
    (output/"RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":
    main()
