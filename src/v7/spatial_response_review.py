"""Post-treatment association only, for five previously defined spatial proxies."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .perturbation_review import condition_contrast
from .spatial_control import sha256_file


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    output=Path(args.output)
    output.mkdir(parents=True,exist_ok=False)
    raw=Path("data/raw/ST/GEO/GSE238264/GSE238264_series_matrix.txt.gz")
    metadata=Path("results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.GSE238264.tsv")
    primitive=Path("results/v7/spatial_discovery/niche_primitives.tsv")
    candidate=Path("results/v7/spatial_discovery/response_blind_shared_architectures.tsv")
    fields={}
    with gzip.open(raw,"rt") as stream:
        for line in stream:
            if line.startswith("!Sample_geo_accession"):
                fields["gsm"]=next(csv.reader([line],delimiter="\t"))[1:]
            if line.startswith("!Sample_characteristics_ch1") and 'phenotype:' in line:
                fields["response"]=[v.removeprefix("phenotype: ") for v in next(csv.reader([line],delimiter="\t"))[1:]]
    source=pd.DataFrame({"gsm":fields["gsm"],"source_response":fields["response"]})
    meta=pd.read_csv(metadata,sep="\t")
    meta["gsm"]=meta.source_file.str.extract(r"(GSM\d+)")[0]
    meta=meta.merge(source,on="gsm",validate="one_to_one")
    response_map={"responder":1,"non-responder":0,"non_responder":0}
    meta["y"]=meta.source_response.map(response_map)
    if meta.y.isna().any() or not meta.y.eq(meta.response_binary.map(response_map)).all():
        raise ValueError("GEO phenotype conflicts with inherited response")
    meta["patient_id"]="GSE238264::"+meta.patient_id
    proxies=pd.read_csv(primitive,sep="\t")
    proxies=proxies[proxies.dataset_id.eq("GSE238264")]
    locked=pd.read_csv(candidate,sep="\t")
    if set(proxies.primitive_id)!=set(locked.primitive_id):
        raise ValueError("candidate selection changed after clinical read")
    patient=proxies.groupby(["patient_id","primitive_id"],as_index=False).estimate.median()
    patient=patient.merge(meta[["patient_id","gsm","y"]],on="patient_id",validate="many_to_one")
    if patient.patient_id.nunique()!=len(meta):
        raise ValueError("spatial/response patient identity does not close")
    rows=[]
    for name,f in patient.groupby("primitive_id",sort=True):
        result=condition_contrast(f.loc[f.y.eq(1),"estimate"],f.loc[f.y.eq(0),"estimate"])
        result["independence"]="patient_4_responders_3_nonresponders"
        result["claim"]="post_treatment_proxy_association_not_prediction_or_longitudinal_repair"
        rows.append(dict(primitive_id=name,**result))
    pd.DataFrame(rows).to_csv(output/"post_treatment_spatial_response.tsv",sep="\t",index=False)
    patient.to_csv(output/"patient_identity_and_proxy.tsv",sep="\t",index=False)
    meta[["patient_id","gsm","source_response","response_binary","y"]].to_csv(output/"response_source_audit.tsv",sep="\t",index=False)
    receipt=dict(status="COMPUTATION_COMPLETE",n_patients=len(meta),n_candidates=len(rows),
                 exposure="five_predefined_candidates_frozen_before_clinical_read; development_support_only",
                 stage_closure="OPEN_BASELINE_LONGITUDINAL_AND_INDEPENDENT_VALIDATION",
                 endpoint="pathologic_response_in_inherited_manifest; GEO_phenotype_labels_verified",
                 multiple_testing="exploratory_five_candidates_no_significance_claim",
                 inputs={str(p):sha256_file(p) for p in [raw,metadata,primitive,candidate]},
                 code_sha256=sha256_file(__file__),
                 outputs={p.name:sha256_file(p) for p in output.iterdir() if p.is_file()})
    (output/"RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":
    main()
