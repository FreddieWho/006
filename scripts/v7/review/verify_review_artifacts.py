#!/usr/bin/env python3
"""Verify completed review receipts and preserve their exact source versions."""
import hashlib
import io
import json
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path("results/v7")
RUNS=[
    "clinical_anchor/review_rerun_2026-09-11",
    "clinical_anchor/post_spatial_review_2026-09-11",
    "spatial_discovery/review_factors_2026-09-11",
    "spatial_discovery/review_permutation_matrix_2026-09-11",
    "repair/review_rerun_2026-09-11",
    "repair/direction_review_2026-09-11",
    "perturbation/review_rerun_2026-09-11",
]
ARCHIVES=[ROOT/"spatial_discovery"/d/"source_at_start.tar.gz" for d in
          ["review_permutation_2026-09-11","review_permutation_fast_2026-09-11"]]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    historical={}
    for archive in ARCHIVES:
        with tarfile.open(archive) as tar:
            for item in tar:
                if item.isfile():
                    content=tar.extractfile(item).read()
                    historical[(item.name,digest(content))]=content
    checks=[]
    for run in RUNS:
        path=ROOT/run
        receipt=json.loads((path/"RUN_RECEIPT.json").read_text())
        if receipt["status"]!="COMPUTATION_COMPLETE":
            raise ValueError(f"unfinished run: {run}")
        for name,expected in receipt["outputs"].items():
            p=Path(name)
            p=p if p.is_file() else path/p
            if not p.is_file() or digest(p.read_bytes())!=expected:
                raise ValueError(f"output hash mismatch: {p}")
        code=receipt.get("code_hashes",{})
        if not code:
            name=("spatial_response_review" if "post_spatial" in run else
                  "direction_review" if "direction_review" in run else "perturbation_review")
            code={f"src/v7/{name}.py":receipt["code_sha256"]}
        archive_path=path/"REPLAY_SOURCE.tar.gz"
        if code and not archive_path.exists():
            with tarfile.open(archive_path,"w:gz") as tar:
                for name,expected in code.items():
                    content=Path(name).read_bytes()
                    if digest(content)!=expected:
                        content=historical.get((name,expected))
                    if content is None or digest(content)!=expected:
                        raise ValueError(f"missing execution-time code: {name}")
                    item=tarfile.TarInfo(name);item.size=len(content);tar.addfile(item,io.BytesIO(content))
        if code:
            with tarfile.open(archive_path) as tar:
                for name,expected in code.items():
                    if digest(tar.extractfile(name).read())!=expected:
                        raise ValueError(f"source snapshot mismatch: {name}")
        checks.append(dict(run=run,receipt_sha256=digest((path/"RUN_RECEIPT.json").read_bytes()),
                           outputs_verified=len(receipt["outputs"]),execution_source_files_verified=len(code),
                           source_snapshot=str(archive_path) if code else "current_source_digest_in_receipt",
                           source_snapshot_sha256=digest(archive_path.read_bytes()) if code else None,
                           five_documents_present=all((path/f).is_file() for f in
                               ["PHASE_SUMMARY.md","OUTPUT_MANIFEST.yaml","DECISION_LOG.md","EVIDENCE_AND_CONFLICTS.md","NEXT_PHASE_READINESS.yaml"])))
    tests=ET.parse(ROOT/"review_2026-09-11/targeted_tests.xml").getroot().find("testsuite").attrib
    result=dict(status="ARTIFACT_CHECKS_PASS",scientific_stage_closure="OPEN",runs=checks,
                targeted_tests=tests,independent_third_party_review="NOT_RUN_FOR_NEW_IMPLEMENTATION")
    (ROOT/"review_2026-09-11/artifact_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(status=result["status"],runs=len(checks),
                          output_hashes=sum(r["outputs_verified"] for r in checks),
                          execution_source_files=sum(r["execution_source_files_verified"] for r in checks),
                          all_five_documents=all(r["five_documents_present"] for r in checks),tests=tests["tests"])))


if __name__=="__main__":
    main()
