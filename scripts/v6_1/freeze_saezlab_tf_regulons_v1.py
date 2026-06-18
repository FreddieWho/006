#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

URLS = {
    "collectri_human": "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&format=tsv&fields=sources,references,curation_effort,dorothea_level",
    "dorothea_human_ABC": "https://omnipathdb.org/interactions?datasets=dorothea&genesymbols=yes&format=tsv&fields=sources,references,curation_effort,dorothea_level&dorothea_levels=A,B,C",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_tsv(url: str) -> pd.DataFrame:
    req = urllib.request.Request(url, headers={"User-Agent": "step2.9-tf-regulon-freeze/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        data = resp.read()
    return pd.read_csv(io.BytesIO(data), sep="\t")


def clean_gene(x: object) -> str:
    return str(x).strip().upper()


def valid_symbol(x: str) -> bool:
    if not x or x == "NAN":
        return False
    if x.startswith("HSA-") or x.startswith("MIMAT") or x.startswith("COMPLEX:"):
        return False
    return bool(re.match(r"^[A-Z0-9][A-Z0-9._-]*$", x))


def mode(row: pd.Series) -> int:
    stim = str(row.get("consensus_stimulation", row.get("is_stimulation", ""))).lower() == "true"
    inhib = str(row.get("consensus_inhibition", row.get("is_inhibition", ""))).lower() == "true"
    if inhib and not stim:
        return -1
    return 1


def normalize(df: pd.DataFrame, resource: str) -> pd.DataFrame:
    out = pd.DataFrame({
        "tf": df["source_genesymbol"].map(clean_gene),
        "target": df["target_genesymbol"].map(clean_gene),
        "direction": df.apply(mode, axis=1),
        "source_uniprot": df.get("source", pd.Series([""] * len(df))).astype(str),
        "target_uniprot": df.get("target", pd.Series([""] * len(df))).astype(str),
        "sources": df.get("sources", pd.Series([""] * len(df))).fillna("").astype(str),
        "references": df.get("references", pd.Series([""] * len(df))).fillna("").astype(str),
        "dorothea_level": df.get("dorothea_level", pd.Series([""] * len(df))).fillna("").astype(str),
        "curation_effort": pd.to_numeric(df.get("curation_effort", pd.Series([pd.NA] * len(df))), errors="coerce"),
        "resource": resource,
    })
    out = out[out["tf"].map(valid_symbol) & out["target"].map(valid_symbol)]
    out = out[out["tf"] != out["target"]]
    out["signed_weight"] = out["direction"].astype(float)
    out["evidence_level"] = out["dorothea_level"].replace("", "not_provided")
    return out.drop_duplicates(["tf", "target", "direction", "resource", "evidence_level"]).reset_index(drop=True)


def to_step_registry(consensus: pd.DataFrame) -> pd.DataFrame:
    grouped = consensus.groupby(["tf", "target"], dropna=False).agg(
        signed_weight=("signed_weight", "mean"),
        direction=("direction", lambda x: -1 if x.mean() < 0 else 1),
        resource=("resource", lambda x: ";".join(sorted(set(map(str, x))))),
        source_resource=("sources", lambda x: ";".join(sorted({v for item in x for v in str(item).split(";") if v}))),
        evidence_level=("evidence_level", lambda x: ";".join(sorted(set(map(str, x))))),
        references=("references", lambda x: ";".join(sorted({v for item in x for v in str(item).split(";") if v}))[:32000]),
        curation_effort=("curation_effort", "max"),
    ).reset_index()
    grouped["direction"] = grouped["direction"].astype(int)
    grouped["feature_id"] = "tf_activity__" + grouped["tf"].str.replace(r"[^0-9A-Za-z]+", "_", regex=True).str.strip("_")
    grouped["feature_family"] = "tf_activity"
    grouped["feature_name"] = grouped["tf"]
    grouped["gene_symbol"] = grouped["target"]
    grouped["used_in_main"] = True
    grouped["cell_context"] = "sample;cellstate"
    grouped["notes"] = "saezlab_omnipath_collectri_plus_dorothea_ABC_consensus"
    return grouped[[
        "feature_id", "feature_family", "feature_name", "gene_symbol", "direction", "signed_weight",
        "resource", "source_resource", "evidence_level", "references", "curation_effort", "used_in_main", "cell_context", "notes"
    ]].sort_values(["feature_id", "gene_symbol", "direction"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default="data/tf_regulon/consensus_tf_regulon_v1")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()
    normalized = []
    raw_stats = {}
    for name, url in URLS.items():
        raw = download_tsv(url)
        raw_path = out_dir / f"{name}_raw.csv.gz"
        raw.to_csv(raw_path, index=False)
        norm = normalize(raw, name)
        norm_path = out_dir / f"{name}_normalized.csv.gz"
        norm.to_csv(norm_path, index=False)
        normalized.append(norm)
        raw_stats[name] = {"url": url, "raw_rows": int(len(raw)), "normalized_rows": int(len(norm)), "raw_columns": list(raw.columns)}
    consensus = pd.concat(normalized, ignore_index=True)
    registry = to_step_registry(consensus)
    registry_path = out_dir / "saezlab_tf_regulon_consensus_v1.csv.gz"
    registry.to_csv(registry_path, index=False)
    files = sorted(out_dir.glob("*.csv.gz"))
    manifest = {
        "created_at": created_at,
        "resource_family": "saezlab_tf_regulon",
        "retrieval_method": "direct OmniPath API; saezlab/OmniPath resources",
        "resources_retrieved": list(URLS),
        "raw_stats": raw_stats,
        "consensus_rows": int(len(registry)),
        "consensus_tfs": int(registry["feature_name"].nunique()),
        "consensus_targets": int(registry["gene_symbol"].nunique()),
        "policy": {"trrust_included": False, "scenic_included": False, "dorothea_levels": "A/B/C only", "response_blind": True},
        "sha256": {p.name: sha256_file(p) for p in files},
    }
    (out_dir / "tf_regulon_resource_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# SaezLab TF Regulon Resource Validation Report", "",
        f"- created_at: `{created_at}`",
        f"- resources_retrieved: `{';'.join(URLS)}`",
        f"- consensus rows: `{len(registry)}`",
        f"- consensus TFs: `{registry['feature_name'].nunique()}`",
        f"- consensus targets: `{registry['gene_symbol'].nunique()}`",
        "- TRRUST included: `false`",
        "- SCENIC included: `false`",
        "- source: `OmniPath API / saezlab resources`", "",
        "## Raw Stats", "",
    ]
    for name, stat in raw_stats.items():
        report.append(f"- {name}: raw_rows=`{stat['raw_rows']}`, normalized_rows=`{stat['normalized_rows']}`")
    report += ["", "## Normalized Schema", "", "- `feature_id`, `feature_family`, `feature_name`, `gene_symbol`, `direction`, `signed_weight`, `resource`, `source_resource`, `evidence_level`, `references`, `curation_effort`, `used_in_main`, `cell_context`, `notes`"]
    (out_dir / "tf_regulon_resource_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {registry_path}")
    print(f"Consensus rows: {len(registry)}")
    print(f"Consensus TFs: {registry['feature_name'].nunique()}")


if __name__ == "__main__":
    main()
