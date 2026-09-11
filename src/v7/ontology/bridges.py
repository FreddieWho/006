"""Response-blind Stage 2 bridges for paired bulk and perturbation data.

The functions in this module only parse source data and establish technical
cross-modal identity.  They do not score biology or inspect clinical outcomes.
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Mapping

import numpy as np
import pandas as pd


GSE207422_METADATA_COLUMNS = ("Sample", "Patient", "Resource")
GSE207422_SAFE_COLUMNS = ("sample_id", "patient_id", "resource", "timepoint")
GSE193736_SAMPLE_PATTERN = re.compile(
    r"^(?P<lineage>CD4|CD8)_(?P<perturbation>[^_]+)_"
    r"(?P<rest_stim>rest|stim)_(?P<replicate>[1-9][0-9]*)$"
)
VERSIONED_ENSEMBL_PATTERN = re.compile(r"^ENSG[0-9]+\.[0-9]+$")


class BridgeContractError(ValueError):
    """Raised when a bridge input violates its dataset-specific contract."""


def _require_columns(frame: pd.DataFrame, required: tuple[str, ...], label: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise BridgeContractError(f"{label} is missing required columns: {missing}")


def _clean_identifier(series: pd.Series, label: str) -> pd.Series:
    clean = series.astype("string").str.strip()
    if clean.isna().any() or clean.eq("").any():
        raise BridgeContractError(f"{label} contains a missing or empty identifier")
    return clean.astype(str)


def _resource_to_timepoint(resource: object) -> str:
    value = str(resource).strip().lower().replace("_", "-")
    if value.startswith("pre"):
        return "pre_treatment"
    if value.startswith("post"):
        return "post_treatment"
    return "unknown"


def read_gse207422_metadata(path: str | Path) -> pd.DataFrame:
    """Read only the three response-blind GSE207422 metadata columns.

    The returned table deliberately cannot carry clinical outcome fields.
    Blank spreadsheet rows are discarded before identifiers are validated.
    """

    source = Path(path)
    raw = pd.read_excel(
        source,
        usecols=list(GSE207422_METADATA_COLUMNS),
        dtype={column: "string" for column in GSE207422_METADATA_COLUMNS},
    )
    _require_columns(raw, GSE207422_METADATA_COLUMNS, source.name)
    # The workbooks append prose footnotes below a blank row.  Those rows have
    # no Patient or Resource and are not sample metadata, so exclude them before
    # any values are propagated into the safe table.
    raw = (
        raw.loc[:, list(GSE207422_METADATA_COLUMNS)]
        .dropna(subset=list(GSE207422_METADATA_COLUMNS), how="any")
        .copy()
    )
    if raw.empty:
        raise BridgeContractError(f"{source.name} contains no metadata rows")

    result = raw.rename(
        columns={"Sample": "sample_id", "Patient": "patient_id", "Resource": "resource"}
    )
    result["sample_id"] = _clean_identifier(result["sample_id"], "sample_id")
    result["patient_id"] = _clean_identifier(result["patient_id"], "patient_id")
    result["resource"] = _clean_identifier(result["resource"], "resource")
    if result["sample_id"].duplicated().any():
        duplicates = result.loc[result["sample_id"].duplicated(False), "sample_id"].unique()
        raise BridgeContractError(f"duplicate sample_id values: {duplicates.tolist()}")
    result["timepoint"] = result["resource"].map(_resource_to_timepoint)
    return result.loc[:, list(GSE207422_SAFE_COLUMNS)].reset_index(drop=True)


def build_gse207422_crosswalk(
    scrna_metadata: pd.DataFrame,
    bulk_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Join scRNA and bulk samples by patient and classify timepoint agreement."""

    required = GSE207422_SAFE_COLUMNS
    _require_columns(scrna_metadata, required, "scRNA metadata")
    _require_columns(bulk_metadata, required, "bulk metadata")

    sc = scrna_metadata.loc[:, list(required)].rename(
        columns={
            "sample_id": "sc_sample_id",
            "resource": "sc_resource",
            "timepoint": "sc_timepoint",
        }
    )
    bulk = bulk_metadata.loc[:, list(required)].rename(
        columns={
            "sample_id": "bulk_sample_id",
            "resource": "bulk_resource",
            "timepoint": "bulk_timepoint",
        }
    )
    for label, frame, sample_column in (
        ("scRNA", sc, "sc_sample_id"),
        ("bulk", bulk, "bulk_sample_id"),
    ):
        if frame[sample_column].duplicated().any():
            raise BridgeContractError(f"{label} metadata contains duplicate sample IDs")
        if frame["patient_id"].duplicated().any():
            raise BridgeContractError(
                f"{label} metadata has multiple samples for one patient; "
                "a unique patient-timepoint mapping is required"
            )

    joined = sc.merge(bulk, on="patient_id", how="inner", validate="one_to_one")
    if joined.empty:
        raise BridgeContractError("GSE207422 scRNA and bulk metadata have no shared patients")
    exact = (
        joined["sc_timepoint"].eq(joined["bulk_timepoint"])
        & joined["sc_timepoint"].ne("unknown")
    )
    joined["match_status"] = np.where(exact, "exact_matched", "timepoint_mismatch")
    joined["bridge_role"] = np.where(exact, "paired", "support_only")
    columns = [
        "patient_id",
        "sc_sample_id",
        "bulk_sample_id",
        "sc_resource",
        "bulk_resource",
        "sc_timepoint",
        "bulk_timepoint",
        "match_status",
        "bridge_role",
    ]
    return joined.loc[:, columns].sort_values(
        ["patient_id", "sc_sample_id", "bulk_sample_id"], kind="stable"
    ).reset_index(drop=True)


def read_gse207422_bulk_log2tpm(path: str | Path) -> pd.DataFrame:
    """Read the GSE207422 gene-by-sample table as a sample-by-gene matrix."""

    source = Path(path)
    frame = pd.read_csv(source, sep="\t")
    if "Gene" not in frame.columns:
        raise BridgeContractError(f"{source.name} is missing the Gene column")
    genes = _clean_identifier(frame["Gene"], "Gene")
    if genes.duplicated().any():
        raise BridgeContractError("GSE207422 bulk matrix contains duplicate gene identifiers")
    sample_columns = [column for column in frame.columns if column != "Gene"]
    if not sample_columns:
        raise BridgeContractError("GSE207422 bulk matrix has no sample columns")
    if len(sample_columns) != len(set(map(str, sample_columns))):
        raise BridgeContractError("GSE207422 bulk matrix contains duplicate sample identifiers")

    numeric = frame.loc[:, sample_columns].apply(pd.to_numeric, errors="coerce")
    values = numeric.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise BridgeContractError("GSE207422 bulk log2TPM contains non-numeric or non-finite values")
    if (values < 0).any():
        raise BridgeContractError("GSE207422 bulk log2TPM contains negative values")

    numeric.index = genes
    matrix = numeric.T
    matrix.index = pd.Index([str(value).strip() for value in matrix.index], name="sample_id")
    matrix.columns = pd.Index(genes, name="gene")
    if (matrix.index == "").any():
        raise BridgeContractError("GSE207422 bulk matrix contains an empty sample identifier")
    return matrix


def parse_gse193736_sample_id(sample_id: str) -> Mapping[str, str | int]:
    """Parse one ``CD4_LTBR_rest_1``-style GSE193736 sample identifier."""

    match = GSE193736_SAMPLE_PATTERN.fullmatch(str(sample_id).strip())
    if match is None:
        raise BridgeContractError(f"invalid GSE193736 sample identifier: {sample_id!r}")
    fields: dict[str, str | int] = match.groupdict()
    fields["replicate"] = int(fields["replicate"])
    return fields


def read_gse193736_bulk_counts(
    path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read GSE193736 raw counts and its complete 24-sample factorial design.

    Returns a sample-by-versioned-Ensembl count matrix and a sample metadata
    table with ``lineage``, ``perturbation``, ``rest_stim`` and ``replicate``.
    """

    source = Path(path)
    frame = pd.read_csv(source)
    if "gene_id" not in frame.columns:
        raise BridgeContractError(f"{source.name} is missing the gene_id column")
    ignored = {"gene_id", "Unnamed: 0"}
    sample_columns = [column for column in frame.columns if column not in ignored]
    if len(sample_columns) != 24:
        raise BridgeContractError(
            f"GSE193736 bulk counts must contain 24 sample columns, found {len(sample_columns)}"
        )
    if len(sample_columns) != len(set(sample_columns)):
        raise BridgeContractError("GSE193736 bulk counts contain duplicate sample identifiers")

    sample_rows = []
    for sample_id in sample_columns:
        parsed = parse_gse193736_sample_id(sample_id)
        sample_rows.append({"sample_id": sample_id, **parsed})
    sample_metadata = pd.DataFrame.from_records(
        sample_rows,
        columns=["sample_id", "lineage", "perturbation", "rest_stim", "replicate"],
    )
    expected = {
        (lineage, perturbation, rest_stim, replicate)
        for lineage in ("CD4", "CD8")
        for perturbation in ("LTBR", "tNGFR")
        for rest_stim in ("rest", "stim")
        for replicate in (1, 2, 3)
    }
    observed = set(
        sample_metadata.loc[:, ["lineage", "perturbation", "rest_stim", "replicate"]]
        .itertuples(index=False, name=None)
    )
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise BridgeContractError(
            f"GSE193736 sample design is incomplete or unexpected; missing={missing}, extra={extra}"
        )

    genes = _clean_identifier(frame["gene_id"], "gene_id")
    invalid_genes = genes.loc[~genes.str.fullmatch(VERSIONED_ENSEMBL_PATTERN)]
    if not invalid_genes.empty:
        raise BridgeContractError(
            "GSE193736 gene_id must be a versioned Ensembl identifier; "
            f"examples={invalid_genes.head(5).tolist()}"
        )
    if genes.duplicated().any():
        raise BridgeContractError("GSE193736 bulk counts contain duplicate gene_id values")

    numeric = frame.loc[:, sample_columns].apply(pd.to_numeric, errors="coerce")
    values = numeric.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise BridgeContractError("GSE193736 bulk counts contain non-numeric or non-finite values")
    if (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise BridgeContractError("GSE193736 bulk counts must be non-negative integers")

    numeric.index = genes
    matrix = numeric.astype(np.int64).T
    matrix.index = pd.Index(sample_columns, name="sample_id")
    matrix.columns = pd.Index(genes, name="gene_id")
    return matrix, sample_metadata
