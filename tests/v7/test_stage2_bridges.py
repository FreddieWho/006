"""Synthetic contract tests for the Stage 2 bulk/perturbation bridges."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from src.v7.ontology.bridges import (
    BridgeContractError,
    build_gse207422_crosswalk,
    parse_gse193736_sample_id,
    read_gse193736_bulk_counts,
    read_gse207422_bulk_log2tpm,
    read_gse207422_metadata,
)


def _safe_metadata(rows: list[tuple[str, str, str]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows, columns=["sample_id", "patient_id", "resource"])
    frame["timepoint"] = frame["resource"].map(
        lambda value: "pre_treatment" if str(value).lower().startswith("pre") else "post_treatment"
    )
    return frame


def _gse193736_columns() -> list[str]:
    return [
        f"{lineage}_{perturbation}_{rest_stim}_{replicate}"
        for lineage in ("CD4", "CD8")
        for perturbation in ("LTBR", "tNGFR")
        for rest_stim in ("rest", "stim")
        for replicate in (1, 2, 3)
    ]


class GSE207422BridgeTest(unittest.TestCase):
    def test_metadata_reader_requests_only_safe_columns(self) -> None:
        source = pd.DataFrame(
            {
                "Sample": ["BD_immune05", "spreadsheet footer"],
                "Patient": ["P05", None],
                "Resource": ["Pre-treatment biopsy", None],
            }
        )
        with patch("src.v7.ontology.bridges.pd.read_excel", return_value=source) as reader:
            result = read_gse207422_metadata("metadata.xlsx")
        self.assertEqual(
            reader.call_args.kwargs["usecols"], ["Sample", "Patient", "Resource"]
        )
        self.assertEqual(
            result.columns.tolist(), ["sample_id", "patient_id", "resource", "timepoint"]
        )
        self.assertEqual(result.loc[0, "timepoint"], "pre_treatment")
        self.assertEqual(len(result), 1)

    def test_crosswalk_has_two_exact_pairs_and_one_support_only_pair(self) -> None:
        scrna = _safe_metadata(
            [
                ("BD_immune05", "P05", "Pre-treatment biopsy"),
                ("BD_immune07", "P07", "Post-treatment surgery"),
                ("BD_immune08", "P08", "Pre-treatment biopsy"),
                ("BD_immune09", "P09", "Post-treatment surgery"),
            ]
        )
        bulk = _safe_metadata(
            [
                ("R05", "P05", "Pre_biopsy"),
                ("R07", "P07", "Pre_biopsy"),
                ("R08", "P08", "Pre_biopsy"),
                ("R16", "P16", "Pre-treatment biopsy"),
            ]
        )
        crosswalk = build_gse207422_crosswalk(scrna, bulk).set_index("patient_id")
        self.assertEqual(crosswalk.index.tolist(), ["P05", "P07", "P08"])
        self.assertEqual(crosswalk.loc["P05", "match_status"], "exact_matched")
        self.assertEqual(crosswalk.loc["P08", "bridge_role"], "paired")
        self.assertEqual(crosswalk.loc["P07", "match_status"], "timepoint_mismatch")
        self.assertEqual(crosswalk.loc["P07", "bridge_role"], "support_only")
        self.assertFalse(any("response" in column.lower() for column in crosswalk.columns))

    def test_bulk_log2tpm_is_transposed_to_sample_by_gene(self) -> None:
        source = pd.DataFrame(
            {
                "Gene": ["CXCL13", "TCF7"],
                "sample_a": [1.0, 2.0],
                "sample_b": [3.0, 4.0],
            }
        )
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bulk.tsv.gz"
            source.to_csv(path, sep="\t", index=False)
            matrix = read_gse207422_bulk_log2tpm(path)
        self.assertEqual(matrix.index.tolist(), ["sample_a", "sample_b"])
        self.assertEqual(matrix.columns.tolist(), ["CXCL13", "TCF7"])
        self.assertEqual(matrix.index.name, "sample_id")
        self.assertEqual(matrix.columns.name, "gene")
        self.assertEqual(matrix.loc["sample_b", "TCF7"], 4.0)


class GSE193736BridgeTest(unittest.TestCase):
    def test_sample_id_parser(self) -> None:
        parsed = parse_gse193736_sample_id("CD8_tNGFR_stim_3")
        self.assertEqual(
            parsed,
            {"lineage": "CD8", "perturbation": "tNGFR", "rest_stim": "stim", "replicate": 3},
        )
        with self.assertRaises(BridgeContractError):
            parse_gse193736_sample_id("CD8_tNGFR_unknown_3")

    def test_bulk_counts_parse_complete_factorial_design(self) -> None:
        columns = _gse193736_columns()
        source = pd.DataFrame(
            {
                "Unnamed: 0": [1, 2],
                "gene_id": ["ENSG00000000003.14", "ENSG00000000419.12"],
                **{
                    column: np.array([offset, offset + 1], dtype=int)
                    for offset, column in enumerate(columns)
                },
            }
        )
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "counts.csv.gz"
            source.to_csv(path, index=False)
            matrix, metadata = read_gse193736_bulk_counts(path)
        self.assertEqual(matrix.shape, (24, 2))
        self.assertEqual(matrix.index.name, "sample_id")
        self.assertEqual(matrix.columns.name, "gene_id")
        self.assertEqual(
            metadata.columns.tolist(),
            ["sample_id", "lineage", "perturbation", "rest_stim", "replicate"],
        )
        self.assertEqual(set(metadata["lineage"]), {"CD4", "CD8"})
        self.assertEqual(set(metadata["perturbation"]), {"LTBR", "tNGFR"})
        self.assertEqual(set(metadata["rest_stim"]), {"rest", "stim"})
        self.assertEqual(set(metadata["replicate"]), {1, 2, 3})

    def test_bulk_counts_reject_unversioned_ensembl(self) -> None:
        columns = _gse193736_columns()
        source = pd.DataFrame(
            {"gene_id": ["ENSG00000000003"], **{column: [0] for column in columns}}
        )
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "counts.csv"
            source.to_csv(path, index=False)
            with self.assertRaisesRegex(BridgeContractError, "versioned Ensembl"):
                read_gse193736_bulk_counts(path)


if __name__ == "__main__":
    unittest.main()
