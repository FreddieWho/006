"""Lightweight invariants for the metadata-only Stage 1 registry."""

from pathlib import Path
import csv
from collections import Counter
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "scripts/v7/registry/build_stage1_registry.py"


def rows(name: str, directory: Path) -> list[dict[str, str]]:
    with (directory / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Stage1RegistryTest(unittest.TestCase):
    def test_build_and_core_invariants(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            result = subprocess.run(
                [sys.executable, str(BUILDER), "--output-dir", str(output)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            logical = rows("multimodal_logical_units.tsv", output)
            physical = rows("spatial_physical_units.tsv", output)
            crosswalk = rows("patient_block_section_crosswalk.tsv", output)
            roles = rows("data_role_assignment.tsv", output)
            with (ROOT / "config/v7/data_source_registry.tsv").open(encoding="utf-8", newline="") as handle:
                configured = {r["source_id"] for r in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(len(logical), len({r["logical_unit_id"] for r in logical}))
            self.assertEqual(len(physical), len({r["physical_unit_id"] for r in physical}))
            self.assertEqual(
                {r["physical_unit_id"] for r in crosswalk},
                {r["physical_unit_id"] for r in physical},
            )
            self.assertEqual(len(roles), 8 * len({r["source_id"] for r in roles}))
            self.assertTrue({r["logical_unit_id"] for r in roles}.issubset({r["logical_unit_id"] for r in logical}))
            self.assertTrue(
                {r["source_id"] for r in physical}.issubset(configured),
                sorted({r["source_id"] for r in physical} - configured),
            )
            self.assertTrue(any(r["dataset_id"] == "GSE291246" and r["cancer"] == "basal_cell_carcinoma" for r in physical))
            self.assertTrue(any(r["dataset_id"] == "GSE238264" and r["response"] == "responder" for r in physical))
            gse291 = [r for r in physical if r["dataset_id"] == "GSE291246"]
            self.assertEqual(Counter(r["counts_layer"] for r in gse291), Counter({"h5ad_or_spot_matrix": 17, "manifest_only_no_expression_qc": 18}))
            self.assertTrue(all("GSE175540" not in r["source_path"] for r in physical if r["dataset_id"] == "GEO::GSE274557"))
            duplicate = rows("duplicate_lineage.tsv", output)
            self.assertTrue(any(r["relation_type"] == "same_patient_nested_spatial_sections" and "GSE291246" in r["member_id"] for r in duplicate))
            gse291_role = [r for r in roles if r["source_id"] == "006_GSE291246"]
            self.assertEqual({r["treatment_context"] for r in gse291_role}, {"post_PD1|treatment_naive"})

    def test_missing_core_input_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            result = subprocess.run(
                [sys.executable, str(BUILDER), "--output-dir", str(Path(temp) / "out"), "--spatial-root", str(Path(temp) / "missing")],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("STAGE1_BLOCKED_MISSING_INPUT", result.stderr)


if __name__ == "__main__":
    unittest.main()
