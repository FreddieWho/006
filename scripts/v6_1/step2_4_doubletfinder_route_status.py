from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    out = root / "results" / "v6_1" / "step2" / run_id / "03_qc"
    out.mkdir(parents=True, exist_ok=True)
    route_csv = out / "doubletfinder_route_status.csv"
    route_md = out / "doubletfinder_route_report.md"

    def r_check(pkg: str) -> bool:
        cmd = ["Rscript", "-e", f"cat(requireNamespace('{pkg}', quietly=TRUE), '\\n')"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        return p.returncode == 0 and "TRUE" in p.stdout

    has_df = r_check("DoubletFinder")
    has_seu = r_check("Seurat")
    has_sce = r_check("SingleCellExperiment")

    status = "ready_for_pilot" if (has_df and has_seu) else "not_run"
    reason = ""
    if status != "ready_for_pilot":
        reason = "missing_required_R_packages"

    prior_cap = 0.05
    route_csv.write_text(
        "run_id,created_at,route,status,has_doubletfinder,has_seurat,has_singlecellexperiment,expected_doublet_rate_cap,not_run_reason\n"
        f"{run_id},{now_iso()},DoubletFinder,{status},{has_df},{has_seu},{has_sce},{prior_cap},{reason}\n",
        encoding="utf-8",
    )
    route_md.write_text(
        "\n".join(
            [
                "# DoubletFinder Route Status",
                "",
                f"- run_id: {run_id}",
                f"- created_at: {now_iso()}",
                f"- has_DoubletFinder: {has_df}",
                f"- has_Seurat: {has_seu}",
                f"- has_SingleCellExperiment: {has_sce}",
                f"- expected_doublet_rate_cap: {prior_cap}",
                f"- status: {status}",
                f"- not_run_reason: {reason or 'none'}",
                "",
                "Note: this run records route readiness/status only. Full DoubletFinder execution is a dedicated branch.",
            ]
        ),
        encoding="utf-8",
    )
    return {"run_id": run_id, "status": status, "csv": str(route_csv), "md": str(route_md)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(main(args.run_id, args.root_dir))
