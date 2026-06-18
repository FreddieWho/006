from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "v6_1" / "step1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_COLLECTION_PATH = ROOT / "docs" / "data_collection.csv"
MVP_PATIENT_METADATA_PATH = ROOT / "mvp" / "outputs" / "patient_metadata_master.csv"
MVP_SIDECARS_DIR = ROOT / "mvp" / "sidecars"
MANIFEST_OUT = ROOT / "docs" / "v6_1" / "output_manifest_v6_1.yaml"
LOCAL_SOURCE_RECORDS = OUT_DIR / "source_records_local.tsv"


# ---------------------------------------------------------------------------
# Input asset discovery (replaces h5ad-only cohort_to_h5ad_path)
# ---------------------------------------------------------------------------
class InputAsset:
    """Discovered local input asset for a cohort."""

    def __init__(
        self,
        cohort_id: str,
        asset_type: str,
        path: Path,
        subset_id: Optional[str] = None,
        priority: int = 1,
    ):
        self.cohort_id = cohort_id
        self.asset_type = asset_type
        self.path = path
        self.subset_id = subset_id or cohort_id
        self.priority = priority

    def __repr__(self) -> str:
        return f"InputAsset({self.cohort_id}, {self.asset_type}, {self.path})"


def cohort_to_input_assets(cohort_id: str) -> List[InputAsset]:
    """Discover all local input assets for a cohort."""
    cid = cohort_id.lower()
    assets: List[InputAsset] = []

    # 1. canonical h5ad
    h5ad_path = ROOT / "data" / "processed" / "srt" / "raw" / f"{cid}.h5ad"
    if h5ad_path.exists():
        assets.append(InputAsset(cohort_id, "h5ad", h5ad_path))

    has_h5ad = any(a.asset_type == "h5ad" for a in assets)

    # 2. mtx10x (single-sample) — skip if canonical h5ad already present
    if not has_h5ad:
        mtx_candidates = list((ROOT / "data" / "imm" / cohort_id).glob("*.mtx.gz"))
        if not mtx_candidates:
            mtx_candidates = list((ROOT / "data" / "combo" / cohort_id).glob("*.mtx.gz"))
        for mtx in mtx_candidates:
            assets.append(InputAsset(cohort_id, "mtx10x", mtx.parent, priority=2))
            break  # one mtx10x asset per cohort (directory-based)

    # 3. txt_counts (dense gene-major) — skip if canonical h5ad already present
    if not has_h5ad:
        txt_candidates = (
            list((ROOT / "data" / "imm" / cohort_id).glob("*counts*.txt.gz"))
            + list((ROOT / "data" / "combo" / cohort_id).glob("*counts*.txt.gz"))
            + list((ROOT / "data" / "imm" / cohort_id).glob("*matrix*.txt.gz"))
            + list((ROOT / "data" / "combo" / cohort_id).glob("*matrix*.txt.gz"))
        )
        for txt in txt_candidates:
            assets.append(InputAsset(cohort_id, "txt_counts", txt, priority=2))
            break

    # 4. h5ad_gz — skip if canonical h5ad already present
    if not has_h5ad:
        h5ad_gz_candidates = (
            list((ROOT / "data" / "imm" / cohort_id).glob("*.h5ad.gz"))
            + list((ROOT / "data" / "combo" / cohort_id).glob("*.h5ad.gz"))
        )
        for gz in h5ad_gz_candidates:
            assets.append(InputAsset(cohort_id, "h5ad_gz", gz, priority=2))
            break

    # 5. bulk matrix
    bulk_matrix_candidates = (
        list((ROOT / "data" / "combo" / cohort_id).glob("*bulk*counts*.txt.gz"))
        + list((ROOT / "data" / "combo" / cohort_id).glob("*bulk*tpm*.txt.gz"))
    )
    for bm in bulk_matrix_candidates:
        assets.append(InputAsset(cohort_id, "bulk_matrix", bm, priority=2))
        break

    # 6. bulk metadata
    bulk_meta_candidates = (
        list((ROOT / "data" / "combo" / cohort_id).glob("*bulk*.xlsx"))
        + list((ROOT / "data" / "combo" / cohort_id).glob("*bulk*.csv"))
        + list((ROOT / "data" / "combo" / cohort_id).glob("*bulk*.txt.gz"))
    )
    for bm in bulk_meta_candidates:
        if bm not in [a.path for a in assets]:
            assets.append(InputAsset(cohort_id, "bulk_metadata", bm, priority=3))
            break

    return assets


# ---------------------------------------------------------------------------
# Non-h5ad sample metadata readers
# ---------------------------------------------------------------------------
def read_mtx10x_sample_metadata(cohort_id: str, mtx_dir: Path) -> pd.DataFrame:
    """Extract sample-level metadata from a 10x mtx directory.

    Returns one row per barcode (cell) since 10x barcodes are the finest
    granularity available without external metadata.
    """
    barcodes_path = None
    for pattern in ["*barcodes*.tsv.gz", "*barcodes*.tsv"]:
        candidates = list(mtx_dir.glob(pattern))
        if candidates:
            barcodes_path = candidates[0]
            break
    if barcodes_path is None:
        return pd.DataFrame()

    barcodes = pd.read_csv(barcodes_path, header=None, names=["barcode"])
    rows: List[Dict[str, str]] = []
    for bc in barcodes["barcode"]:
        rows.append(
            {
                "cohort_id": cohort_id,
                "sample_id": normalize_identifier(bc),
                "patient_id": "unknown",
                "timepoint_raw": "",
                "treatment_regimen_raw": "",
                "response_raw": "",
                "response_standard": "",
                "tissue_source": "",
                "source_path": str(barcodes_path),
                "ambiguity_note": "mtx10x_no_external_metadata",
            }
        )
    return pd.DataFrame(rows)


def read_txt_counts_sample_metadata(cohort_id: str, counts_path: Path) -> pd.DataFrame:
    """Extract sample-level metadata from a dense txt counts header.

    Reads the first line to obtain cell/sample ids.
    """
    with gzip.open(counts_path, "rt") as fh:
        header = fh.readline().strip().split("\t")
    # header may or may not have a gene-name placeholder
    first_data = fh.readline().strip().split("\t")
    n_cols = len(first_data) - 1
    if len(header) == n_cols + 1:
        sample_ids = header[1:]
    else:
        sample_ids = header

    rows: List[Dict[str, str]] = []
    for sid in sample_ids:
        rows.append(
            {
                "cohort_id": cohort_id,
                "sample_id": normalize_identifier(sid),
                "patient_id": "unknown",
                "timepoint_raw": "",
                "treatment_regimen_raw": "",
                "response_raw": "",
                "response_standard": "",
                "tissue_source": "",
                "source_path": str(counts_path),
                "ambiguity_note": "txt_counts_no_external_metadata",
            }
        )
    return pd.DataFrame(rows)


def read_h5ad_gz_sample_metadata(cohort_id: str, h5ad_gz_path: Path) -> pd.DataFrame:
    """Extract sample-level metadata from a gzipped h5ad.

    Decompresses to a temporary file, reads obs, then cleans up.
    """
    import tempfile
    import shutil

    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    with gzip.open(h5ad_gz_path, "rb") as fi, open(tmp_path, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    try:
        return extract_h5ad_sample_metadata(cohort_id, tmp_path)
    finally:
        tmp_path.unlink()


def read_bulk_matrix_sample_metadata(cohort_id: str, counts_path: Path) -> pd.DataFrame:
    """Extract sample-level metadata from a bulk counts matrix header."""
    with gzip.open(counts_path, "rt") as fh:
        header = fh.readline().strip().split("\t")
    # skip gene id / name columns
    sample_ids = header[2:] if len(header) > 2 else header[1:]
    rows: List[Dict[str, str]] = []
    for sid in sample_ids:
        rows.append(
            {
                "cohort_id": cohort_id,
                "sample_id": normalize_identifier(sid),
                "patient_id": "unknown",
                "timepoint_raw": "",
                "treatment_regimen_raw": "",
                "response_raw": "",
                "response_standard": "",
                "tissue_source": "",
                "source_path": str(counts_path),
                "ambiguity_note": "bulk_matrix_no_external_metadata",
            }
        )
    return pd.DataFrame(rows)


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "na", "n/a"}:
        return ""
    return text


def normalize_identifier(value: object) -> str:
    text = clean_str(value)
    if not text:
        return "unknown"
    if re.fullmatch(r"\d+\.0", text):
        return text.split(".")[0]
    return text


def decode_scalar(x: object) -> str:
    if isinstance(x, bytes):
        return x.decode("utf-8", "ignore").strip()
    return clean_str(x)


def unique_join(series: Iterable[object], sep: str = ";") -> str:
    vals = sorted({clean_str(v) for v in series if clean_str(v)})
    return sep.join(vals) if vals else "unknown"


def bool_to_str(v: bool) -> str:
    return "True" if bool(v) else "False"


def is_true(v: object) -> bool:
    return str(v).strip().lower() in {"true", "1", "t", "yes", "y"}


def classify_cancer_group(disease_text: str) -> str:
    text = disease_text.lower()
    if "hcc" in text or "liver" in text:
        return "hcc_liver"
    if text in {"", "unknown"}:
        return "unknown"
    if "/" in disease_text or ";" in disease_text:
        return "mixed"
    return "non_hcc"


def has_spatial_modality(modality_text: str) -> bool:
    t = modality_text.lower()
    return "spatial" in t or re.search(r"\bst\b", t) is not None


def has_bulk_modality(modality_text: str) -> bool:
    return "bulk" in modality_text.lower()


def has_scrna_modality(modality_text: str) -> bool:
    return "scrna" in modality_text.lower()


def has_perturb(modality_text: str, group_text: str, cohort_id: str) -> bool:
    t = f"{modality_text};{group_text};{cohort_id}".lower()
    return any(k in t for k in ["perturb", "l1000", "scperturb", "xatlas", "gse193736", "gse90063", "gse133344", "gse306429"])


def split_treatment_options(t: str) -> List[str]:
    raw = clean_str(t)
    if not raw:
        return []
    parts = [p.strip() for p in raw.split("/") if p.strip()]
    return parts


def classify_treatment_context(
    data_modality: str,
    treatment_raw: str,
    prior_treatment_raw: str,
    timepoint_schema: str,
    group_text: str,
) -> Tuple[str, str, str]:
    modality = data_modality.lower()
    treatment_options = split_treatment_options(treatment_raw)
    treatment_norm = ";".join(treatment_options).lower()
    prior = prior_treatment_raw.lower()
    timepoint = timepoint_schema.lower()

    if has_perturb(modality, group_text, ""):
        return "external_anchor_only", "perturb prior dataset", ""

    only_bulk_like = ("bulk" in modality) and ("scrna" not in modality)
    only_spatial_like = has_spatial_modality(modality) and ("scrna" not in modality)
    if only_bulk_like:
        return "external_anchor_only", "bulk/clinical anchor without scRNA", ""
    if only_spatial_like:
        return "external_anchor_only", "spatial/tissue anchor without scRNA", ""

    has_pd1_axis = any(
        k in treatment_norm
        for k in [
            "pd1",
            "pd-1",
            "pdl1",
            "pd-l1",
            "anti-pd-1",
            "anti-pd-l1",
            "atezolizumab",
            "nivolumab",
            "pembrolizumab",
        ]
    )
    has_combo = any("+" in x for x in treatment_options)

    has_locoregional_current = any(
        k in treatment_norm for k in ["tace", "haic", "rfa", "radio", "radiotherapy", "locoregional"]
    )
    overmixed_regimen = len(treatment_options) >= 4
    post_only_schema = ("post" in timepoint) and ("pre" not in timepoint) and ("on" not in timepoint)

    ambiguity_notes: List[str] = []
    if overmixed_regimen:
        ambiguity_notes.append("overmixed_regimen_options")
    if ";none" in treatment_norm or treatment_norm.endswith("none") or "none;" in treatment_norm:
        ambiguity_notes.append("contains_none_treatment_option")
    if any(k in prior for k in ["tace", "haic", "rfa", "radio", "radiotherapy"]):
        ambiguity_notes.append("prior_locoregional_therapy")

    if has_locoregional_current or post_only_schema:
        reason = "current regimen/timepoint indicates high confounding context"
        return "high_confounding_support", reason, ";".join(sorted(set(ambiguity_notes)))

    if has_pd1_axis and has_combo:
        reason = "explicit PD-1-centered combination regimen"
        return "PD1X_extension", reason, ";".join(sorted(set(ambiguity_notes)))

    if has_pd1_axis and not has_combo:
        reason = "PD-1/PD-L1 or clean ICI monotherapy axis"
        return "PD1_ICI_anchor", reason, ";".join(sorted(set(ambiguity_notes)))

    if has_scrna_modality(modality):
        reason = "scRNA cohort with non-PD1 or unclear regimen"
        return "high_confounding_support", reason, ";".join(sorted(set(ambiguity_notes)))

    return "external_anchor_only", "supportive dataset outside scRNA main training", ";".join(sorted(set(ambiguity_notes)))


def infer_source_type(resource_text: str, cohort_id: str) -> str:
    t = resource_text.lower()
    cid = cohort_id.lower()
    if "inhouse" in t or cid in {"task01", "task02", "imbrave150"}:
        return "inhouse"
    if "ncbi.nlm.nih.gov/geo" in t:
        return "public_geo"
    if any(k in t for k in ["mendeley", "singlecell.broadinstitute", "trace.ncbi.nlm.nih.gov", "lambrechtslab", "cell.com", "nature.com", "pubmed"]):
        return "public_repository"
    if t in {"", "unknown"}:
        return "unknown"
    return "public_other"


def infer_response_quality(outcome_text: str, has_os_pfs_text: str) -> str:
    outcome = outcome_text.lower()
    survival = has_os_pfs_text.lower()
    if any(k in outcome for k in ["mrecist", "recist", "patho", "pcr", "mri", "clinical", "expansion", "r/nr", "cb", "ncb"]):
        return "high"
    if outcome not in {"", "unknown", "none"}:
        return "medium"
    if survival not in {"", "unknown", "none"}:
        return "medium"
    return "low"


def map_response_binary(response_raw: str) -> str:
    r = clean_str(response_raw).strip().upper()
    if r in {"", "UNKNOWN", "NA", "N/A", "UN", "NE", "NOT AVAILABLE"}:
        return "unknown"
    if r in {"CR", "PR", "R", "RESPONDER", "MPR", "PCR", "CB", "ICB", "BENEFIT"}:
        return "responder"
    if r in {"PD", "NR", "NON-RESPONDER", "NON_RESPONDER", "NCB", "NONBENEFIT", "NON-BENEFIT", "SD"}:
        return "non_responder"
    if any(k in r for k in ["RESPONDER", "BENEFIT"]) and "NON" not in r:
        return "responder"
    if any(k in r for k in ["NON", "PROGRESS", "RESIST"]) or r == "NR":
        return "non_responder"
    return "unknown"


def response_order(response_raw: str) -> str:
    r = clean_str(response_raw).upper()
    mapping = {"CR": "4", "PR": "3", "R": "3", "SD": "2", "NON CR/PD": "2", "PD": "1", "NR": "1"}
    return mapping.get(r, "unknown")


def normalize_timepoint(timepoint_raw: str, timepoint_std: str = "") -> str:
    raw = clean_str(timepoint_std).lower().replace("_", "-") or clean_str(timepoint_raw).lower().replace("_", "-")
    if raw in {"pre", "pre-tx", "pre-treatment", "pretreatment", "t0", "baseline"}:
        return "pre"
    if raw in {"on-treatment", "on treatment", "t1", "b1", "b2", "b3"}:
        return "on_treatment"
    if raw in {"post", "post-tx", "post-treatment", "post treatment", "t2"}:
        return "post"
    if "long" in raw and "post" in raw:
        return "long_post"
    if "pre" in raw:
        return "pre"
    if "post" in raw:
        return "post"
    if "follow" in raw and "up" in raw:
        return "post"
    if any(k in raw for k in ["on", "during", "week", "cycle", "b1", "b2", "b3"]):
        return "on_treatment"
    return "unknown"


def cohort_to_h5ad_path(cohort_id: str) -> Optional[Path]:
    cid = cohort_id.lower()
    candidates = [
        ROOT / "data" / "processed" / "srt" / "raw" / f"{cid}.h5ad",
        ROOT / "data" / "processed" / "srt" / "raw" / f"{cid.replace('_', '')}.h5ad",
        ROOT / "mvp" / "data" / f"{cohort_id}.h5ad",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


# ---------------------------------------------------------------------------
# Subset expansion for multi-subset cohorts
# ---------------------------------------------------------------------------
SUBSET_SPEC = {
    "GSE123813": [
        {"subset_id": "GSE123813_bcc", "match_fn": lambda row: row.get("disease", "").lower().strip() == "bcc"},
        {"subset_id": "GSE123813_scc", "match_fn": lambda row: row.get("disease", "").lower().strip() == "scc"},
    ],
    "GSE140228": [
        {"subset_id": "GSE140228_droplet", "match_fn": lambda row: True},
        {"subset_id": "GSE140228_smartseq2", "match_fn": lambda row: True},
    ],
    "GSE207422": [
        {"subset_id": "GSE207422_sc", "match_fn": lambda row: "scrna" in row.get("data_type", "").lower()},
        {"subset_id": "GSE207422_bulk", "match_fn": lambda row: "bulk" in row.get("data_type", "").lower() and "scrna" not in row.get("data_type", "").lower()},
    ],
    "GSE235863": [
        {"subset_id": "GSE235863_9pt_cd45", "match_fn": lambda row: "cd45" in row.get("celltype", "").lower()},
        {"subset_id": "GSE235863_5pt_cd8", "match_fn": lambda row: "cd8" in row.get("celltype", "").lower()},
        {"subset_id": "GSE235863_bulk", "match_fn": lambda row: "bulk" in row.get("data_type", "").lower() and "scrna" not in row.get("data_type", "").lower()},
    ],
    "lambrecht_brca": [
        {"subset_id": "lambrecht_brca_cohort1_cells", "match_fn": lambda row: True},
        {"subset_id": "lambrecht_brca_cohort2_cells", "match_fn": lambda row: True},
        {"subset_id": "lambrecht_brca_cohort1_tcell", "match_fn": lambda row: True},
        {"subset_id": "lambrecht_brca_cohort1_myeloid", "match_fn": lambda row: True},
        {"subset_id": "lambrecht_brca_cohort1_dc", "match_fn": lambda row: True},
    ],
}


def expand_data_collection_subsets(dc: pd.DataFrame) -> pd.DataFrame:
    """Duplicate rows for cohorts that need per-subset registry entries."""
    rows = []
    for _, row in dc.iterrows():
        cid = clean_str(row.get("id", ""))
        specs = SUBSET_SPEC.get(cid)
        if not specs:
            rows.append(row.to_dict())
            continue
        matched = False
        for spec in specs:
            if spec["match_fn"](row.to_dict()):
                d = row.to_dict()
                d["id"] = spec["subset_id"]
                rows.append(d)
                matched = True
        if not matched:
            rows.append(row.to_dict())
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Bulk metadata extraction helpers
# ---------------------------------------------------------------------------
def extract_bulk_sample_metadata(cohort_id: str) -> pd.DataFrame:
    """Return sample-level metadata DataFrame for bulk-only cohorts."""
    rows: List[Dict[str, str]] = []
    if cohort_id == "GSE207422_bulk":
        meta_path = ROOT / "data" / "combo" / "GSE207422" / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"
        if meta_path.exists():
            df = pd.read_excel(meta_path, sheet_name="sheet1", dtype=str).fillna("")
            for _, r in df.iterrows():
                rows.append(
                    {
                        "cohort_id": cohort_id,
                        "sample_id": normalize_identifier(r.get("Sample", "")),
                        "patient_id": normalize_identifier(r.get("Patient", "")),
                        "timepoint_raw": clean_str(r.get("Resource", "")),
                        "treatment_regimen_raw": clean_str(r.get("PD1 Antibody", "")),
                        "response_raw": clean_str(r.get("RECIST", "")),
                        "response_standard": "",
                        "tissue_source": "tumor",
                        "source_path": str(meta_path),
                        "ambiguity_note": "",
                    }
                )
    elif cohort_id == "GSE235863_bulk":
        counts_path = ROOT / "data" / "combo" / "GSE235863" / "GSE235863_bulk_rna_seq_raw_counts.txt.gz"
        if counts_path.exists():
            with gzip.open(counts_path, "rt") as fh:
                header = fh.readline().strip().split("\t")
            for sid in header[2:]:
                rows.append(
                    {
                        "cohort_id": cohort_id,
                        "sample_id": normalize_identifier(sid),
                        "patient_id": "unknown",
                        "timepoint_raw": "post",
                        "treatment_regimen_raw": "unknown",
                        "response_raw": "unknown",
                        "response_standard": "",
                        "tissue_source": "tumor",
                        "source_path": str(counts_path),
                        "ambiguity_note": "bulk_metadata_not_available_locally",
                    }
                )
    return pd.DataFrame(rows)


class ObsColumnReader:
    def __init__(self, obs_group: h5py.Group, col: str):
        self.col = col
        self.obj = obs_group[col]
        self.is_categorical = isinstance(self.obj, h5py.Group) and {"categories", "codes"}.issubset(set(self.obj.keys()))
        if self.is_categorical:
            self.categories = [decode_scalar(x) for x in self.obj["categories"][:]]

    def read_slice(self, start: int, end: int) -> List[str]:
        if self.is_categorical:
            codes = self.obj["codes"][start:end]
            out: List[str] = []
            for c in codes:
                ci = int(c)
                if ci < 0 or ci >= len(self.categories):
                    out.append("")
                else:
                    out.append(clean_str(self.categories[ci]))
            return out

        values = self.obj[start:end]
        return [decode_scalar(x) for x in values]


def pick_column(candidates: List[str], available: List[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in available}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def _obs_nrows(obs: h5py.Group) -> int:
    """Return number of rows in obs, robust to missing _index dataset."""
    idx_name = obs.attrs.get("_index", "_index")
    if idx_name in obs:
        return len(obs[idx_name])
    # fallback: use first available column
    for k in obs.keys():
        return len(obs[k])
    return 0


def _preview_values(obs: h5py.Group, col: str, n: int = 5000) -> List[str]:
    reader = ObsColumnReader(obs, col)
    size = _obs_nrows(obs)
    end = min(size, n)
    return [clean_str(v) for v in reader.read_slice(0, end)]


def _looks_like_cell_barcode(text: str) -> bool:
    t = text.strip()
    if not t:
        return False
    if re.search(r"[ACGT]{12,}", t):
        return True
    if any(k in t for k in ["AAAC", "TTTG", "GGTT", "barcode"]):
        return True
    if len(t) >= 35 and "_" in t:
        return True
    return False


def choose_sample_column(obs: h5py.Group, available_cols: List[str]) -> Optional[str]:
    # Prefer true sample-level columns; avoid cell-barcode-like sample_id fields.
    candidates = ["sample_id", "sample", "sample_id_raw", "sample_prefix", "donor_id", "patient_id", "patient", "Patient_ID", "orig.ident", "lib"]
    for cand in candidates:
        col = pick_column([cand], available_cols)
        if col is None:
            continue
        preview = [v for v in _preview_values(obs, col, n=5000) if v]
        if not preview:
            continue
        uniq_ratio = len(set(preview)) / max(1, len(preview))
        barcode_ratio = sum(1 for v in preview[:500] if _looks_like_cell_barcode(v)) / max(1, min(500, len(preview)))
        # If the would-be sample column behaves like per-cell identifier, skip it.
        if cand in {"sample_id", "sample", "sample_id_raw"} and uniq_ratio > 0.95 and barcode_ratio > 0.2:
            continue
        return col
    return None


def extract_h5ad_sample_metadata(cohort_id: str, h5ad_path: Path) -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            return pd.DataFrame()
        obs = f["obs"]
        cols = list(obs.keys())

        sample_col = choose_sample_column(obs, cols)
        patient_col = pick_column(["patient_id", "patient", "Patient_ID", "donor_id", "donor", "subject_id"], cols)
        time_col = pick_column(["timepoint", "Timepoint", "tp", "timepoint_raw"], cols)
        treatment_col = pick_column(["treatment", "Treatment", "treatment_raw", "sample_treatment"], cols)
        response_col = pick_column(["response", "Response", "sample_treatment_Resp", "path_response_raw", "Path_response", "response_raw"], cols)
        response_std_col = pick_column(["resp_standard", "response_standard", "response_criteria", "criteria"], cols)
        tissue_col = pick_column(["sample_tissue", "tissue", "sample_disease"], cols)

        if sample_col is None:
            return pd.DataFrame()

        readers: Dict[str, ObsColumnReader] = {"sample_id": ObsColumnReader(obs, sample_col)}
        if patient_col:
            readers["patient_id"] = ObsColumnReader(obs, patient_col)
        if time_col:
            readers["timepoint_raw"] = ObsColumnReader(obs, time_col)
        if treatment_col:
            readers["treatment_regimen_raw"] = ObsColumnReader(obs, treatment_col)
        if response_col:
            readers["response_raw"] = ObsColumnReader(obs, response_col)
        if response_std_col:
            readers["response_standard"] = ObsColumnReader(obs, response_std_col)
        if tissue_col:
            readers["tissue_source"] = ObsColumnReader(obs, tissue_col)

        n = _obs_nrows(obs)
        chunk = 200000
        agg: Dict[str, Dict[str, set]] = {}

        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            values = {k: r.read_slice(start, end) for k, r in readers.items()}
            for i in range(end - start):
                sid = normalize_identifier(values["sample_id"][i])
                if sid == "unknown":
                    continue
                rec = agg.setdefault(
                    sid,
                    {
                        "patient_id": set(),
                        "timepoint_raw": set(),
                        "treatment_regimen_raw": set(),
                        "response_raw": set(),
                        "response_standard": set(),
                        "tissue_source": set(),
                    },
                )
                for field in rec.keys():
                    if field in values:
                        v = clean_str(values[field][i])
                        if v and v.upper() not in {"NA", "NAN", "NONE"}:
                            rec[field].add(v)

        out_rows: List[Dict[str, str]] = []
        for sid, rec in agg.items():
            pid_vals = sorted(rec["patient_id"]) if rec["patient_id"] else []
            pid = pid_vals[0] if pid_vals else sid.split(".")[0]

            ambiguity = []
            for key in ["patient_id", "timepoint_raw", "treatment_regimen_raw", "response_raw", "response_standard"]:
                if len(rec[key]) > 1:
                    ambiguity.append(f"multi_{key}")

            out_rows.append(
                {
                    "cohort_id": cohort_id,
                    "sample_id": sid,
                    "patient_id": normalize_identifier(pid),
                    "timepoint_raw": ";".join(sorted(rec["timepoint_raw"])),
                    "treatment_regimen_raw": ";".join(sorted(rec["treatment_regimen_raw"])),
                    "response_raw": ";".join(sorted(rec["response_raw"])),
                    "response_standard": ";".join(sorted(rec["response_standard"])),
                    "tissue_source": ";".join(sorted(rec["tissue_source"])),
                    "source_path": str(h5ad_path),
                    "ambiguity_note": ";".join(ambiguity),
                }
            )

        return pd.DataFrame(out_rows)


def first_non_empty(vals: List[str]) -> str:
    for v in vals:
        if clean_str(v):
            return clean_str(v)
    return ""


def merge_values(vals: List[str]) -> Tuple[str, str]:
    uniq = sorted({clean_str(v) for v in vals if clean_str(v)})
    if not uniq:
        return "", ""
    if len(uniq) == 1:
        return uniq[0], ""
    return uniq[0], ";".join(uniq)


def write_manifest(outputs: List[str], known_missing: List[str]) -> None:
    lines = [
        'manifest_version: "v6.1"',
        "execution_records:",
        '  - step_id: "v6_1_step1_data_governance"',
        '    step_alias: "cohort_registry_and_metadata_master_build"',
        "    inputs:",
        '      - "docs/data_collection.csv"',
        '      - "mvp/outputs/patient_metadata_master.csv"',
        '      - "mvp/sidecars/*/sample_summary.csv"',
        '      - "data/processed/srt/raw/*.h5ad"',
        '      - "data/** (local existing assets only)"',
        "    outputs:",
    ]
    for p in outputs:
        lines.append(f'      - "{p}"')

    lines += [
        "    major_fields:",
        '      - "patient_id / sample_id / timepoint_normalized / treatment_context / response_binary / response_standard"',
        '      - "cohort-level role assignment and exclusion reasons"',
        '      - "data leakage risk checks"',
        "    known_missing:",
    ]
    for item in known_missing:
        lines.append(f'      - "{item}"')

    lines += [
        "    downstream_dependencies:",
        '      - "Step2 immune-state measurement (without modeling)"',
        '      - "Step3 strong baseline package"',
        '      - "Step4 module discovery and PD-1 anchor logic"',
        "",
    ]

    MANIFEST_OUT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    dc = pd.read_csv(DATA_COLLECTION_PATH, dtype=str).fillna("")
    if "Unnamed: 26" in dc.columns:
        dc = dc.drop(columns=["Unnamed: 26"])
    dc.columns = [c.strip() for c in dc.columns]
    dc["id"] = dc["id"].map(clean_str)
    dc = dc[dc["id"] != ""].copy()

    pm = pd.read_csv(MVP_PATIENT_METADATA_PATH, dtype=str).fillna("") if MVP_PATIENT_METADATA_PATH.exists() else pd.DataFrame()

    # Cohort registry
    dc = expand_data_collection_subsets(dc)

    cohort_rows: List[Dict[str, object]] = []
    local_source_records: List[Dict[str, str]] = []

    for cohort_id, g in dc.groupby("id", sort=True):
        modalities = unique_join(g["data_type"])
        diseases = unique_join(g["disease"])
        sample_type = unique_join(g["tissue"])
        platform = unique_join(g["protocal"])
        resource = unique_join(g["resource"])
        treatment_raw = unique_join(g["treatment"])
        prior_raw = unique_join(g["prior_treatment"])
        timepoint_schema = unique_join(g["timepoint"])
        response_label_type = unique_join(g["has_outcome"])
        has_os_pfs = unique_join(g["has_os_pfs"])
        group_text = unique_join(g["group"])

        context, context_reason, ambiguity = classify_treatment_context(modalities, treatment_raw, prior_raw, timepoint_schema, group_text)

        raw_paths: List[str] = []
        assets = cohort_to_input_assets(cohort_id)
        h5ad: Optional[Path] = None
        for asset in assets:
            raw_paths.append(str(asset.path))
            local_source_records.append(
                {
                    "cohort_id": cohort_id,
                    "source_type": f"local_{asset.asset_type}",
                    "source_path": str(asset.path),
                    "purpose": "sample-level metadata extraction (patient/sample/timepoint/treatment/response/response_standard)",
                    "status": "used",
                }
            )
            if asset.asset_type == "h5ad":
                h5ad = asset.path

        mvp_data_h5ad = ROOT / "mvp" / "data" / f"{cohort_id}.h5ad"
        if mvp_data_h5ad.exists() and str(mvp_data_h5ad) not in raw_paths:
            raw_paths.append(str(mvp_data_h5ad))

        for p in [ROOT / "data" / "combo" / cohort_id, ROOT / "data" / "imm" / cohort_id, ROOT / "data" / "perturb_seq" / cohort_id]:
            if p.exists():
                raw_paths.append(str(p))

        # Fallback for subset cohorts (e.g. GSE123813_bcc → data/imm/GSE123813)
        if "_" in cohort_id:
            parent = cohort_id.split("_")[0]
            for p in [ROOT / "data" / "combo" / parent, ROOT / "data" / "imm" / parent, ROOT / "data" / "perturb_seq" / parent]:
                if p.exists() and str(p) not in raw_paths:
                    raw_paths.append(str(p))

        processed_candidates = sorted((ROOT / "mvp" / "outputs" / "_phase3_tmp").glob(f"{cohort_id}_*"))
        processed_paths = [str(p) for p in processed_candidates]
        if h5ad:
            processed_paths.append(str(h5ad))

        sidecar_path = ROOT / "mvp" / "sidecars" / cohort_id
        clinical_candidates = []
        if not pm.empty and (pm["cohort_id"].map(clean_str) == cohort_id).any():
            clinical_candidates.append(str(MVP_PATIENT_METADATA_PATH) + f" (filter: cohort_id={cohort_id})")
        if (ROOT / "data" / "raw" / "clinical" / cohort_id).exists():
            clinical_candidates.append(str(ROOT / "data" / "raw" / "clinical" / cohort_id))

        source_type = infer_source_type(resource, cohort_id)
        response_quality = infer_response_quality(response_label_type, has_os_pfs)

        spatial_available = has_spatial_modality(modalities)
        bulk_available = has_bulk_modality(modalities)
        perturb_available = has_perturb(modalities, group_text, cohort_id)
        scrna_available = has_scrna_modality(modalities)

        usable_pd1_anchor = scrna_available and context == "PD1_ICI_anchor"
        usable_pd1x = scrna_available and context == "PD1X_extension"
        usable_hcc = ("hcc" in diseases.lower()) and scrna_available
        usable_shared = scrna_available and context in {"PD1_ICI_anchor", "PD1X_extension"}
        usable_external = context == "external_anchor_only" or bulk_available or spatial_available or perturb_available
        usable_spatial = spatial_available

        exclusion = []
        if context in {"high_confounding_support", "external_anchor_only"}:
            exclusion.append(f"excluded_from_main_training_{context}")
        if not raw_paths and not sidecar_path.exists():
            exclusion.append("not_locally_ingested")

        notes = [context_reason]
        if ambiguity:
            notes.append(f"ambiguity:{ambiguity}")
        if clean_str(prior_raw) and prior_raw.lower() not in {"none", "unknown"}:
            notes.append(f"prior_treatment={prior_raw}")

        cohort_rows.append(
            {
                "cohort_id": cohort_id,
                "cohort_name": cohort_id,
                "disease": diseases,
                "cancer_group": classify_cancer_group(diseases),
                "source_type": source_type,
                "data_modality": modalities,
                "sample_type": sample_type,
                "platform": platform,
                "center_or_study": cohort_id,
                "raw_data_path": ";".join(sorted(set(raw_paths))) if raw_paths else "unknown",
                "processed_data_path": ";".join(sorted(set(processed_paths))) if processed_paths else "unknown",
                "sidecar_path": str(sidecar_path) if sidecar_path.exists() else "unknown",
                "clinical_metadata_path": ";".join(clinical_candidates) if clinical_candidates else "unknown",
                "treatment_regimen_raw": treatment_raw,
                "treatment_context": context,
                "timepoint_schema": timepoint_schema,
                "response_label_type": response_label_type,
                "response_label_quality": response_quality,
                "paired_available": bool_to_str("paired" in unique_join(g["tp_paired"]).lower()),
                "spatial_available": bool_to_str(spatial_available),
                "bulk_available": bool_to_str(bulk_available),
                "perturb_available": bool_to_str(perturb_available),
                "usable_for_PD1_anchor": bool_to_str(usable_pd1_anchor),
                "usable_for_PD1X_extension": bool_to_str(usable_pd1x),
                "usable_for_HCC_specific": bool_to_str(usable_hcc),
                "usable_for_shared_module": bool_to_str(usable_shared),
                "usable_for_external_anchor": bool_to_str(usable_external),
                "usable_for_spatial_adjudication": bool_to_str(usable_spatial),
                "exclusion_reason": ";".join(sorted(set(exclusion))),
                "notes": "; ".join(sorted(set(notes))),
            }
        )

    cohort_registry = pd.DataFrame(cohort_rows).sort_values("cohort_id").reset_index(drop=True)
    cohort_registry_path = OUT_DIR / "cohort_registry_v6_1.csv"
    cohort_registry.to_csv(cohort_registry_path, index=False)

    # Extract local sample metadata from discovered assets
    asset_frames: List[pd.DataFrame] = []
    for cid in cohort_registry["cohort_id"].tolist():
        assets = cohort_to_input_assets(cid)
        for asset in assets:
            try:
                ext = pd.DataFrame()
                if asset.asset_type == "h5ad":
                    ext = extract_h5ad_sample_metadata(cid, asset.path)
                elif asset.asset_type == "mtx10x":
                    ext = read_mtx10x_sample_metadata(cid, asset.path)
                elif asset.asset_type == "txt_counts":
                    ext = read_txt_counts_sample_metadata(cid, asset.path)
                elif asset.asset_type == "h5ad_gz":
                    ext = read_h5ad_gz_sample_metadata(cid, asset.path)
                elif asset.asset_type == "bulk_matrix":
                    ext = read_bulk_matrix_sample_metadata(cid, asset.path)
                if not ext.empty:
                    asset_frames.append(ext)
            except Exception as e:
                local_source_records.append(
                    {
                        "cohort_id": cid,
                        "source_type": f"local_{asset.asset_type}",
                        "source_path": str(asset.path),
                        "purpose": "sample-level metadata extraction",
                        "status": f"error:{e}",
                    }
                )

    h5ad_sample = pd.concat(asset_frames, ignore_index=True) if asset_frames else pd.DataFrame()

    # Bulk-only metadata extraction (no h5ad)
    bulk_frames = []
    for cid in cohort_registry["cohort_id"].tolist():
        if "bulk" in cid.lower():
            try:
                ext = extract_bulk_sample_metadata(cid)
                if not ext.empty:
                    bulk_frames.append(ext)
                    local_source_records.append(
                        {
                            "cohort_id": cid,
                            "source_type": "local_bulk_metadata",
                            "source_path": str(ext["source_path"].iloc[0]),
                            "purpose": "sample-level metadata extraction for bulk cohort",
                            "status": "used",
                        }
                    )
            except Exception as e:
                local_source_records.append(
                    {
                        "cohort_id": cid,
                        "source_type": "local_bulk_metadata",
                        "source_path": "unknown",
                        "purpose": "sample-level metadata extraction for bulk cohort",
                        "status": f"error:{e}",
                    }
                )
    if bulk_frames:
        h5ad_sample = pd.concat([h5ad_sample, *bulk_frames], ignore_index=True)

    # Source A: mvp patient metadata
    source_a_rows: List[Dict[str, str]] = []
    if not pm.empty:
        for _, r in pm.iterrows():
            source_a_rows.append(
                {
                    "cohort_id": clean_str(r.get("cohort_id", "")),
                    "sample_id": normalize_identifier(r.get("sample_id", "")),
                    "patient_id": normalize_identifier(r.get("patient_id", "")),
                    "original_sample_id": clean_str(r.get("sample_id", "")),
                    "timepoint_raw": clean_str(r.get("timepoint_raw", "")),
                    "timepoint_standardized": clean_str(r.get("timepoint_standardized", "")),
                    "treatment_regimen_raw": clean_str(r.get("therapy_class", "")),
                    "response_raw": clean_str(r.get("response_raw", "")),
                    "response_binary": clean_str(r.get("response_binary", "")),
                    "response_standard": "",
                    "tissue_source": clean_str(r.get("tissue_source", "")),
                    "source_priority": "2",
                    "source_name": "mvp_patient_metadata_master",
                    "ambiguity_note": "",
                }
            )

    source_a = pd.DataFrame(source_a_rows)

    # Source B: sidecar summaries
    source_b_rows: List[Dict[str, str]] = []
    for summary_path in sorted(MVP_SIDECARS_DIR.glob("*/sample_summary.csv")):
        sid = summary_path.parent.name
        sdf = pd.read_csv(summary_path, dtype=str).fillna("")
        for _, row in sdf.iterrows():
            sample_id = normalize_identifier(row.get("sample_id", ""))
            patient_id = normalize_identifier(row.get("patient_id", ""))
            if patient_id == "unknown" and sample_id != "unknown":
                patient_id = normalize_identifier(sample_id.split(".")[0])
            source_b_rows.append(
                {
                    "cohort_id": sid,
                    "sample_id": sample_id,
                    "patient_id": patient_id,
                    "original_sample_id": clean_str(row.get("sample_id", "")) or sample_id,
                    "timepoint_raw": clean_str(row.get("timepoint", "")),
                    "timepoint_standardized": "",
                    "treatment_regimen_raw": "",
                    "response_raw": "",
                    "response_binary": "",
                    "response_standard": "",
                    "tissue_source": "",
                    "source_priority": "3",
                    "source_name": "mvp_sidecar_sample_summary",
                    "ambiguity_note": "",
                }
            )

    source_b = pd.DataFrame(source_b_rows)

    # Source C: h5ad extracted sample metadata
    source_c = pd.DataFrame()
    if not h5ad_sample.empty:
        source_c = h5ad_sample.copy()
        source_c["original_sample_id"] = source_c["sample_id"]
        source_c["timepoint_standardized"] = ""
        source_c["response_binary"] = ""
        source_c["source_priority"] = "1"
        source_c["source_name"] = "local_h5ad_obs"

    all_sources = pd.concat([df for df in [source_c, source_a, source_b] if not df.empty], ignore_index=True)
    all_sources = all_sources[all_sources["cohort_id"].map(clean_str) != ""].copy()
    all_sources["sample_id"] = all_sources["sample_id"].map(normalize_identifier)
    all_sources["patient_id"] = all_sources["patient_id"].map(normalize_identifier)

    grouped_rows: List[Dict[str, str]] = []
    for (cohort_id, sample_id), g in all_sources.groupby(["cohort_id", "sample_id"], dropna=False):
        if sample_id == "unknown":
            # retain unknown sample ids only when no identifiable rows exist
            if len(g) > 1:
                continue

        g = g.sort_values("source_priority")
        patient_v, patient_amb = merge_values(g["patient_id"].tolist())
        if not patient_v:
            patient_v = normalize_identifier(sample_id.split(".")[0]) if sample_id != "unknown" else "unknown"

        time_raw_v, time_raw_amb = merge_values(g["timepoint_raw"].tolist())
        time_std_v, _ = merge_values(g["timepoint_standardized"].tolist())
        treat_v, treat_amb = merge_values(g["treatment_regimen_raw"].tolist())
        resp_raw_v, resp_raw_amb = merge_values(g["response_raw"].tolist())
        resp_std_v, resp_std_amb = merge_values(g["response_standard"].tolist())
        tissue_v, tissue_amb = merge_values(g["tissue_source"].tolist())

        original_sample = first_non_empty(g["original_sample_id"].tolist()) or sample_id

        ambiguity = [x for x in [patient_amb, time_raw_amb, treat_amb, resp_raw_amb, resp_std_amb, tissue_amb] if x]
        source_names = sorted(set(g["source_name"].tolist()))

        grouped_rows.append(
            {
                "cohort_id": cohort_id,
                "sample_id": sample_id,
                "patient_id": patient_v,
                "original_sample_id": original_sample,
                "timepoint_raw": time_raw_v,
                "timepoint_standardized": time_std_v,
                "treatment_regimen_raw": treat_v,
                "response_raw": resp_raw_v,
                "response_standard": resp_std_v,
                "response_binary": first_non_empty(g["response_binary"].tolist()),
                "tissue_source": tissue_v,
                "source_name": ";".join(source_names),
                "source_ambiguity": ";".join(ambiguity),
            }
        )

    merged_sample_base = pd.DataFrame(grouped_rows)

    # Final sample metadata master
    cohort_map = cohort_registry.set_index("cohort_id").to_dict(orient="index")

    final_rows: List[Dict[str, object]] = []
    for _, r in merged_sample_base.iterrows():
        cohort_id = clean_str(r["cohort_id"])
        if cohort_id not in cohort_map:
            continue
        c = cohort_map[cohort_id]

        disease = clean_str(c.get("disease", "unknown")) or "unknown"
        sample_type = clean_str(c.get("sample_type", "unknown")) or "unknown"
        treatment_context = clean_str(c.get("treatment_context", "unknown")) or "unknown"

        treatment_raw = clean_str(r.get("treatment_regimen_raw", "")) or clean_str(c.get("treatment_regimen_raw", "")) or "unknown"
        timepoint_raw = clean_str(r.get("timepoint_raw", "")) or "unknown"
        timepoint_normalized = normalize_timepoint(timepoint_raw, clean_str(r.get("timepoint_standardized", "")))

        response_raw = clean_str(r.get("response_raw", "")) or "unknown"
        response_binary = clean_str(r.get("response_binary", "")).replace("-", "_")
        if response_binary not in {"responder", "non_responder", "unknown", "not_applicable"}:
            response_binary = map_response_binary(response_raw)
        cohort_response_label_type = clean_str(c.get("response_label_type", "")).lower()
        if response_binary == "unknown" and cohort_response_label_type in {"none", "", "unknown"}:
            response_binary = "not_applicable"

        response_standard = clean_str(r.get("response_standard", "")) or clean_str(c.get("response_label_type", "")) or "unknown"

        sc_available = has_scrna_modality(clean_str(c.get("data_modality", "")))
        bulk_available = is_true(c.get("bulk_available", "False"))
        spatial_available = is_true(c.get("spatial_available", "False"))

        usable_main = (
            sc_available
            and treatment_context in {"PD1_ICI_anchor", "PD1X_extension"}
            and normalize_identifier(r.get("sample_id", "")) != "unknown"
            and timepoint_normalized != "unknown"
        )
        usable_sens = sc_available and (treatment_context == "high_confounding_support")

        exclusion = []
        if normalize_identifier(r.get("sample_id", "")) == "unknown":
            exclusion.append("missing_sample_id")
        if timepoint_normalized == "unknown":
            exclusion.append("unknown_timepoint")
        if treatment_context == "high_confounding_support":
            exclusion.append("high_confounding_context")
        if treatment_context == "external_anchor_only":
            exclusion.append("external_anchor_only")
        if not sc_available:
            exclusion.append("not_scRNA")

        if clean_str(r.get("source_ambiguity", "")):
            exclusion.append(f"source_ambiguity:{clean_str(r.get('source_ambiguity',''))}")

        final_rows.append(
            {
                "cohort_id": cohort_id,
                "patient_id": normalize_identifier(r.get("patient_id", "")),
                "sample_id": normalize_identifier(r.get("sample_id", "")),
                "original_sample_id": clean_str(r.get("original_sample_id", "")) or normalize_identifier(r.get("sample_id", "")),
                "disease": disease,
                "tissue_source": clean_str(r.get("tissue_source", "")) or sample_type,
                "sample_type": sample_type,
                "timepoint_raw": timepoint_raw,
                "timepoint_normalized": timepoint_normalized,
                "treatment_regimen_raw": treatment_raw,
                "treatment_context": treatment_context,
                "response_standard": response_standard,
                "response_raw": response_raw,
                "response_binary": response_binary,
                "response_ordered": response_order(response_raw),
                "survival_available": bool_to_str(clean_str(c.get("response_label_quality", "")) in {"high", "medium"}),
                "spatial_available": bool_to_str(spatial_available),
                "bulk_available": bool_to_str(bulk_available),
                "scRNA_available": bool_to_str(sc_available),
                "usable_in_main_analysis": bool_to_str(usable_main),
                "usable_in_sensitivity_analysis": bool_to_str(usable_sens),
                "exclusion_reason": ";".join(sorted(set(exclusion))),
            }
        )

    sample_metadata = pd.DataFrame(final_rows)
    sample_metadata = sample_metadata.drop_duplicates(subset=["cohort_id", "patient_id", "sample_id"]).sort_values(["cohort_id", "patient_id", "sample_id"]).reset_index(drop=True)

    cohort_treatment_map = cohort_registry.set_index("cohort_id")["treatment_regimen_raw"].to_dict()
    treatment_missing_mask = sample_metadata["treatment_regimen_raw"].isin(["", "unknown"])
    sample_metadata.loc[treatment_missing_mask, "treatment_regimen_raw"] = (
        sample_metadata.loc[treatment_missing_mask, "cohort_id"].map(cohort_treatment_map).fillna("unknown")
    )
    sample_metadata["treatment_regimen_raw"] = sample_metadata["treatment_regimen_raw"].replace("", "unknown")

    allowed_tp = {"pre", "on_treatment", "post", "long_post", "unknown"}
    sample_metadata.loc[~sample_metadata["timepoint_normalized"].isin(allowed_tp), "timepoint_normalized"] = "unknown"

    allowed_rb = {"responder", "non_responder", "unknown", "not_applicable"}
    sample_metadata.loc[~sample_metadata["response_binary"].isin(allowed_rb), "response_binary"] = "unknown"

    sample_metadata_path = OUT_DIR / "sample_metadata_master_v6_1.csv"
    sample_metadata.to_csv(sample_metadata_path, index=False)

    # Patient metadata
    patient_rows: List[Dict[str, object]] = []
    for (cohort_id, patient_id), g in sample_metadata.groupby(["cohort_id", "patient_id"], dropna=False):
        disease = unique_join(g["disease"].tolist())
        context_mode = g["treatment_context"].mode().iloc[0] if len(g["treatment_context"].mode()) else "unknown"

        has_pre = (g["timepoint_normalized"] == "pre").any()
        has_on = (g["timepoint_normalized"] == "on_treatment").any()
        has_post = g["timepoint_normalized"].isin(["post", "long_post"]).any()

        score = {"CR": 4, "PR": 3, "R": 3, "SD": 2, "PD": 1, "NR": 1}
        best_raw = "unknown"
        best_score = -1
        for rr in g["response_raw"].tolist():
            k = clean_str(rr).upper()
            if k in score and score[k] > best_score:
                best_score = score[k]
                best_raw = k

        if best_raw == "unknown":
            binaries = list(g["response_binary"].unique())
            if "responder" in binaries:
                best_binary = "responder"
            elif "non_responder" in binaries:
                best_binary = "non_responder"
            else:
                best_binary = "unknown"
        else:
            best_binary = "responder" if best_raw in {"CR", "PR", "R"} else "non_responder"

        std_mode = g[g["response_standard"].map(clean_str) != ""]["response_standard"].mode()
        response_standard_primary = std_mode.iloc[0] if len(std_mode) else "unknown"

        modalities = int((g["scRNA_available"] == "True").any()) + int((g["bulk_available"] == "True").any()) + int((g["spatial_available"] == "True").any())

        patient_rows.append(
            {
                "cohort_id": cohort_id,
                "patient_id": patient_id,
                "disease": disease,
                "treatment_context_primary": context_mode,
                "has_pre_sample": bool_to_str(has_pre),
                "has_on_treatment_sample": bool_to_str(has_on),
                "has_post_sample": bool_to_str(has_post),
                "has_paired_pre_post": bool_to_str(has_pre and has_post),
                "best_response_raw": best_raw,
                "best_response_binary": best_binary,
                "response_standard_primary": response_standard_primary,
                "survival_available": bool_to_str((g["survival_available"] == "True").any()),
                "number_of_samples": int(g["sample_id"].nunique()),
                "number_of_modalities": modalities,
                "use_for_PD1_anchor": bool_to_str(((g["usable_in_main_analysis"] == "True") & (g["treatment_context"] == "PD1_ICI_anchor")).any()),
                "use_for_PD1X_extension": bool_to_str(((g["usable_in_main_analysis"] == "True") & (g["treatment_context"] == "PD1X_extension")).any()),
                "use_for_HCC_specific": bool_to_str(("hcc" in disease.lower()) and ((g["scRNA_available"] == "True").any())),
                "use_for_external_anchor": bool_to_str((g["treatment_context"] == "external_anchor_only").any()),
                "notes": "",
            }
        )

    patient_metadata = pd.DataFrame(patient_rows).sort_values(["cohort_id", "patient_id"]).reset_index(drop=True)
    patient_metadata_path = OUT_DIR / "patient_metadata_master_v6_1.csv"
    patient_metadata.to_csv(patient_metadata_path, index=False)

    # Treatment context flags
    reason_map = {
        "PD1_ICI_anchor": "PD-1/PD-L1 monotherapy or clean ICI anchor axis",
        "PD1X_extension": "PD-1-centered combination context",
        "high_confounding_support": "high confounding context (regimen/timepoint complexity)",
        "external_anchor_only": "external anchor dataset outside main scRNA training",
    }

    flag_rows: List[Dict[str, str]] = []
    for _, r in sample_metadata.iterrows():
        ambiguities = []
        regimen = clean_str(r["treatment_regimen_raw"]).lower()
        if regimen in {"", "unknown"}:
            ambiguities.append("missing_regimen")
        if ";" in regimen or regimen.count("/") >= 2:
            ambiguities.append("mixed_regimen_tokens")
        if clean_str(r["timepoint_normalized"]) == "unknown":
            ambiguities.append("unknown_timepoint")
        if "source_ambiguity" in clean_str(r["exclusion_reason"]):
            ambiguities.append("source_ambiguity")

        if "missing_regimen" in ambiguities or len(ambiguities) >= 3:
            confidence = "low"
        elif ambiguities:
            confidence = "medium"
        else:
            confidence = "high"

        flag_rows.append(
            {
                "cohort_id": r["cohort_id"],
                "patient_id": r["patient_id"],
                "sample_id": r["sample_id"],
                "treatment_regimen_raw": r["treatment_regimen_raw"],
                "treatment_context": r["treatment_context"],
                "confidence_level": confidence,
                "reason_for_assignment": reason_map.get(r["treatment_context"], "rule-based assignment"),
                "ambiguity_note": ";".join(sorted(set(ambiguities))),
            }
        )

    treatment_flags = pd.DataFrame(flag_rows)
    treatment_flags_path = OUT_DIR / "treatment_context_flags.csv"
    treatment_flags.to_csv(treatment_flags_path, index=False)

    # Role assignment
    role_rows = []
    for _, r in cohort_registry.iterrows():
        roles = []
        if is_true(r["usable_for_PD1_anchor"]):
            roles += ["main_scRNA_training", "PD1_anchor_training"]
        if is_true(r["usable_for_PD1X_extension"]):
            roles += ["main_scRNA_training", "PD1X_extension_support"]
        if is_true(r["usable_for_HCC_specific"]):
            roles.append("HCC_specific_deep_dive")
        if is_true(r["usable_for_spatial_adjudication"]):
            roles.append("spatial_adjudication")
        if is_true(r["bulk_available"]):
            roles.append("bulk_external_anchor")
        if is_true(r["perturb_available"]):
            roles.append("perturb_prior")
        if r["treatment_context"] == "high_confounding_support":
            roles.append("sensitivity_only")
        if not roles:
            roles = ["excluded"]

        ordered = []
        for role in [
            "main_scRNA_training",
            "PD1_anchor_training",
            "PD1X_extension_support",
            "HCC_specific_deep_dive",
            "spatial_adjudication",
            "bulk_external_anchor",
            "perturb_prior",
            "sensitivity_only",
            "excluded",
        ]:
            if role in roles and role not in ordered:
                ordered.append(role)

        role_rows.append(
            {
                "cohort_id": r["cohort_id"],
                "assigned_roles": ";".join(ordered),
                "primary_role": ordered[0],
                "rationale": r["notes"],
            }
        )

    role_df = pd.DataFrame(role_rows).sort_values("cohort_id").reset_index(drop=True)
    role_path = OUT_DIR / "dataset_role_assignment.csv"
    role_df.to_csv(role_path, index=False)

    # leakage risk log
    main_mask = sample_metadata["usable_in_main_analysis"] == "True"
    main_samples = sample_metadata[main_mask].copy()
    patient_multi = main_samples.groupby(["cohort_id", "patient_id"])["sample_id"].nunique().reset_index(name="n")
    n_multi = int((patient_multi["n"] > 1).sum())

    response_known = sample_metadata[sample_metadata["response_binary"].isin(["responder", "non_responder"])]
    dominance = []
    if not response_known.empty:
        by = response_known.groupby(["cohort_id", "response_binary"]).size().reset_index(name="n")
        tot = response_known.groupby("cohort_id").size().reset_index(name="total")
        m = by.merge(tot, on="cohort_id", how="left")
        for _, x in m.iterrows():
            dominance.append(f"- `{x['cohort_id']}` / `{x['response_binary']}`: {int(x['n'])}/{int(x['total'])} ({x['n']/x['total']:.2%})")
    else:
        dominance.append("- 暂无可用于 responder/non_responder 分析的标签。")

    bulk_patients = set(sample_metadata[(sample_metadata["bulk_available"] == "True") & (sample_metadata["patient_id"] != "unknown")]["patient_id"].tolist())
    scrna_patients = set(sample_metadata[(sample_metadata["scRNA_available"] == "True") & (sample_metadata["usable_in_main_analysis"] == "True") & (sample_metadata["patient_id"] != "unknown")]["patient_id"].tolist())

    post_only_cohorts = cohort_registry[
        cohort_registry["timepoint_schema"].str.lower().str.contains("post")
        & ~cohort_registry["timepoint_schema"].str.lower().str.contains("pre")
    ]["cohort_id"].tolist()

    mixed_context = cohort_registry[cohort_registry["notes"].str.contains("ambiguity", na=False)]["cohort_id"].tolist()

    risk_md = f"""# data_leakage_risk_log（v6.1 Step1）

## 1. 同一 patient 多样本跨 train/test 泄漏风险
- 主分析可用 scRNA 样本数：**{len(main_samples)}**
- 具有 >=2 个主分析样本的 patient 数：**{n_multi}**
- 建议：后续统一采用 patient-group split，禁止 sample-level 随机切分。

## 2. cohort/center 对 response 组别支配风险
{chr(10).join(dominance)}

## 3. pre/post 泄漏风险
- `has_paired_pre_post=True` 的 patient 数：**{int((patient_metadata['has_paired_pre_post']=='True').sum())}**
- 建议：同一 patient 的 pre/post 必须进入同一折。

## 4. bulk external anchor 与 scRNA training 重叠
- bulk patient 数：**{len(bulk_patients)}**
- 主分析 scRNA patient 数：**{len(scrna_patients)}**
- 字符串级 patient_id 交集：**{len(bulk_patients.intersection(scrna_patients))}**

## 5. response 受后续治疗影响风险
- post-only cohort：{', '.join(post_only_cohorts) if post_only_cohorts else 'none'}

## 6. therapy_context 混杂风险
- 规则标记为 ambiguity 的 cohort：{', '.join(sorted(mixed_context)) if mixed_context else 'none'}
"""
    risk_path = OUT_DIR / "data_leakage_risk_log.md"
    risk_path.write_text(risk_md, encoding="utf-8")

    # missing data request
    not_ingested = cohort_registry[cohort_registry["exclusion_reason"].str.contains("not_locally_ingested", na=False)]["cohort_id"].tolist()
    unknown_response_n = int((sample_metadata["response_binary"] == "unknown").sum())
    unknown_tp_n = int((sample_metadata["timepoint_normalized"] == "unknown").sum())

    missing_md = f"""# missing_data_request（v6.1 Step1）

## 概览
- registry cohort 数：**{len(cohort_registry)}**
- sample metadata 行数：**{len(sample_metadata)}**
- response_binary=unknown 行数：**{unknown_response_n}**
- timepoint_normalized=unknown 行数：**{unknown_tp_n}**

## 需要补齐（优先级）
1. 逐 sample 的治疗线次（line-of-therapy）与既往治疗时间关系。
2. `response` 为 NA/NE/UN 的 cohort 的原始临床对照表。
3. cohort 级混合治疗拆分到 sample 级的 mapping。

## 本地未入库/未标准化 cohort
{', '.join(not_ingested) if not_ingested else 'none'}

## Source 记录
- 本步骤本地 source 记录：`results/v6_1/step1/source_records_local.tsv`
- 联网 source 记录：如执行成功将写入 `results/v6_1/step1/source_records.tsv`
"""
    missing_path = OUT_DIR / "missing_data_request.md"
    missing_path.write_text(missing_md, encoding="utf-8")

    # summary
    context_counts = cohort_registry["treatment_context"].value_counts().to_dict()
    summary_md = f"""# step1_summary（v6.1，修订版）

## 本轮修订点
1. 修订 treatment_context 判定：先依据当前治疗方案（`treatment`），既往治疗（`prior_treatment`）仅作为 ambiguity/confound 标记。
2. 新增本地 h5ad `obs` 逐 sample 回填，显著降低 response `unknown`。
3. 在 sample 表新增 `response_standard` 列（如 mRECIST/RECIST/PATHOLOGICAL_RESPONSE）。
4. TASK01/TASK02 明确从 `/home/huyudi/006/data/processed/srt/raw/` 纳入。

## 执行范围
- 仅做数据治理（不做模型训练/机制解释/绘图）。

## 输出文件
1. cohort_registry_v6_1.csv
2. sample_metadata_master_v6_1.csv
3. patient_metadata_master_v6_1.csv
4. treatment_context_flags.csv
5. dataset_role_assignment.csv
6. data_leakage_risk_log.md
7. missing_data_request.md
8. source_records_local.tsv

## 关键统计
- cohort 数：**{len(cohort_registry)}**
- sample 行数：**{len(sample_metadata)}**
- patient 行数：**{len(patient_metadata)}**
- response_binary=unknown：**{unknown_response_n}**
- treatment_context 分布：**{json.dumps(context_counts, ensure_ascii=False)}**
"""
    summary_path = OUT_DIR / "step1_summary.md"
    summary_path.write_text(summary_md, encoding="utf-8")

    # local source records
    pd.DataFrame(local_source_records).drop_duplicates().to_csv(LOCAL_SOURCE_RECORDS, sep="\t", index=False)

    outputs = [
        str(cohort_registry_path.relative_to(ROOT)),
        str(sample_metadata_path.relative_to(ROOT)),
        str(patient_metadata_path.relative_to(ROOT)),
        str(treatment_flags_path.relative_to(ROOT)),
        str(role_path.relative_to(ROOT)),
        str(risk_path.relative_to(ROOT)),
        str(missing_path.relative_to(ROOT)),
        str(summary_path.relative_to(ROOT)),
        str(LOCAL_SOURCE_RECORDS.relative_to(ROOT)),
    ]
    known_missing = [
        f"not_locally_ingested_cohorts={len(not_ingested)}",
        f"unknown_response_rows={unknown_response_n}",
        f"unknown_timepoint_rows={unknown_tp_n}",
    ]
    write_manifest(outputs, known_missing)


if __name__ == "__main__":
    main()
