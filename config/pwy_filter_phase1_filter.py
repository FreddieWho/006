"""
Configuration for Phase 1: Baseline GO BP Filtering

Parameters with literature references for:
- Filter 1: Gene count
- Filter 2: Information Content (IC)
- Filter 3: Semantic similarity (deduplication)
"""

from pathlib import Path

# Import common settings
import sys
sys.path.append(str(Path(__file__).parent.parent))
from config.pwy_filter_common import RESULTS_DIR

# ============================================================================
# JSON Field Mapping
# ============================================================================
# User's go_bp_full.json uses different field names
GENE_FIELD = 'geneSymbols'  # Field name for gene list
NAME_FIELD = 'systematicName'  # Field name for pathway name (optional)
GO_ID_FIELD = 'exactSource'  # Field name for GO ID

# ============================================================================
# Output Directories
# ============================================================================
PHASE1_DIR = RESULTS_DIR / "phase1_baseline"
PHASE1_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Filter 1: Gene Count
# ============================================================================
# Reference: Subramanian et al. PNAS 2005
#   "Gene Set Enrichment Analysis: A Knowledge-Based Approach"
#   Original text: "We considered gene sets containing between 15 and 500 genes"
#
# Reference: Liberzon et al. Cell Syst 2015
#   "The Molecular Signatures Database (MSigDB) Hallmark Gene Set Collection"
#   MSigDB guideline: 15-300 genes for focused analysis

MIN_GENES = 15
MAX_GENES = 300

# ============================================================================
# Filter 2: Information Content (IC) Threshold
# ============================================================================
# Reference: Yu et al. Bioinformatics 2010
#   "clusterProfiler: an R package for comparing biological themes"
#   GOSemSim package uses IC to measure term specificity
#
# Reference: Pesquita et al. PLoS Comput Biol 2009
#   "Semantic Similarity in Biomedical Ontologies"
#   IC = -log(P(term)), where P is frequency
#   IC > 4 means term annotates < 1.8% of genes (e^-4 ≈ 0.018)
#   This ensures term specificity

IC_THRESHOLD = 4.0

# ============================================================================
# Filter 3: Semantic Similarity (for deduplication)
# ============================================================================
# Reference: Yu et al. Bioinformatics 2010
#   Original implementation in simplify() function
#   Quote: "cutoff = 0.7 was used to remove redundant GO terms"
#   Terms with similarity > 0.7 are considered redundant
#   When redundant, keep the one with higher IC (more specific)

SEMANTIC_SIMILARITY_THRESHOLD = 0.7

# ============================================================================
# Similarity Calculation Method
# ============================================================================
# Method: Jaccard similarity of gene sets
# Reference: Jiang & Conrath, ROCLING 1997 (semantic similarity measures)
# Jaccard = |A ∩ B| / |A ∪ B|
# Widely used for gene set comparison due to simplicity and interpretability

SIMILARITY_METHOD = 'jaccard'  # Options: 'jaccard', 'overlap', 'dice'

# ============================================================================
# Output Files
# ============================================================================
# Main outputs (TSV format for cross-language compatibility)
FILTERED_PATHWAYS_TSV = PHASE1_DIR / "go_bp_filtered.tsv"
FILTERING_SUMMARY = PHASE1_DIR / "filtering_summary.tsv"
REMOVED_PATHWAYS = PHASE1_DIR / "removed_pathways.tsv"
IC_SCORES = PHASE1_DIR / "ic_scores.tsv"

# Auxiliary outputs (for downstream phases)
FILTERED_PATHWAYS_JSON = PHASE1_DIR / "go_bp_filtered.json"  # JSON for easy loading
GENE_UNIVERSE = PHASE1_DIR / "gene_universe.txt"  # All unique genes