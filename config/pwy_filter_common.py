"""
Common configuration for pathway filter pipeline

Shared paths, directories, and global settings used across all phases.
"""

from pathlib import Path
import logging
from datetime import datetime

# ============================================================================
# Project Structure
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
SRC_DIR = PROJECT_ROOT / "src" / "pathway_filter"
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_filter"
LOGS_DIR = PROJECT_ROOT / "logs" / "pathway_filter"

# Input data directories
PATHWAY_DIR = DATA_DIR / "pathway"
PUB_DATA_DIR = DATA_DIR / "pub"

# Create all directories
for dir_path in [RESULTS_DIR, LOGS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Input Files (User Provided)
# ============================================================================
# GO ontology and pathways
GO_OBO_FILE = PATHWAY_DIR / "go-basic.obo"
GO_BP_JSON = PATHWAY_DIR / "go_bp_full.json"
HALLMARK_JSON = PATHWAY_DIR / "hallmark.json"

# TCGA data (manually downloaded)
TCGA_DIR = PUB_DATA_DIR / "tcga"
TCGA_EXPRESSION = TCGA_DIR / "tcga_expression.tsv"
TCGA_METADATA = TCGA_DIR / "tcga_metadata.tsv"

# ICGC data (manually downloaded)
ICGC_DIR = PUB_DATA_DIR / "icgc"
ICGC_EXPRESSION = ICGC_DIR / "icgc_all_projects_donors_exp.tsv"
ICGC_PHENOTYPE = ICGC_DIR / "icgc_all_projects_donors_phenotype.tsv"

# ============================================================================
# Logging Configuration
# ============================================================================
def setup_logger(phase_name, level=logging.INFO):
    """
    Setup logger for a specific phase
    
    Args:
        phase_name: Name of the phase (e.g., 'phase1')
        level: Logging level (default: INFO)
    
    Returns:
        logger: Configured logger
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOGS_DIR / f"{phase_name}_{timestamp}.log"
    
    # Create logger
    logger = logging.getLogger(phase_name)
    logger.setLevel(level)
    logger.handlers = []  # Clear existing handlers
    
    # File handler
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(level)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    
    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

# ============================================================================
# TCGA Configuration
# ============================================================================
# Reference: Liu et al. Cell 2018
#   "An Integrated TCGA Pan-Cancer Clinical Data Resource"
# Selected solid tumors (excluding hematologic malignancies)
TCGA_SOLID_TUMORS = [
    'BLCA',  # Bladder Urothelial Carcinoma
    'BRCA',  # Breast invasive carcinoma
    'CESC',  # Cervical squamous cell carcinoma
    'CHOL',  # Cholangiocarcinoma
    'COAD',  # Colon adenocarcinoma
    'ESCA',  # Esophageal carcinoma
    'GBM',   # Glioblastoma multiforme
    'HNSC',  # Head and Neck squamous cell carcinoma
    'KICH',  # Kidney Chromophobe
    'KIRC',  # Kidney renal clear cell carcinoma
    'KIRP',  # Kidney renal papillary cell carcinoma
    'LIHC',  # Liver hepatocellular carcinoma
    'LUAD',  # Lung adenocarcinoma
    'LUSC',  # Lung squamous cell carcinoma
    'PAAD',  # Pancreatic adenocarcinoma
    'PRAD',  # Prostate adenocarcinoma
    'READ',  # Rectum adenocarcinoma
    'SKCM',  # Skin Cutaneous Melanoma
    'STAD',  # Stomach adenocarcinoma
    'UCEC',  # Uterine Corpus Endometrial Carcinoma
]

# ============================================================================
# ICGC Configuration
# ============================================================================
# Reference: ICGC/TCGA Pan-Cancer Analysis of Whole Genomes, Nature 2020
# Projects available in user's manual download
ICGC_PROJECTS = {
    'BRCA-UK': 'BRCA',
    'COAD-US': 'COAD',
    'LIRI-JP': 'LIHC',
    'PACA-AU': 'PAAD',
    'PRAD-CA': 'PRAD',
    'PRAD-UK': 'PRAD',
    # Add more as user provides
}

# ============================================================================
# General Settings
# ============================================================================
# Random seed for reproducibility
RANDOM_SEED = 25

# Number of parallel jobs
N_JOBS = 20  # Use all available cores

# Figure settings
FIGURE_DPI = 300
FIGURE_FORMAT = 'pdf'
