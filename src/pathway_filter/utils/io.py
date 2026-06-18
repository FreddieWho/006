"""
Input/Output utility functions

Handles loading and saving of pathway data in various formats:
- GO OBO format
- JSON format
- TSV format
- GMT format (for GSEA/GSVA)
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd

# GO processing
try:
    import obonet
    import networkx as nx
except ImportError:
    raise ImportError("Please install: pip install obonet networkx")

logger = logging.getLogger(__name__)


def load_go_obo(obo_file: Path) -> nx.MultiDiGraph:
    """
    Load GO ontology from OBO format
    
    Args:
        obo_file: Path to go-basic.obo file
    
    Returns:
        nx.MultiDiGraph: GO ontology graph
    
    Reference:
        Gene Ontology Consortium, Nucleic Acids Res 2021
        OBO format specification: http://www.geneontology.org/docs/
    """
    logger.info(f"Loading GO ontology from {obo_file}")
    
    if not obo_file.exists():
        raise FileNotFoundError(f"GO OBO file not found: {obo_file}")
    
    graph = obonet.read_obo(obo_file)
    logger.info(f"  Loaded {len(graph)} GO terms")
    
    return graph


def load_go_pathways(json_file: Path) -> Dict[str, Dict[str, Any]]:
    """
    Load GO BP pathways from JSON format
    
    Args:
        json_file: Path to go_bp_full.json
    
    Returns:
        dict: {pathway_id: {name, genes, ...}}
    
    Expected JSON structure:
        {
            "GO:0000001": {
                "name": "pathway name",
                "genes": ["GENE1", "GENE2", ...],
                ...
            },
            ...
        }
    """
    logger.info(f"Loading GO BP pathways from {json_file}")
    
    if not json_file.exists():
        raise FileNotFoundError(f"GO BP JSON file not found: {json_file}")
    
    with open(json_file, 'r', encoding='utf-8') as f:
        pathways = json.load(f)
    
    logger.info(f"  Loaded {len(pathways)} GO BP pathways")
    
    return pathways


def load_hallmark(json_file: Path) -> Dict[str, Dict[str, Any]]:
    """
    Load Hallmark gene sets from JSON format
    
    Args:
        json_file: Path to hallmark.json
    
    Returns:
        dict: {pathway_id: {name, genes, ...}}
    
    Reference:
        Liberzon et al. Cell Syst 2015
        "The Molecular Signatures Database Hallmark Gene Set Collection"
    """
    logger.info(f"Loading Hallmark gene sets from {json_file}")
    
    if not json_file.exists():
        raise FileNotFoundError(f"Hallmark JSON file not found: {json_file}")
    
    with open(json_file, 'r', encoding='utf-8') as f:
        hallmark = json.load(f)
    
    logger.info(f"  Loaded {len(hallmark)} Hallmark gene sets")
    
    return hallmark


def save_pathways_tsv(pathways: Dict[str, Dict[str, Any]], 
                      output_file: Path,
                      include_genes: bool = False):
    """
    Save pathways to TSV format
    
    Args:
        pathways: Dictionary of pathways
        output_file: Output TSV file path
        include_genes: If True, include gene list as semicolon-separated string
    
    Output columns:
        - pathway_id
        - name
        - n_genes
        - genes (optional, semicolon-separated)
        - ic_score (if available)
        - ... (other metadata)
    """
    logger.info(f"Saving pathways to {output_file}")
    
    records = []
    for pathway_id, info in pathways.items():
        record = {
            'pathway_id': pathway_id,
            'name': info.get('name', 'Unknown'),
            'n_genes': len(info.get('genes', [])),
        }
        
        # Add genes if requested
        if include_genes:
            genes = info.get('genes', [])
            record['genes'] = ';'.join(genes)
        
        # Add other metadata
        for key in ['ic_score', 'detectability', 'reproducibility', 'source']:
            if key in info:
                record[key] = info[key]
        
        records.append(record)
    
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  Saved {len(records)} pathways")


def save_pathways_json(pathways: Dict[str, Dict[str, Any]], 
                       output_file: Path):
    """
    Save pathways to JSON format
    
    Args:
        pathways: Dictionary of pathways
        output_file: Output JSON file path
    """
    logger.info(f"Saving pathways to {output_file}")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(pathways, f, indent=2, ensure_ascii=False)
    
    logger.info(f"  Saved {len(pathways)} pathways")


def save_pathways_gmt(pathways: Dict[str, Dict[str, Any]], 
                      output_file: Path):
    """
    Save pathways to GMT format (for GSEA/GSVA)
    
    Args:
        pathways: Dictionary of pathways
        output_file: Output GMT file path
    
    GMT format:
        pathway_name\tdescription\tGENE1\tGENE2\t...
    
    Reference:
        Subramanian et al. PNAS 2005
        GMT format specification: 
        https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    """
    logger.info(f"Saving pathways to GMT format: {output_file}")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        for pathway_id, info in pathways.items():
            name = info.get('name', pathway_id)
            genes = info.get('genes', [])
            
            # GMT format: name, description, genes...
            line = f"{name}\t{pathway_id}\t" + '\t'.join(genes) + '\n'
            f.write(line)
    
    logger.info(f"  Saved {len(pathways)} pathways")


def load_expression_matrix(file_path: Path, 
                          sep: str = '\t',
                          index_col: int = 0,
                          comment: str = None) -> pd.DataFrame:
    """
    Load expression matrix from file
    
    Args:
        file_path: Path to expression file
        sep: Separator (default: tab)
        index_col: Index column (default: 0, gene names)
        comment: Comment character to skip
    
    Returns:
        pd.DataFrame: Expression matrix (genes x samples)
    """
    logger.info(f"Loading expression matrix from {file_path}")
    
    if not file_path.exists():
        raise FileNotFoundError(f"Expression file not found: {file_path}")
    
    df = pd.read_csv(file_path, sep=sep, index_col=index_col, comment=comment)
    logger.info(f"  Loaded matrix: {df.shape[0]} genes x {df.shape[1]} samples")
    
    return df


def load_metadata(file_path: Path, 
                 sep: str = '\t') -> pd.DataFrame:
    """
    Load sample metadata from file
    
    Args:
        file_path: Path to metadata file
        sep: Separator (default: tab)
    
    Returns:
        pd.DataFrame: Metadata
    """
    logger.info(f"Loading metadata from {file_path}")
    
    if not file_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {file_path}")
    
    df = pd.read_csv(file_path, sep=sep)
    logger.info(f"  Loaded {len(df)} samples")
    
    return df