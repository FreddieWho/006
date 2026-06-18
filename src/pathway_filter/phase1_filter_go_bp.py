"""
Phase 1: Baseline GO BP Filtering

Apply three filters to GO Biological Process pathways:
  1. Gene count: [MIN_GENES, MAX_GENES]
  2. Information Content (IC): >= IC_THRESHOLD
  3. Semantic deduplication: Remove redundant pathways (similarity > threshold)

Input:
  - data/pathway/go-basic.obo: GO ontology
  - data/pathway/go_bp_full.json: Full GO BP pathways

Output:
  - results/pathway_filter/phase1_baseline/go_bp_filtered.tsv: Filtered pathways (main)
  - results/pathway_filter/phase1_baseline/go_bp_filtered.json: Filtered pathways (JSON)
  - results/pathway_filter/phase1_baseline/filtering_summary.tsv: Filtering statistics
  - results/pathway_filter/phase1_baseline/removed_pathways.tsv: Removed pathways with reasons
  - results/pathway_filter/phase1_baseline/ic_scores.tsv: IC scores for all terms
  - results/pathway_filter/phase1_baseline/gene_universe.txt: All unique genes

References:
  [1] Subramanian et al. PNAS 2005
      "Gene Set Enrichment Analysis: A Knowledge-Based Approach"
  [2] Yu et al. Bioinformatics 2010
      "clusterProfiler: an R package for comparing biological themes among gene clusters"
  [3] Pesquita et al. PLoS Comput Biol 2009
      "Semantic Similarity in Biomedical Ontologies"
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import logging
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Any
import numpy as np
import pandas as pd
from tqdm import tqdm

# GO processing
import networkx as nx

# Import configurations
from config.pwy_filter_common import (
    GO_OBO_FILE, GO_BP_JSON, setup_logger
)
from config.pwy_filter_phase1_filter import (
    MIN_GENES, MAX_GENES,
    IC_THRESHOLD,
    SEMANTIC_SIMILARITY_THRESHOLD,
    SIMILARITY_METHOD,
    GENE_FIELD,  
    GO_ID_FIELD,  
    PHASE1_DIR,
    FILTERED_PATHWAYS_TSV,
    FILTERED_PATHWAYS_JSON,
    FILTERING_SUMMARY,
    REMOVED_PATHWAYS,
    IC_SCORES,
    GENE_UNIVERSE,
)

# Import utilities
from src.pathway_filter.utils.io import (
    load_go_obo,
    load_go_pathways,
    save_pathways_tsv,
    save_pathways_json,
)

# Setup logger
logger = setup_logger('phase1_filter')


class GOPathwayFilter:
    """
    Filter GO BP pathways based on baseline criteria
    
    Implements three sequential filters with literature-backed thresholds
    """
    
    def __init__(self):
        """Initialize filter"""
        self.go_graph: nx.MultiDiGraph = None
        self.pathways: Dict[str, Dict[str, Any]] = None
        self.ic_scores: Dict[str, float] = {}
        self.removed: Dict[str, Dict[str, Any]] = {}
        
        logger.info("="*80)
        logger.info("Phase 1: Baseline GO BP Filtering")
        logger.info("="*80)
    
    def load_data(self):
        """Load GO ontology and pathways"""
        logger.info("\n" + "-"*80)
        logger.info("Loading input data")
        logger.info("-"*80)
        
        # Load GO ontology
        self.go_graph = load_go_obo(GO_OBO_FILE)
        
        # Load GO BP pathways
        raw_pathways = load_go_pathways(GO_BP_JSON)
        
        # ⭐ Normalize pathway structure
        logger.info("Normalizing pathway data structure...")
        self.pathways = {}
        for pathway_id, pathway_info in raw_pathways.items():
            # Extract genes based on actual field name
            if isinstance(pathway_info, dict):
                genes = pathway_info.get(GENE_FIELD, [])
            elif isinstance(pathway_info, list):
                genes = pathway_info
            else:
                logger.warning(f"Unexpected format for {pathway_id}: {type(pathway_info)}")
                genes = []
            
            # Create normalized structure
            self.pathways[pathway_id] = {
                'name': pathway_info.get('systematicName', pathway_id) if isinstance(pathway_info, dict) else pathway_id,
                'genes': genes,
                'go_id': pathway_info.get(GO_ID_FIELD, '') if isinstance(pathway_info, dict) else '',
                'original_data': pathway_info if isinstance(pathway_info, dict) else {}
            }
        
        initial_count = len(self.pathways)
        
        # Log gene count distribution
        gene_counts = [len(p['genes']) for p in self.pathways.values()]
        logger.info(f"Gene count distribution:")
        logger.info(f"  Min: {min(gene_counts)}")
        logger.info(f"  Max: {max(gene_counts)}")
        logger.info(f"  Mean: {np.mean(gene_counts):.1f}")
        logger.info(f"  Median: {np.median(gene_counts):.1f}")
        
        logger.info(f"\nInitial pathway count: {initial_count}")
        
        return initial_count
    
    def filter_gene_count(self) -> int:
        """
        Filter 1: Gene count [MIN_GENES, MAX_GENES]
        
        Reference: Subramanian et al. PNAS 2005 [1]
        Reference: Liberzon et al. Cell Syst 2015
        
        Returns:
            Number of pathways retained
        """
        logger.info("\n" + "-"*80)
        logger.info(f"Filter 1: Gene count [{MIN_GENES}, {MAX_GENES}]")
        logger.info(f"Reference: Subramanian et al. PNAS 2005")
        logger.info("-"*80)
        
        filtered = {}
        removed_too_small = 0
        removed_too_large = 0
        
        for pathway_id, pathway_info in self.pathways.items():
            genes = pathway_info.get('genes', [])
            gene_count = len(genes)
            
            if gene_count < MIN_GENES:
                removed_too_small += 1
                self.removed[pathway_id] = {
                    **pathway_info,
                    'filter': 'gene_count',
                    'reason': f'Too few genes ({gene_count} < {MIN_GENES})'
                }
            elif gene_count > MAX_GENES:
                removed_too_large += 1
                self.removed[pathway_id] = {
                    **pathway_info,
                    'filter': 'gene_count',
                    'reason': f'Too many genes ({gene_count} > {MAX_GENES})'
                }
            else:
                filtered[pathway_id] = pathway_info
        
        logger.info(f"  Retained: {len(filtered)} pathways")
        logger.info(f"  Removed: {removed_too_small + removed_too_large} pathways")
        logger.info(f"    - Too small (< {MIN_GENES}): {removed_too_small}")
        logger.info(f"    - Too large (> {MAX_GENES}): {removed_too_large}")
        
        self.pathways = filtered
        return len(filtered)
    
    def compute_ic_scores(self):
        """
        Compute Information Content (IC) for each GO term
        
        Reference: Pesquita et al. PLoS Comput Biol 2009 [3]
        """
        logger.info("\n" + "-"*80)
        logger.info("Computing Information Content (IC) scores")
        logger.info("Reference: Pesquita et al. PLoS Comput Biol 2009")
        logger.info("-"*80)
        
        # ⭐ Check if we have any pathways left
        if len(self.pathways) == 0:
            logger.warning("No pathways left to compute IC scores!")
            return {}
        
        # Collect all genes annotated to each term
        term_to_genes = defaultdict(set)
        
        logger.info("  Propagating gene annotations to ancestors...")
        for pathway_id, pathway_info in tqdm(self.pathways.items(), 
                                             desc="  Processing terms"):
            genes = set(pathway_info.get('genes', []))
            
            # Add genes to this term
            term_to_genes[pathway_id].update(genes)
            
            # Get GO ID for this pathway
            go_id = pathway_info.get('go_id', '')
            
            # Propagate to ancestors if GO ID is available
            if go_id and go_id in self.go_graph:
                try:
                    ancestors = nx.ancestors(self.go_graph, go_id)
                    # Map ancestors back to pathway IDs
                    for ancestor_go_id in ancestors:
                        # Find pathway with this GO ID
                        for pid, pinfo in self.pathways.items():
                            if pinfo.get('go_id', '') == ancestor_go_id:
                                term_to_genes[pid].update(genes)
                except nx.NetworkXError:
                    pass
        
        # Get total unique genes
        all_genes = set()
        for genes in term_to_genes.values():
            all_genes.update(genes)
        total_genes = len(all_genes)
        
        logger.info(f"  Total unique genes in GO BP: {total_genes}")
        
        # ⭐ Handle edge case
        if total_genes == 0:
            logger.warning("No genes found! Using gene count as proxy for IC")
            # Use inverse of gene count as IC proxy
            for term_id, pathway_info in self.pathways.items():
                n_genes = len(pathway_info.get('genes', []))
                self.ic_scores[term_id] = np.log(500 / (n_genes + 1))  # Arbitrary scaling
        else:
            # Compute IC for each term
            logger.info("  Computing IC scores...")
            for term_id in self.pathways.keys():
                genes = term_to_genes.get(term_id, set())
                
                freq = len(genes) / total_genes
                ic = -np.log(freq + 1e-10)
                
                self.ic_scores[term_id] = ic
        
        ic_values = list(self.ic_scores.values())
        
        # ⭐ Check if we have IC scores
        if len(ic_values) == 0:
            logger.warning("No IC scores computed!")
            return {}
        
        logger.info(f"  IC score range: [{min(ic_values):.2f}, {max(ic_values):.2f}]")
        logger.info(f"  IC score mean: {np.mean(ic_values):.2f}")
        logger.info(f"  IC score median: {np.median(ic_values):.2f}")
        
        # Save IC scores
        ic_df = pd.DataFrame([
            {'pathway_id': pid, 'ic_score': ic}
            for pid, ic in self.ic_scores.items()
        ])
        ic_df = ic_df.sort_values('ic_score', ascending=False)
        ic_df.to_csv(IC_SCORES, sep='\t', index=False)
        logger.info(f"  Saved IC scores to {IC_SCORES}")
        
        return self.ic_scores
    
    def filter_ic(self) -> int:
        """
        Filter 2: IC >= IC_THRESHOLD
        
        Reference: Yu et al. Bioinformatics 2010 [2]
        
        Returns:
            Number of pathways retained
        """
        logger.info("\n" + "-"*80)
        logger.info(f"Filter 2: IC >= {IC_THRESHOLD}")
        logger.info("Reference: Yu et al. Bioinformatics 2010")
        logger.info("-"*80)
        
        if not self.ic_scores:
            self.compute_ic_scores()
        
        # ⭐ Handle case with no IC scores
        if not self.ic_scores:
            logger.warning("No IC scores available, skipping IC filter")
            return len(self.pathways)
        
        filtered = {}
        removed_count = 0
        
        for pathway_id, pathway_info in self.pathways.items():
            ic = self.ic_scores.get(pathway_id, 0)
            
            if ic >= IC_THRESHOLD:
                pathway_info['ic_score'] = ic
                filtered[pathway_id] = pathway_info
            else:
                removed_count += 1
                self.removed[pathway_id] = {
                    **pathway_info,
                    'ic_score': ic,
                    'filter': 'ic_threshold',
                    'reason': f'IC too low ({ic:.2f} < {IC_THRESHOLD})'
                }
        
        logger.info(f"  Retained: {len(filtered)} pathways")
        logger.info(f"  Removed: {removed_count} pathways (IC < {IC_THRESHOLD})")
        
        self.pathways = filtered
        return len(filtered)
    
    def compute_semantic_similarity(self, term1: str, term2: str) -> float:
        """
        Compute semantic similarity using Jaccard similarity
        
        Reference: Jiang & Conrath, ROCLING 1997
        
        Args:
            term1: Pathway ID
            term2: Pathway ID
        
        Returns:
            Similarity score [0, 1]
        """
        genes1 = set(self.pathways[term1].get('genes', []))
        genes2 = set(self.pathways[term2].get('genes', []))
        
        if len(genes1) == 0 or len(genes2) == 0:
            return 0.0
        
        intersection = len(genes1 & genes2)
        union = len(genes1 | genes2)
        
        if SIMILARITY_METHOD == 'jaccard':
            similarity = intersection / union if union > 0 else 0.0
        elif SIMILARITY_METHOD == 'overlap':
            similarity = intersection / min(len(genes1), len(genes2))
        elif SIMILARITY_METHOD == 'dice':
            similarity = 2 * intersection / (len(genes1) + len(genes2))
        else:
            similarity = intersection / union if union > 0 else 0.0
        
        return similarity
    
    def semantic_deduplication(self) -> int:
        """
        Filter 3: Remove redundant pathways
        
        Reference: Yu et al. Bioinformatics 2010 [2]
        
        Returns:
            Number of pathways retained
        """
        logger.info("\n" + "-"*80)
        logger.info(f"Filter 3: Semantic deduplication (similarity > {SEMANTIC_SIMILARITY_THRESHOLD})")
        logger.info(f"Method: {SIMILARITY_METHOD.capitalize()} similarity")
        logger.info("Reference: Yu et al. Bioinformatics 2010")
        logger.info("-"*80)
        
        pathway_ids = list(self.pathways.keys())
        n_pathways = len(pathway_ids)
        
        # ⭐ Handle edge cases
        if n_pathways == 0:
            logger.warning("No pathways to deduplicate")
            return 0
        
        if n_pathways == 1:
            logger.info("Only 1 pathway, skipping deduplication")
            return 1
        
        logger.info(f"  Computing pairwise similarities for {n_pathways} pathways...")
        logger.info(f"  Total comparisons: {n_pathways * (n_pathways - 1) // 2:,}")
        
        to_remove = set()
        redundancy_info = []
        
        with tqdm(total=n_pathways, desc="  Deduplicating") as pbar:
            for i in range(n_pathways):
                if pathway_ids[i] in to_remove:
                    pbar.update(1)
                    continue
                
                for j in range(i+1, n_pathways):
                    if pathway_ids[j] in to_remove:
                        continue
                    
                    sim = self.compute_semantic_similarity(pathway_ids[i], pathway_ids[j])
                    
                    if sim > SEMANTIC_SIMILARITY_THRESHOLD:
                        ic_i = self.ic_scores.get(pathway_ids[i], 0)
                        ic_j = self.ic_scores.get(pathway_ids[j], 0)
                        
                        if ic_i >= ic_j:
                            removed_id = pathway_ids[j]
                            kept_id = pathway_ids[i]
                        else:
                            removed_id = pathway_ids[i]
                            kept_id = pathway_ids[i]
                        
                        to_remove.add(removed_id)
                        
                        redundancy_info.append({
                            'removed_id': removed_id,
                            'removed_name': self.pathways[removed_id].get('name', 'Unknown'),
                            'kept_id': kept_id,
                            'kept_name': self.pathways[kept_id].get('name', 'Unknown'),
                            'similarity': sim,
                            'removed_ic': self.ic_scores.get(removed_id, 0),
                            'kept_ic': self.ic_scores.get(kept_id, 0),
                        })
                        
                        if removed_id == pathway_ids[i]:
                            break
                
                pbar.update(1)
        
        # Remove redundant pathways
        filtered = {}
        for pid, info in self.pathways.items():
            if pid in to_remove:
                redundant_with = None
                for r in redundancy_info:
                    if r['removed_id'] == pid:
                        redundant_with = r['kept_id']
                        break
                
                self.removed[pid] = {
                    **info,
                    'filter': 'semantic_deduplication',
                    'reason': f'Redundant with {redundant_with} (similarity > {SEMANTIC_SIMILARITY_THRESHOLD})'
                }
            else:
                filtered[pid] = info
        
        logger.info(f"  Retained: {len(filtered)} pathways")
        logger.info(f"  Removed: {len(to_remove)} redundant pathways")
        
        if redundancy_info:
            redundancy_df = pd.DataFrame(redundancy_info)
            redundancy_file = PHASE1_DIR / "redundancy_details.tsv"
            redundancy_df.to_csv(redundancy_file, sep='\t', index=False)
            logger.info(f"  Saved redundancy details to {redundancy_file}")
        
        self.pathways = filtered
        return len(filtered)
    
    def save_results(self):
        """Save filtered pathways and summary"""
        logger.info("\n" + "-"*80)
        logger.info("Saving results")
        logger.info("-"*80)
        
        # Add IC scores
        for pathway_id, pathway_info in self.pathways.items():
            if 'ic_score' not in pathway_info:
                pathway_info['ic_score'] = self.ic_scores.get(pathway_id, 0)
        
        # Save filtered pathways
        save_pathways_tsv(self.pathways, FILTERED_PATHWAYS_TSV, include_genes=True)
        logger.info(f"  Main output: {FILTERED_PATHWAYS_TSV}")
        
        save_pathways_json(self.pathways, FILTERED_PATHWAYS_JSON)
        logger.info(f"  JSON output: {FILTERED_PATHWAYS_JSON}")
        
        # Save removed pathways
        if self.removed:
            removed_records = []
            for pid, info in self.removed.items():
                removed_records.append({
                    'pathway_id': pid,
                    'name': info.get('name', 'Unknown'),
                    'n_genes': len(info.get('genes', [])),
                    'ic_score': info.get('ic_score', 0),
                    'filter': info.get('filter', 'unknown'),
                    'reason': info.get('reason', 'unknown'),
                })
            
            removed_df = pd.DataFrame(removed_records)
            removed_df.to_csv(REMOVED_PATHWAYS, sep='\t', index=False)
            logger.info(f"  Removed pathways: {REMOVED_PATHWAYS}")
        
        # Filtering summary
        initial_count = len(self.pathways) + len(self.removed)
        
        filter_counts = defaultdict(int)
        for info in self.removed.values():
            filter_counts[info.get('filter', 'unknown')] += 1
        
        summary_data = {
            'stage': [
                'initial',
                'after_gene_count',
                'after_ic_threshold',
                'after_deduplication',
            ],
            'n_pathways': [
                initial_count,
                initial_count - filter_counts.get('gene_count', 0),
                initial_count - filter_counts.get('gene_count', 0) - filter_counts.get('ic_threshold', 0),
                len(self.pathways),
            ],
            'removed': [
                0,
                filter_counts.get('gene_count', 0),
                filter_counts.get('ic_threshold', 0),
                filter_counts.get('semantic_deduplication', 0),
            ]
        }
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(FILTERING_SUMMARY, sep='\t', index=False)
        logger.info(f"  Filtering summary: {FILTERING_SUMMARY}")
        
        # Gene universe
        all_genes = set()
        for info in self.pathways.values():
            all_genes.update(info.get('genes', []))
        
        if all_genes:
            with open(GENE_UNIVERSE, 'w') as f:
                for gene in sorted(all_genes):
                    f.write(gene + '\n')
            logger.info(f"  Gene universe ({len(all_genes)} genes): {GENE_UNIVERSE}")
    
    def print_summary(self):
        """Print final summary"""
        logger.info("\n" + "="*80)
        logger.info("Phase 1 Summary")
        logger.info("="*80)
        
        logger.info(f"\nFinal pathway count: {len(self.pathways)}")
        
        if len(self.pathways) > 0:
            # Gene count distribution
            gene_counts = [len(info.get('genes', [])) for info in self.pathways.values()]
            logger.info(f"\nGene count statistics:")
            logger.info(f"  Min: {min(gene_counts)}")
            logger.info(f"  Max: {max(gene_counts)}")
            logger.info(f"  Mean: {np.mean(gene_counts):.1f}")
            logger.info(f"  Median: {np.median(gene_counts):.1f}")
            
            # IC score distribution
            ic_values = [info.get('ic_score', 0) for info in self.pathways.values()]
            logger.info(f"\nIC score statistics:")
            logger.info(f"  Min: {min(ic_values):.2f}")
            logger.info(f"  Max: {max(ic_values):.2f}")
            logger.info(f"  Mean: {np.mean(ic_values):.2f}")
            logger.info(f"  Median: {np.median(ic_values):.2f}")
            
            # Total unique genes
            all_genes = set()
            for info in self.pathways.values():
                all_genes.update(info.get('genes', []))
            logger.info(f"\nTotal unique genes: {len(all_genes)}")
        
        logger.info("\n" + "="*80)


def main():
    """Main function for Phase 1"""
    try:
        filter_obj = GOPathwayFilter()
        
        filter_obj.load_data()
        filter_obj.filter_gene_count()
        filter_obj.compute_ic_scores()
        filter_obj.filter_ic()
        filter_obj.semantic_deduplication()
        filter_obj.save_results()
        filter_obj.print_summary()
        
        logger.info("\n" + "="*80)
        logger.info("Phase 1 completed successfully!")
        logger.info(f"Output directory: {PHASE1_DIR}")
        logger.info("="*80)
        
        return 0
        
    except Exception as e:
        logger.error(f"\nError in Phase 1: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    exit(main())