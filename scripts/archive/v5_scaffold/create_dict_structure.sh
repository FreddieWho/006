#!/bin/bash

echo "Creating project directories..."

# 创建目录
mkdir -p data/raw/scRNA_seq
mkdir -p data/raw/clinical
mkdir -p data/processed
mkdir -p data/pathways
mkdir -p data/features/pathway_scores
mkdir -p data/features/cell_embeddings
mkdir -p data/graphs/cell_graphs
mkdir -p scripts/00_setup
mkdir -p scripts/03_pathway_filtering
mkdir -p src/models
mkdir -p experiments/configs
mkdir -p notebooks
mkdir -p results/models/checkpoints
mkdir -p figures/main
mkdir -p reports/pathway_filtering
mkdir -p logs/training
mkdir -p docs

echo "✓ Directories created!"

# 创建__init__.py
touch src/__init__.py
touch src/models/__init__.py
touch scripts/__init__.py

echo "✓ __init__.py files created!"

# 创建README.md
cat > README.md << 'EOF'
# PD-1 Response Prediction

Single-cell RNA-seq analysis for immunotherapy response prediction.
 
