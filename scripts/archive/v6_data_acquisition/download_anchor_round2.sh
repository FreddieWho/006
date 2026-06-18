#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/huyudi/006"
GSE301="${ROOT}/data/external_anchor_candidates/GSE301741_hnscc_pembrolizumab_scRNA"
GSE286="${ROOT}/data/external_anchor_candidates/GSE286827_hnscc_durvalumab_scRNA"

mkdir -p "${GSE301}" "${GSE286}"

curl -L --fail --retry 5 -C - \
  -o "${GSE301}/GSE301741_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301741/suppl/GSE301741_RAW.tar"

curl -L --fail --retry 5 -C - \
  -o "${GSE301}/GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301741/suppl/GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds"

curl -L --fail --retry 5 -C - --output-dir "${GSE286}" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_Bcells.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_CD4T.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_CD8T.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_Malignant.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_TAMs.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_TumorSpecific_CD8T.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_immunecells_broad.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_Bcells.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_CD4T_TCR.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_CD8T_TCR.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_Malignant.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_TAMs.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_TumorSpecific_CD8T.rds.gz" \
  -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_immunecells_broad.rds.gz"
