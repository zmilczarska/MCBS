#!/bin/bash
set -euo pipefail

# Args: path to annotation BED file
bedfile="$1"
cd "$HOME/modelowanie"

# Validate input file
if [ ! -f "$bedfile" ]; then
    echo "[ERROR] BED file not found: $bedfile" >&2
    exit 1
fi

celltype_name=$(basename "$bedfile" .bed)
outdir="TDEP-sLDSC/results/ldsc_${celltype_name}"

mkdir -p "$outdir"

echo "[INFO] Starting LDSC for cell type: $celltype_name"
echo "[INFO] BED file: $(realpath "$bedfile")"
echo "[INFO] Output directory: $(realpath "$outdir")"
echo "[INFO] Preview BED file:"
head "$bedfile"

# Sequentially compute LD scores for chromosomes 1-22
for chr in {1..22}; do
    echo "[INFO] Processing chromosome $chr"
    bash TDEP-sLDSC/scripts/compute_chrom_ldscores.sh "$(realpath "$bedfile")" "$(realpath "$outdir")" "$chr"
done

echo "[INFO] Finished LDSC for $celltype_name"
