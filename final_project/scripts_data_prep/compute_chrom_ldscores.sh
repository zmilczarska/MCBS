#!/bin/bash
#SBATCH --job-name=ldsc_chr
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=8:00:00
#SBATCH --output=TDEP-sLDSC/results/compute_chrom_ldscores_%j.out
#SBATCH --error=TDEP-sLDSC/results/compute_chrom_ldscores_%j.err

set -euo pipefail

# Arguments:
# $1 = BED file
# $2 = output directory
# $3 = chromosome number

bedfile="$1"
workdir="$2"
chr="$3"
cd "$HOME/modelowanie"

# Reference files
plink_ref="TDEP-sLDSC/data/sldsc_ref/1000G_EUR_Phase3_plink/1000G.EUR.QC"
hm_snp="TDEP-sLDSC/data/sldsc_ref/hm_snp.txt"

echo "[INFO] Using BED file: $bedfile"
echo "[INFO] Working directory: $workdir"
echo "[INFO] Chromosome: $chr"

# Validate inputs
if [ ! -f "$bedfile" ]; then
    echo "[ERROR] BED file not found: $bedfile" >&2
    exit 1
fi

if [ ! -f "${plink_ref}.${chr}.bim" ]; then
    echo "[ERROR] Plink BIM file not found: ${plink_ref}.${chr}.bim" >&2
    exit 1
fi

mkdir -p "$workdir"

annot_file="${workdir}/baseline.${chr}.annot.gz"
echo "[INFO] Annotation file: $annot_file"

# Activate conda environment with ldsc
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate ldsc

echo "[INFO] Creating annotation file with make_annot.py"
if ! python3 "$HOME/modelowanie/ldsc-python3/make_annot.py" \
    --bed-file "$bedfile" \
    --bimfile "${plink_ref}.${chr}.bim" \
    --annot-file "$annot_file"; then
    echo "[ERROR] make_annot.py failed" >&2
    exit 1
fi

echo "Running make_annot.py:"
echo "  --bed-file $bedfile"
echo "  --bimfile ${plink_ref}.${chr}.bim"
echo "  --annot-file $annot_file"
ls -lh "$annot_file"

echo "[INFO] Calculating LD scores with ldsc.py"
if ! python3 "$HOME/modelowanie/ldsc-python3/ldsc.py" \
    --l2 \
    --bfile "${plink_ref}.${chr}" \
    --ld-wind-cm 1 \
    --annot "$annot_file" \
    --thin-annot \
    --out "${workdir}/baseline.${chr}" \
    --print-snps "$hm_snp"; then
    echo "[ERROR] ldsc.py failed" >&2
    exit 1
fi

echo "[INFO] Finished chromosome $chr"
