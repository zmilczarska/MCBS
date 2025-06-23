#!/bin/bash
#SBATCH --job-name=ldsc_chr1-7
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=14GB  # 2GB per chromosome
#SBATCH --time=08:00:00
#SBATCH --output=results/ldsc_chr1-7_hipp.out
#SBATCH --error=results/ldsc_chr1-7_hipp.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PROJECT_DIR="$HOME/modelowanie/TDEP-sLDSC"
cd "$PROJECT_DIR" || exit 1

OUTPUT_DIR="$PROJECT_DIR/results/ldsc_Hipp"
mkdir -p "$OUTPUT_DIR"

SUMSTATS="$PROJECT_DIR/data/scz2022_ldsc.sumstats.gz"
BASELINE_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_EUR_Phase3_baseline"
CELLTYPE_DIR="$PROJECT_DIR/data/Hippocampal_CA1_3"
WEIGHTS_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_Phase3_weights_hm3_no_MHC"
FREQ_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_Phase3_frq"

for CHR in {1..7}; do
    echo "Processing chromosome $CHR"
    OUTPUT_PREFIX="$OUTPUT_DIR/scz2022_chr${CHR}"
    
    python ~/modelowanie/ldsc-python3/ldsc.py \
        --h2 "$SUMSTATS" \
        --ref-ld-chr "$BASELINE_DIR/baseline.$CHR,$CELLTYPE_DIR/baseline.$CHR" \
        --w-ld-chr "$WEIGHTS_DIR/weights.hm3_noMHC.$CHR" \
        --frqfile-chr "$FREQ_DIR/1000G.EUR.QC.$CHR" \
        --overlap-annot \
        --print-coefficients \
        --out "$OUTPUT_PREFIX"
    
    if [ $? -eq 0 ]; then
        echo "Successfully processed chromosome $CHR"
    else
        echo "Failed to process chromosome $CHR"
    fi
done

echo "Finished processing chromosomes 1-7"
