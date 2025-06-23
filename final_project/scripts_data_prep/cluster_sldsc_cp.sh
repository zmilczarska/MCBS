#!/bin/bash
#SBATCH --job-name=tdep_sldsc_full
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --time=8:00:00
#SBATCH --output=results_sldsc/tdep_sldsc_full.out
#SBATCH --error=results_sldsc/tdep_sldsc_full.err


source ~/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PROJECT_DIR="$HOME/modelowanie/TDEP-sLDSC"
cd "$PROJECT_DIR" || exit 1

OUTPUT_DIR="$PROJECT_DIR/results_sldsc"
mkdir -p "$OUTPUT_DIR"

SUMSTATS="$PROJECT_DIR/data/scz2022_ldsc.sumstats.gz"
BASELINE_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_EUR_Phase3_baseline"
WEIGHTS_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_Phase3_weights_hm3_no_MHC"
FREQ_DIR="$PROJECT_DIR/data/sldsc_ref/1000G_Phase3_frq"

# Find all cell type directories
CELLTYPE_DIRS=($(ls -d $PROJECT_DIR/{data,results}/ldsc_v* 2>/dev/null))

if [ ${#CELLTYPE_DIRS[@]} -eq 0 ]; then
    echo "ERROR: No cell type directories found!" >&2
    exit 1
fi

# Function to count annotation columns
count_annot_cols() {
    local file=$1
    if [ -f "$file" ]; then
        zcat -f "$file" | head -1 | tr '\t' '\n' | wc -l
    else
        echo 0
    fi
}

for CELLTYPE_DIR in "${CELLTYPE_DIRS[@]}"; do
    CELLTYPE_NAME=$(basename "$CELLTYPE_DIR")
    CELLTYPE_OUTPUT="$OUTPUT_DIR/$CELLTYPE_NAME"
    mkdir -p "$CELLTYPE_OUTPUT"
    
    echo "Processing: $CELLTYPE_NAME"
    
    for CHR in {1..22}; do
        (
            OUTPUT_PREFIX="$CELLTYPE_OUTPUT/scz2022_chr${CHR}"
            BASELINE_ANN="$BASELINE_DIR/baseline.$CHR.annot.gz"
            CELL_ANN="$CELLTYPE_DIR/baseline.$CHR.annot.gz"
            
            # Check if files exist
            if [ ! -f "$BASELINE_ANN" ]; then
                echo "Missing baseline file: $BASELINE_ANN" > "${OUTPUT_PREFIX}.err"
                continue
            fi
            
            if [ ! -f "$CELL_ANN" ]; then
                echo "Missing cell type file: $CELL_ANN" > "${OUTPUT_PREFIX}.err"
                continue
            fi
            
            # Count annotation columns
            BASE_COLS=$(count_annot_cols "$BASELINE_ANN")
            CELL_COLS=$(count_annot_cols "$CELL_ANN")
            
            if [ "$BASE_COLS" -ne "$CELL_COLS" ]; then
                echo "Annotation mismatch (baseline: $BASE_COLS vs cell: $CELL_COLS columns)" > "${OUTPUT_PREFIX}.err"
                continue
            fi
            
            # Run LDSC
            python ~/modelowanie/ldsc-python3/ldsc.py \
                --h2 "$SUMSTATS" \
                --ref-ld-chr "$BASELINE_DIR/baseline.$CHR,$CELLTYPE_DIR/baseline.$CHR" \
                --w-ld-chr "$WEIGHTS_DIR/weights.hm3_noMHC.$CHR" \
                --frqfile-chr "$FREQ_DIR/1000G.EUR.QC.$CHR" \
                --overlap-annot \
                --print-coefficients \
                --out "$OUTPUT_PREFIX" 2> "${OUTPUT_PREFIX}.tmp_err"
            
            # Check exit status
            if [ $? -eq 0 ]; then
                echo "Successfully processed chr$CHR for $CELLTYPE_NAME"
                rm -f "${OUTPUT_PREFIX}.err" "${OUTPUT_PREFIX}.tmp_err"
            else
                echo "Failed to process chr$CHR for $CELLTYPE_NAME"
                mv "${OUTPUT_PREFIX}.tmp_err" "${OUTPUT_PREFIX}.err"
                echo "Exit status: $?" >> "${OUTPUT_PREFIX}.err"
            fi
        ) &
    done
    wait
    
    echo "Finished processing $CELLTYPE_NAME"
done

echo "All cell types processed. Results saved in $OUTPUT_DIR"
