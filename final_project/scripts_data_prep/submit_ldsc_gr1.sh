#!/bin/bash
#SBATCH --job-name=ldsc_group1
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=8:00:00
#SBATCH --output=../logs/ldsc_group1_%j.out
#SBATCH --error=../logs/ldsc_group1_%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

for bedfile in \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v154.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v155.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v156.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v157.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v158.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v159.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v160.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v161.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v162.bed \
    /home/zmilczarska/modelowanie/TDEP-sLDSC/data/custom/amygdala/v163.bed
do
    celltype_name=$(basename "$bedfile" .bed)
    mkdir -p "/home/zmilczarska/modelowanie/TDEP-sLDSC/results/ldsc_${celltype_name}"
    echo "Processing $celltype_name"
    bash /home/zmilczarska/modelowanie/TDEP-sLDSC/scripts/cell-type-ldscores.sh "$bedfile"
done
