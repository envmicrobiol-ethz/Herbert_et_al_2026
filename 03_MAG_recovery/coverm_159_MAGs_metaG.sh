#!/bin/bash
#SBATCH --job-name=coverm_metaG_159
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=50G
#SBATCH --time=72:00:00
#SBATCH --output=coverm_metaG_159_%j.out
#SBATCH --error=coverm_metaG_159_%j.err

# ============================================================
# CoverM — map all metaG reads against all 159 MAGs
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

MAG_DIR="21_coassembly/MAGs/drep_95/dereplicated_genomes"
OUTDIR="21_coassembly/MAGs/coverm_metaG_159"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

READS=(
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode17_filtered.fastq.gz"
    "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
    "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
)

echo "[$(date)] Running CoverM: metaG reads vs 159 MAGs..."
echo "  MAGs: $(ls ${MAG_DIR}/*.fa | wc -l)"

coverm genome \
    --genome-fasta-files "${MAG_DIR}"/*.fa \
    --single "${READS[@]}" \
    --methods mean covered_bases relative_abundance \
    --min-covered-fraction 0 \
    --mapper minimap2-ont \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/all_159_MAGs_metaG_coverage.tsv"

echo "[$(date)] Done."
echo "Results: ${OUTDIR}/all_159_MAGs_metaG_coverage.tsv"
