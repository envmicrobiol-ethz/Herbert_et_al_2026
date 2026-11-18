#!/bin/bash
#SBATCH --job-name=coverm_metaT_coasm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=72:00:00
#SBATCH --output=coverm_metaT_coasm_%j.out
#SBATCH --error=coverm_metaT_coasm_%j.err

# ============================================================
# CoverM — map metaT reads to novel co-assembly MAGs
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

NOVEL_MAGS="21_coassembly/MAGs/drep_95/dereplicated_genomes"
READS_DIR="[metaT_reads_dir]"
OUTDIR="21_coassembly/MAGs/coverm_metaT"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

echo "[$(date)] Running CoverM: metaT reads vs co-assembly MAGs..."
echo "  MAGs: $(ls ${NOVEL_MAGS}/COASM_*.fa | wc -l)"

coverm genome \
    --genome-fasta-files "${NOVEL_MAGS}"/COASM_*.fa \
    --coupled \
        "${READS_DIR}/WAf15_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WAf15_8_R2.trimmed.fastq.gz" \
    --coupled \
        "${READS_DIR}/WBe15_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WBe15_8_R2.trimmed.fastq.gz" \
    --coupled \
        "${READS_DIR}/WBe218_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WBe218_8_R2.trimmed.fastq.gz" \
    --methods mean covered_bases relative_abundance reads_per_base \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-sr \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/COASM_MAGs_metaT_coverage.tsv"

echo "[$(date)] Done."
echo "Results: ${OUTDIR}/COASM_MAGs_metaT_coverage.tsv"
