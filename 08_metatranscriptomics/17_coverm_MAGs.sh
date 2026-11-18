#!/bin/bash
#SBATCH --job-name=coverm_metaT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=50G
#SBATCH --time=144:00:00
#SBATCH --output=coverm_metaT_%j.out
#SBATCH --error=coverm_metaT_%j.err

# ============================================================
# CoverM — map metaT reads to MAGs
# ============================================================
# Maps paired-end metaT reads from all three samples against
# the dereplicated MAG set and reports coverage metrics.
#
# Uses minimap2-sr preset (short reads, paired-end Illumina)
# rather than minimap2-ont since metaT is short-read.
#
# Usage: sbatch coverm_metaT_MAGs.sh
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

# ── Config ────────────────────────────────────────────────────
READS_DIR="[metaT_reads_dir]"
MAG_DIR="04_drep/drep_output/dereplicated_genomes"
OUTDIR="17_coverm_metaT"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

# ── Build read arguments ──────────────────────────────────────
# CoverM takes paired reads via --coupled
READS_ARGS=()
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    PE1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"
    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — R1 not found: ${PE1}"
        continue
    fi
    echo "[$(date)] Adding sample: ${SAMPLE}"
    READS_ARGS+=(--coupled "${PE1}" "${PE2}")
done

# ── Run CoverM ────────────────────────────────────────────────
echo "[$(date)] Running CoverM (metaT reads vs MAGs)..."

coverm genome \
    --genome-fasta-directory "${MAG_DIR}" \
    --genome-fasta-extension fa \
    "${READS_ARGS[@]}" \
    --methods mean covered_bases relative_abundance reads_per_base \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-sr \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/metaT_MAG_coverages.tsv"

echo "[$(date)] CoverM complete."
echo "[$(date)] Coverage table: ${OUTDIR}/metaT_MAG_coverages.tsv"
