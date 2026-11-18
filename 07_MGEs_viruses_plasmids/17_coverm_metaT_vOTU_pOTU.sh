#!/bin/bash
#SBATCH --job-name=coverm_metaT_vp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=coverm_metaT_vp_%j.out
#SBATCH --error=coverm_metaT_vp_%j.err

# ============================================================
# CoverM — map metaT reads to vOTUs and pOTUs
# ============================================================
# Uses minimap2-sr preset for short-read Illumina metaT data.
# Reports mean coverage, covered bases, relative abundance
# and reads_per_base (length-normalised, useful for expression).
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

READS_DIR="[metaT_reads_dir]"
VOTUS_FASTA="07_viruses/vOTUs/vOTUs.fasta"
POTUS_DIR="08_plasmids/pOTUs/dereplicated_genomes"
OUTDIR="17_coverm_metaT_vp"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

# ── Build read arguments ──────────────────────────────────────
READS_ARGS=()
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    PE1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — trimmed R1 not found"
        continue
    fi
    READS_ARGS+=(--coupled "${PE1}" "${PE2}")
done

# ── vOTUs ─────────────────────────────────────────────────────
echo "[$(date)] Running CoverM: metaT reads vs vOTUs..."

coverm genome \
    --genome-fasta-files "${VOTUS_FASTA}" \
    "${READS_ARGS[@]}" \
    --methods mean covered_bases relative_abundance reads_per_base \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-sr \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/metaT_vOTU_coverages.tsv"

echo "[$(date)] vOTU coverage done: ${OUTDIR}/metaT_vOTU_coverages.tsv"

# ── pOTUs — concatenate first then run coverm ─────────────────
echo "[$(date)] Concatenating pOTU FASTAs..."
cat "${POTUS_DIR}"/*.fa > "${OUTDIR}/pOTUs_combined.fasta"

echo "[$(date)] Running CoverM: metaT reads vs pOTUs..."

coverm genome \
    --genome-fasta-files "${OUTDIR}/pOTUs_combined.fasta" \
    "${READS_ARGS[@]}" \
    --methods mean covered_bases relative_abundance reads_per_base \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-sr \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/metaT_pOTU_coverages.tsv"

echo "[$(date)] pOTU coverage done: ${OUTDIR}/metaT_pOTU_coverages.tsv"
echo "[$(date)] All done."
