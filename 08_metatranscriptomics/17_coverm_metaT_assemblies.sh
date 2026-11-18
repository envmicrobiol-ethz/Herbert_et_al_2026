#!/bin/bash
#SBATCH --job-name=coverm_metaT_asm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=300:00:00
#SBATCH --output=coverm_metaT_asm_%j.out
#SBATCH --error=coverm_metaT_asm_%j.err

# ============================================================
# CoverM — map metaT reads to all Wildstrubel assemblies
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

READS_DIR="[metaT_reads_dir]"
ASSEMBLIES_DIR="../Wildstrubel_assemblies"
OUTDIR="17_coverm_metaT_assemblies"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

# ── Build read arguments ──────────────────────────────────────
READS_ARGS=()
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    PE1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — R1 not found"
        continue
    fi
    echo "[$(date)] Adding sample: ${SAMPLE}"
    READS_ARGS+=(--coupled "${PE1}" "${PE2}")
done

# ── Run CoverM per assembly ───────────────────────────────────
for ASM in WIL1_before WIL2_before WIL2_after WILB_D1 WILB_D3 WILB_D10 WILD_D0; do
    FASTA="${ASSEMBLIES_DIR}/${ASM}.fasta"
    OUT="${OUTDIR}/${ASM}_metaT_coverage.tsv"

    if [[ ! -f "${FASTA}" ]]; then
        echo "[$(date)] [SKIP] ${ASM} — assembly not found"
        continue
    fi
    if [[ -f "${OUT}" ]]; then
        echo "[$(date)] [SKIP] ${ASM} — already done"
        continue
    fi

    echo "[$(date)] CoverM: metaT reads vs ${ASM}..."
    coverm contig \
        --reference "${FASTA}" \
        "${READS_ARGS[@]}" \
        --methods mean covered_bases rpkm \
        --min-covered-fraction 0 \
        --min-read-percent-identity 90 \
        --mapper minimap2-sr \
        --threads "${THREADS}" \
        --output-file "${OUT}"

    echo "[$(date)] Done: ${ASM} — $(tail -n +2 ${OUT} | wc -l) contigs"
done

echo "[$(date)] All done. Results: ${OUTDIR}/"
