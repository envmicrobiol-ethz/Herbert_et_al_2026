#!/bin/bash
#SBATCH --job-name=coasm_transcriptm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=50G
#SBATCH --time=300:00:00
#SBATCH --output=coasm_transcriptm_%j.out
#SBATCH --error=coasm_transcriptm_%j.err

# ============================================================
# TranscriptM for full co-assembly
# Covers all 211,503 contigs including 205,413 unbinned
# Following same approach as 15_transcriptm_assemblies
# WBe218_8 excluded to save time
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1

READS_DIR="[metaT_reads_dir]"
COASM_FA="21_coassembly/assembly/assembly.fasta"
OUTDIR="15_transcriptm_coasm"
DB="[transcriptm_db]/"

# Co-assembly contigs are named contig_XXXXX — rename with COASM__ prefix
# to avoid duplicate sequence name issues in coverm (same as assembly approach)
RENAMED_DIR="${OUTDIR}/prokka_ready_coasm"
mkdir -p "${RENAMED_DIR}"

if [[ ! -f "${RENAMED_DIR}/COASM.fasta" ]]; then
    echo "[$(date)] Renaming co-assembly contigs..."
    awk '/^>/{print ">COASM__"substr($0,2); next}{print}' \
        "${COASM_FA}" > "${RENAMED_DIR}/COASM.fasta"
    echo "[$(date)] Done: $(grep -c '^>' ${RENAMED_DIR}/COASM.fasta) contigs"
fi

mkdir -p "${OUTDIR}"

# ── Run TranscriptM per sample ────────────────────────────────
for SAMPLE in WAf15_8 WBe15_8; do

    PE1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"

    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — R1 not found: ${PE1}"
        continue
    fi

    COASM_OUT="${OUTDIR}/${SAMPLE}_coasm"

    if [[ -d "${COASM_OUT}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — output already exists"
        continue
    fi

    echo "[$(date)] TranscriptM: ${SAMPLE} vs co-assembly..."
    transcriptm count \
        -1 "${PE1}" \
        -2 "${PE2}" \
        -db "${DB}" \
        -g "${RENAMED_DIR}" \
        -x fasta \
        -o "${COASM_OUT}" \
        --skip-qc
    echo "[$(date)] Done: ${SAMPLE} co-assembly"

done

echo "[$(date)] Co-assembly TranscriptM complete. Outputs in: ${OUTDIR}/"
echo ""
echo "Next steps:"
echo "  1. Run GreeningLab DIAMOND on Prokka FAA from annotate/ directories"
echo "  2. Filter results to unbinned contigs using 21_coassembly/MAGs/coasm_contig_to_mag.tsv"
echo "  3. Join gene_id to TranscriptM final_table/rawcount_tpm_table.txt"
