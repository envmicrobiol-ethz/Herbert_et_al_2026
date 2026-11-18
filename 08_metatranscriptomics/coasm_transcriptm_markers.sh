#!/bin/bash
#SBATCH --job-name=coasm_tm_markers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=50G
#SBATCH --time=300:00:00
#SBATCH --output=coasm_tm_markers_%j.out
#SBATCH --error=coasm_tm_markers_%j.err

# ============================================================
# TranscriptM for co-assembly marker gene contigs
# Uses pre-built Prodigal ref + GFF (no internal Prokka)
# Only maps to contigs containing GreeningLab/HydDB PASS genes
#
# Prerequisites (must exist):
#   21_coassembly/transcriptm_markers/ref/COASM_markers_ref.fna
#   21_coassembly/transcriptm_markers/ref/COASM_markers_ref.gff
#
# Samples: WAf15_8, WBe15_8 (WBe218_8 excluded to save time)
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1

READS_DIR="[metaT_reads_dir]"
REF="21_coassembly/transcriptm_markers/ref/COASM_markers_ref.fna"
GFF="21_coassembly/transcriptm_markers/ref/COASM_markers_ref.gff"
OUTDIR="21_coassembly/transcriptm_markers"
DB="[transcriptm_db]/"

# Verify prerequisites
if [[ ! -f "${REF}" ]]; then
    echo "ERROR: ${REF} not found — run coasm_prodigal_diamond.sh first"
    exit 1
fi
if [[ ! -f "${GFF}" ]]; then
    echo "ERROR: ${GFF} not found — run coasm_prodigal_diamond.sh first"
    exit 1
fi

echo "[$(date)] Reference: ${REF} ($(grep -c '^>' ${REF}) contigs)"
echo "[$(date)] GFF: ${GFF} ($(grep -c $'\tCDS\t' ${GFF}) CDS)"

for SAMPLE in WAf15_8 WBe15_8; do
    PE1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    SAMPLE_OUT="${OUTDIR}/${SAMPLE}_coasm_markers"

    if [[ -d "${SAMPLE_OUT}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — output already exists"
        continue
    fi

    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — R1 not found: ${PE1}"
        continue
    fi

    echo "[$(date)] TranscriptM: ${SAMPLE}..."
    transcriptm count \
        -1 "${PE1}" \
        -2 "${PE2}" \
        -n 24 \
        -m 480 \
        --ref "${REF}" \
        --gff "${GFF}" \
        -db "${DB}" \
        -o "${SAMPLE_OUT}" \
        --skip-qc

    echo "[$(date)] Done: ${SAMPLE}"
    echo "[$(date)] Genes with coverage: $(tail -n +2 ${SAMPLE_OUT}/final_table/rawcount_tpm_table.txt 2>/dev/null | wc -l)"
done

echo ""
echo "[$(date)] All done."
echo "Results:"
for SAMPLE in WAf15_8 WBe15_8; do
    TPM="${OUTDIR}/${SAMPLE}_coasm_markers/final_table/rawcount_tpm_table.txt"
    [[ -f "${TPM}" ]] && echo "  ${SAMPLE}: $(tail -n +2 ${TPM} | wc -l) genes"
done
echo ""
echo "Join key: gene_id in rawcount_tpm_table.txt"
echo "  → matches Prodigal protein IDs in COASM_markers_ref.faa"
echo "  → link to GreeningLab/HydDB hits via contig_id column"
