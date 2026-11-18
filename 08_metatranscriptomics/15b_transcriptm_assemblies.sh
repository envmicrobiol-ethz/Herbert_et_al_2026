#!/bin/bash
#SBATCH --job-name=transcriptm_assemblies
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=20G
#SBATCH --time=300:00:00
#SBATCH --output=transcriptm_assemblies_%j.out
#SBATCH --error=transcriptm_assemblies_%j.err

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1

# ── Config ────────────────────────────────────────────────────
READS_DIR="[metaT_reads_dir]"
ASSEMBLIES_DIR="[assemblies_dir]"
OUTDIR="15_transcriptm_assemblies"
DB="[transcriptm_db]/"

# Assembly contig names are generic (contig_1, contig_2 etc) across samples
# so we need to rename them to avoid duplicate sequence errors in coverm
RENAMED_DIR="15_transcriptm_assemblies/prokka_ready_assemblies"
mkdir -p "${RENAMED_DIR}"

echo "[$(date)] Renaming assembly contigs to avoid duplicates..."
for FASTA in "${ASSEMBLIES_DIR}"/*.fasta; do
    ASM=$(basename "${FASTA}" .fasta)
    awk -v asm="${ASM}" '/^>/{print ">"asm"__"substr($0,2); next}{print}' \
        "${FASTA}" > "${RENAMED_DIR}/${ASM}.fasta"
done
echo "[$(date)] Assemblies renamed: $(ls ${RENAMED_DIR}/*.fasta | wc -l)"

mkdir -p "${OUTDIR}"

# ── Run TranscriptM per sample ────────────────────────────────
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do

    PE1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    PE2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"

    if [[ ! -f "${PE1}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — R1 not found: ${PE1}"
        continue
    fi

    echo "[$(date)] Sample: ${SAMPLE}"

    ASM_OUT="${OUTDIR}/${SAMPLE}_assemblies"
    if [[ -d "${ASM_OUT}" ]]; then
        echo "[$(date)] [SKIP] ${SAMPLE} assemblies — output already exists"
    else
        echo "[$(date)] TranscriptM: ${SAMPLE} vs all assemblies..."
        transcriptm count \
            -1 "${PE1}" \
            -2 "${PE2}" \
            -db "${DB}" \
            -g "${RENAMED_DIR}" \
            -x fasta \
            -o "${ASM_OUT}" \
            --skip-qc
        echo "[$(date)] Done: ${SAMPLE} vs assemblies"
    fi

done

echo "[$(date)] TranscriptM assemblies complete. Outputs in: ${OUTDIR}/"
