#!/bin/bash
#SBATCH --job-name=transcriptm
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=144:00:00
#SBATCH --output=transcriptm_%j.out
#SBATCH --error=transcriptm_%j.err

source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1

READS_DIR="[metaT_reads_dir]"
MAG_DIR="prokka_ready_mags"
OUTDIR="15_transcriptm"

for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    R1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    OUT="${OUTDIR}/${SAMPLE}_MAGs"
    echo "Running TranscriptM for sample: ${SAMPLE} vs MAGs"
    transcriptm count \
        -1 "$R1" \
        -2 "$R2" \
        -g "$MAG_DIR" \
        -x fasta \
        -o "$OUT" \
        --skip-qc
done

echo "=== ALL TranscriptM runs completed ==="
