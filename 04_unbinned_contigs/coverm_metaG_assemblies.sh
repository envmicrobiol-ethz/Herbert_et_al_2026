#!/bin/bash
#SBATCH --job-name=coverm_metaG_asm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=50G
#SBATCH --time=100:00:00
#SBATCH --output=coverm_metaG_asm_%j.out
#SBATCH --error=coverm_metaG_asm_%j.err

# ============================================================
# CoverM — map all metaG reads against combined assembly
# (renamed assemblies with prefixed contig names)
# Enables cross-sample expression/coverage analysis
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate coverm_0.6.1

RENAMED_ASM_DIR="15_transcriptm_assemblies/renamed_assemblies"
OUTDIR="coverm_metaG_assemblies"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"

# Concatenate renamed assemblies (prefixed contig names, no duplicates)
COMBINED="${OUTDIR}/all_assemblies_combined.fasta"
if [[ ! -f "${COMBINED}" ]]; then
    echo "[$(date)] Concatenating renamed assemblies..."
    cat "${RENAMED_ASM_DIR}"/*.fasta > "${COMBINED}"
    echo "  Total contigs: $(grep -c '^>' ${COMBINED})"
fi

# Sample reads
declare -A READS
READS["WIL1_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
READS["WIL2_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz"
READS["WIL2_after"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz"
READS["WILB_D3"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
READS["WILB_D10"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode17_filtered.fastq.gz"
READS["WILB_D0"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
READS["WILB_D1"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"

READS_ARGS=()
for SAMPLE in WIL1_before WIL2_before WIL2_after WILB_D0 WILB_D1 WILB_D3 WILB_D10; do
    R="${READS[${SAMPLE}]}"
    [[ -f "${R}" ]] && READS_ARGS+=("${R}") || echo "  [SKIP] ${SAMPLE} reads not found"
done

echo "[$(date)] Mapping $(echo ${#READS_ARGS[@]}) samples against combined assembly..."

coverm contig \
    --reference "${COMBINED}" \
    --single "${READS_ARGS[@]}" \
    --methods mean covered_bases rpkm \
    --min-covered-fraction 0 \
    --min-read-percent-identity 85 \
    --mapper minimap2-ont \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/all_assemblies_metaG_coverage.tsv"

echo "[$(date)] Done."
echo "  Contigs with coverage: $(tail -n +2 ${OUTDIR}/all_assemblies_metaG_coverage.tsv | wc -l)"
echo "  Results: ${OUTDIR}/all_assemblies_metaG_coverage.tsv"
