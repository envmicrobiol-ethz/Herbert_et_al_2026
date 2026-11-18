#!/bin/bash
#SBATCH --job-name=read_mapping
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=read_mapping_%j.out
#SBATCH --error=read_mapping_%j.err

# ============================================================
# Step 1: Map reads back to assembly using minimap2
# ============================================================
# minimap2 and samtools are in separate environments/modules:
#   minimap2  → conda activate minimap2_V2.22
#   samtools  → module load (gcc required first)
# Both are active simultaneously in the same shell session.
#
# sample_map.tsv format (tab-separated, no header):
#   sample_name <TAB> /path/to/reads.fastq.gz
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate minimap2_V2.22

module load stack/.2024-06-silent gcc/12.2.0
module load samtools/1.16.1

# ── Config ────────────────────────────────────────────────────
ASSEMBLY="${ASSEMBLY:-../Wildstrubel_assemblies/WILB_D1.fasta}"
SAMPLE_MAP="${SAMPLE_MAP:-sample_map.tsv}"
OUTDIR="${OUTDIR:-01_mapping}"
THREADS="${THREADS:-16}"

# ── Setup ─────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"

ASSEMBLY_NAME=$(basename "${ASSEMBLY}" .fasta)
echo "[$(date)] Starting mapping to assembly: ${ASSEMBLY_NAME}"

echo "[$(date)] Indexing assembly..."
minimap2 -d "${OUTDIR}/${ASSEMBLY_NAME}.mmi" "${ASSEMBLY}" -t "${THREADS}"

# ── Map each sample ───────────────────────────────────────────
while IFS=$'\t' read -r SAMPLE READS; do

    [[ "${SAMPLE}" == \#* ]] && continue
    [[ -z "${SAMPLE}" ]]     && continue

    echo "[$(date)] Mapping sample: ${SAMPLE}"

    BAM="${OUTDIR}/${SAMPLE}_vs_${ASSEMBLY_NAME}.bam"

    minimap2 \
        -ax map-ont \
        -t "${THREADS}" \
        --secondary=no \
        "${OUTDIR}/${ASSEMBLY_NAME}.mmi" \
        "${READS}" \
    | samtools view -bS -F 4 - \
    | samtools sort -@ "${THREADS}" -o "${BAM}" -

    samtools index "${BAM}"

    echo "[$(date)] Done: ${BAM}"

done < "${SAMPLE_MAP}"

echo "[$(date)] All mappings complete."
