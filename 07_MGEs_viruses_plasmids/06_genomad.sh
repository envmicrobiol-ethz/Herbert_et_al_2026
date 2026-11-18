#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8192
#SBATCH --time=24:00:00
#SBATCH --job-name=genomad
#SBATCH --output=genomad_%j.out
#SBATCH --error=genomad_%j.err

# ============================================================
# Step 6: geNomad — identify viral and plasmid contigs
# ============================================================
# Runs geNomad end-to-end on all 8 Wildstrubel assemblies.
# Both viral and plasmid contigs are identified in one run.
# Output is used by:
#   07_checkv_votu.sh  — viral contigs → CheckV → vOTUs
#   08_plasmid_potu.sh — plasmid contigs → dRep → pOTUs
#
# Note: --disable-find-proviruses used to match group convention.
# Remove this flag if you want integrated proviruses identified.
# ============================================================

set -euo pipefail

# ── Environment ───────────────────────────────────────────────
module load stack/2024-06 gcc/12.2.0 python/3.11.6 mmseqs2/14-7e284 aragorn/1.2.41
source [genomad_venv]/bin/activate

# ── Config ────────────────────────────────────────────────────
THREADS=20
ASSEMBLIES_DIR="../Wildstrubel_assemblies"
OUTPUT_BASE="${OUTPUT_BASE:-06_genomad}"
DB_PATH="[genomad_db]"

mkdir -p "${OUTPUT_BASE}"

# ── Run geNomad on each assembly ──────────────────────────────
for FASTA in "${ASSEMBLIES_DIR}"/*.fasta; do
    SAMPLE=$(basename "${FASTA}" .fasta)
    SAMPLE_OUTDIR="${OUTPUT_BASE}/${SAMPLE}"

    echo "[$(date)] Running geNomad for: ${SAMPLE}"

    genomad end-to-end \
        --cleanup \
        "${FASTA}" \
        "${SAMPLE_OUTDIR}" \
        "${DB_PATH}" \
        --splits 8 \
        --enable-score-calibration \
        --disable-find-proviruses \
        --threads "${THREADS}"

    echo "[$(date)] Finished geNomad for: ${SAMPLE}"
done

echo "[$(date)] All assemblies complete."
echo "[$(date)] Viral summaries:   ${OUTPUT_BASE}/<sample>/<sample>_summary/<sample>_virus_summary.tsv"
echo "[$(date)] Plasmid summaries: ${OUTPUT_BASE}/<sample>/<sample>_summary/<sample>_plasmid_summary.tsv"
