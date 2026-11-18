#!/bin/bash
#SBATCH --job-name=mag_annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=20G
#SBATCH --time=100:00:00
#SBATCH --output=mag_annotation_%j.out
#SBATCH --error=mag_annotation_%j.err

# ============================================================
# Step 5: MAG taxonomy (GTDB-Tk) and coverages (CoverM)
# ============================================================
# GTDB-Tk (gtdbtk-2.6.1): GTDBTK_DATA_PATH is baked into the
#   conda environment via `conda env config vars set`, so no
#   manual export is needed here.
#
# CoverM (coverm_0.6.1): uses --mapper minimap2-ont, which
#   calls minimap2 internally — no separate mapping step needed.
#
# Note: archaeal MAGs (methanogens etc.) appear in the separate
#   gtdbtk.ar53.summary.tsv output file, not the bacterial one.
#
# Usage:
#   MAG_DIR=/path/to/04_drep/drep_output/dereplicated_genomes \
#   SAMPLE_MAP=sample_map.tsv \
#   sbatch 05_gtdbtk_coverm.sh
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── Config ────────────────────────────────────────────────────
MAG_DIR="${MAG_DIR:-04_drep/drep_output/dereplicated_genomes}"
SAMPLE_MAP="${SAMPLE_MAP:-sample_map.tsv}"
OUTDIR="${OUTDIR:-05_annotation}"
THREADS="${THREADS:-32}"

mkdir -p "${OUTDIR}/gtdbtk" "${OUTDIR}/coverm"


# ════════════════════════════════════════════════════════════
# GTDB-Tk
# ════════════════════════════════════════════════════════════
# GTDBTK_DATA_PATH is set via `conda env config vars` in the
# gtdbtk-2.6.1 environment — no manual export required.
conda activate gtdb-tk_2.4.0
export GTDBTK_DATA_PATH="[gtdbtk_db]"
echo "[$(date)] Running GTDB-Tk classify_wf..."

gtdbtk classify_wf \
    --genome_dir "${MAG_DIR}" \
    --extension fa \
    --out_dir "${OUTDIR}/gtdbtk" \
    --cpus "${THREADS}" \
    --skip_ani_screen

echo "[$(date)] GTDB-Tk complete."
echo "[$(date)] Bacterial summary: ${OUTDIR}/gtdbtk/gtdbtk.bac120.summary.tsv"
echo "[$(date)] Archaeal summary:  ${OUTDIR}/gtdbtk/gtdbtk.ar53.summary.tsv"
conda deactivate


# ════════════════════════════════════════════════════════════
# CoverM
# ════════════════════════════════════════════════════════════
conda activate coverm_0.6.1
echo "[$(date)] Running CoverM for MAG coverages..."

READS_ARGS=()
while IFS=$'\t' read -r SAMPLE READS; do
    [[ "${SAMPLE}" == \#* ]] && continue
    [[ -z "${SAMPLE}" ]]     && continue
    READS_ARGS+=(--single "${READS}")
done < "${SAMPLE_MAP}"

coverm genome \
    --genome-fasta-directory "${MAG_DIR}" \
    --genome-fasta-extension fa \
    "${READS_ARGS[@]}" \
    --methods mean covered_bases relative_abundance \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-ont \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/coverm/MAG_coverages.tsv"

echo "[$(date)] CoverM complete."
echo "[$(date)] Coverage table: ${OUTDIR}/coverm/MAG_coverages.tsv"
conda deactivate


# ── Quick summary ─────────────────────────────────────────────
echo ""
echo "=== Pipeline complete ==="
echo "MAGs:       ${MAG_DIR}"
echo "Taxonomy:   ${OUTDIR}/gtdbtk/"
echo "Coverages:  ${OUTDIR}/coverm/MAG_coverages.tsv"
echo ""
N_BAC=$(awk 'END{print NR-1}' "${OUTDIR}/gtdbtk/gtdbtk.bac120.summary.tsv" 2>/dev/null || echo "?")
N_ARC=$(awk 'END{print NR-1}' "${OUTDIR}/gtdbtk/gtdbtk.ar53.summary.tsv"  2>/dev/null || echo "?")
echo "Bacterial MAGs: ${N_BAC}"
echo "Archaeal MAGs:  ${N_ARC}"
