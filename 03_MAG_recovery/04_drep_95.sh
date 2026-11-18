#!/bin/bash
#SBATCH --job-name=drep_95
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=drep_95_%j.out
#SBATCH --error=drep_95_%j.err

# ============================================================
# dRep dereplication at 95% ANI (species-level)
# ============================================================
# Re-uses the CheckM2 genome_info.csv from the 99.9% run
# so CheckM2 does not need to be re-run.
#
# Quality thresholds same as before:
#   completeness ≥50%, contamination <10%
#
# Usage: sbatch 04b_drep_95.sh
# Run from MAGs directory
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── Config ────────────────────────────────────────────────────
BINS_DIR="04_drep/all_bins"
GENOME_INFO="04_drep/genome_info.csv"
OUTDIR="04_drep_95"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Verify inputs exist
if [[ ! -d "${BINS_DIR}" ]]; then
    echo "ERROR: bins directory not found: ${BINS_DIR}"
    echo "Expected the all_bins/ directory from the 99.9% dRep run"
    exit 1
fi

if [[ ! -f "${GENOME_INFO}" ]]; then
    echo "ERROR: genome_info.csv not found: ${GENOME_INFO}"
    exit 1
fi

N=$(ls "${BINS_DIR}"/*.fa | wc -l)
echo "[$(date)] Bins found: ${N}"
echo "[$(date)] Using existing CheckM2 quality table: ${GENOME_INFO}"
echo "[$(date)] Running dRep at 95% ANI (species-level)..."

conda activate drep_3.4.5

dRep dereplicate \
    "${OUTDIR}" \
    --genomes "${BINS_DIR}"/*.fa \
    --genomeInfo "${GENOME_INFO}" \
    --completeness 50 \
    --contamination 10 \
    --P_ani 0.9 \
    --S_ani 0.95 \
    --S_algorithm ANIn \
    --coverage_method larger \
    -p "${THREADS}"

N_MAGS=$(ls "${OUTDIR}/dereplicated_genomes/"*.fa 2>/dev/null | wc -l)
echo "[$(date)] dRep 95% complete."
echo "[$(date)] Dereplicated MAGs (species-level, 95% ANI): ${N_MAGS}"
echo "[$(date)] MAGs: ${OUTDIR}/dereplicated_genomes/"
