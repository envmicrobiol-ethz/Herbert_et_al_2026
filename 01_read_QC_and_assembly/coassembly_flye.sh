#!/bin/bash
#SBATCH --job-name=coassembly_flye
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=50G
#SBATCH --time=120:00:00
#SBATCH --output=coassembly_flye_%j.out
#SBATCH --error=coassembly_flye_%j.err

# ============================================================
# Flye co-assembly of all Wildstrubel ONT reads
# ============================================================
# Concatenates all sample reads and runs a single Flye assembly
# to recover MAGs from organisms present across multiple samples
# at low individual-sample coverage.
#
# Usage: sbatch coassembly_flye.sh
# ============================================================

set -euo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate flye_V2.9

# ── Config ────────────────────────────────────────────────────
# Raw filtered reads — update paths if different
READS=(
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
    "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
    "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
    "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
)

OUTDIR="20_coassembly"
COMBINED_READS="${OUTDIR}/all_samples_combined.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUTDIR}"

# ── Concatenate all reads ─────────────────────────────────────
if [[ -s "${COMBINED_READS}" ]]; then
    echo "[$(date)] [SKIP] Combined reads already exist"
else
    echo "[$(date)] Concatenating reads from all samples..."
    cat "${READS[@]}" > "${COMBINED_READS}"
    N=$(zcat "${COMBINED_READS}" | awk 'NR%4==1' | wc -l)
    echo "[$(date)] Total reads: ${N}"
fi

# ── Run Flye co-assembly ──────────────────────────────────────
if [[ -f "${OUTDIR}/assembly/assembly.fasta" ]]; then
    echo "[$(date)] [SKIP] Assembly already exists"
else
    echo "[$(date)] Running Flye co-assembly..."
    flye \
        --nano-hq "${COMBINED_READS}" \
        --threads "${THREADS}" \
        --meta \
        --iterations 3 \
        --out-dir "${OUTDIR}/assembly" \
        > "${OUTDIR}/flye.log" 2>&1

    echo "[$(date)] Flye co-assembly complete."
fi

N_CONTIGS=$(grep -c "^>" "${OUTDIR}/assembly/assembly.fasta" 2>/dev/null || echo 0)
TOTAL_LEN=$(awk '/^>/{next}{total+=length($0)}END{print total}' \
    "${OUTDIR}/assembly/assembly.fasta" 2>/dev/null || echo 0)

echo "[$(date)] Assembly stats:"
echo "  Contigs: ${N_CONTIGS}"
echo "  Total length: ${TOTAL_LEN} bp"
echo "  Assembly: ${OUTDIR}/assembly/assembly.fasta"
echo ""
echo "Next steps:"
echo "  1. Run binning pipeline (scripts 02-04) on co-assembly"
echo "  2. Dereplicate with existing MAGs using dRep"
