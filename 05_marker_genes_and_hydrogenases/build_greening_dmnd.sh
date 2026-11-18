#!/bin/bash
#SBATCH --job-name=build_greening_dmnd
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --output=build_greening_dmnd_%j.out
#SBATCH --error=build_greening_dmnd_%j.err

# ============================================================
# Build DIAMOND databases from all GreeningLab FASTAs
# 2020 database (51 FASTAs) + 2024 extended (8 FASTAs)
# ============================================================

set -eo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate diamond_2.1.10

DB_2020="HMMs/GreeningLab-database-main/Original database (2020)/Greening lab metabolic marker gene databases"
DB_2024="HMMs/GreeningLab-database-main/Extended_database_2024"
DMND_DIR="HMMs/greening_dmnd"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${DMND_DIR}"

# ── 2020 database FASTAs ──────────────────────────────────────
echo "[$(date)] Building DIAMOND databases from 2020 database..."
for FA in "${DB_2020}"/*.fasta; do
    NAME=$(basename "${FA}" .fasta | tr ' ' '_')
    DMND="${DMND_DIR}/${NAME}.dmnd"

    if [[ -f "${DMND}" ]]; then
        echo "  [SKIP] ${NAME}"
        continue
    fi

    N=$(grep -c "^>" "${FA}" 2>/dev/null || echo 0)
    echo "  [$(date '+%H:%M:%S')] Building: ${NAME} (${N} seqs)..."
    diamond makedb \
        --in "${FA}" \
        -d "${DMND_DIR}/${NAME}" \
        --threads "${THREADS}" \
        --quiet
    echo "  Done: $(ls -lh ${DMND} | awk '{print $5}')"
done

# ── 2024 extended database FASTAs ─────────────────────────────
echo ""
echo "[$(date)] Building DIAMOND databases from 2024 extended database..."
for FA in "${DB_2024}"/*.faa; do
    NAME=$(basename "${FA}" .faa | tr ' ' '_')
    DMND="${DMND_DIR}/${NAME}.dmnd"

    if [[ -f "${DMND}" ]]; then
        echo "  [SKIP] ${NAME}"
        continue
    fi

    N=$(grep -c "^>" "${FA}" 2>/dev/null || echo 0)
    echo "  [$(date '+%H:%M:%S')] Building: ${NAME} (${N} seqs)..."
    diamond makedb \
        --in "${FA}" \
        -d "${DMND_DIR}/${NAME}" \
        --threads "${THREADS}" \
        --quiet
    echo "  Done: $(ls -lh ${DMND} | awk '{print $5}')"
done

echo ""
echo "[$(date)] Done."
echo "Total DIAMOND databases: $(ls ${DMND_DIR}/*.dmnd | wc -l)"
ls -lh "${DMND_DIR}/"
