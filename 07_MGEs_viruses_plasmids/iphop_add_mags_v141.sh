#!/bin/bash
#SBATCH -A [slurm_account]
#SBATCH -n 1
#SBATCH --cpus-per-task=40
#SBATCH --time=12-48:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=iphop_add_mags
#SBATCH --output=iphop_add_mags_%j.out
#SBATCH --error=iphop_add_mags_%j.err

set -eo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate iphop_env_1.4.1

MAG_DIR="21_coassembly/MAGs/drep_95/dereplicated_genomes"
GTDB_DIR="23_phylogenetics/gtdbtk_denovo_bac"
BASE_DB="[iphop_db]"
NEW_DB="iphop_results/Jun_2025_pub_rw_w_Wildstrubel_MAGs"
VOTUS="07_viruses/vOTUs/vOTUs.fasta"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

mkdir -p 26_iphop_results

# ════════════════════════════════════════════════════════════
# 1. Add 159 MAGs to iPHoP database
# ════════════════════════════════════════════════════════════
if [[ -d "${NEW_DB}" ]]; then
    echo "[$(date)] [SKIP] Custom database already exists"
else
    echo "[$(date)] Adding 159 MAGs to iPHoP database..."
    echo "  MAGs: $(ls ${MAG_DIR}/*.fa | wc -l)"
    echo "  GTDB-Tk dir: ${GTDB_DIR}"
    echo "  Base DB: ${BASE_DB}"

    iphop add_to_db \
        --fna_dir "${MAG_DIR}" \
        --gtdb_dir "${GTDB_DIR}" \
        --db_dir "${BASE_DB}" \
        --out_dir "${NEW_DB}"

    echo "[$(date)] Custom database built: ${NEW_DB}"
fi

# ════════════════════════════════════════════════════════════
# 2. Re-run host prediction for vOTUs
# ════════════════════════════════════════════════════════════
if [[ -f "iphop_results/vOTUs_with_MAGs/Host_prediction_to_genus_m90.csv" ]]; then
    echo "[$(date)] [SKIP] vOTU host prediction already done"
else
    echo "[$(date)] Running iPHoP predict on vOTUs..."
    iphop predict \
        --fa_file "${VOTUS}" \
        --db_dir "${NEW_DB}" \
        --out_dir "26_iphop_results/vOTUs_with_MAGs" \
        --num_threads "${THREADS}"
    echo "[$(date)] vOTU host prediction complete."
fi

# ════════════════════════════════════════════════════════════
# 3. Re-run host prediction for pOTUs
# ════════════════════════════════════════════════════════════
POTUS_FASTA="iphop_results/pOTUs_combined.fasta"
if [[ ! -f "${POTUS_FASTA}" ]]; then
    cat 08_plasmids/pOTUs/dereplicated_genomes/*.fa > "${POTUS_FASTA}"
fi

if [[ -f "iphop_results/pOTUs_with_MAGs/Host_prediction_to_genus_m90.csv" ]]; then
    echo "[$(date)] [SKIP] pOTU host prediction already done"
else
    echo "[$(date)] Running iPHoP predict on pOTUs..."
    iphop predict \
        --fa_file "${POTUS_FASTA}" \
        --db_dir "${NEW_DB}" \
        --out_dir "iphop_results/pOTUs_with_MAGs" \
        --num_threads "${THREADS}"
    echo "[$(date)] pOTU host prediction complete."
fi

echo ""
echo "[$(date)] Done."
echo "Custom DB:        ${NEW_DB}"
echo "vOTU predictions: iphop_results/vOTUs_with_MAGs/Host_prediction_to_genus_m90.csv"
echo "pOTU predictions: iphop_results/pOTUs_with_MAGs/Host_prediction_to_genus_m90.csv"
