#!/bin/bash
#SBATCH --job-name=prodigal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --output=prodigal_%j.out
#SBATCH --error=prodigal_%j.err

# ============================================================
# Step 9: Predict proteins from vOTU and pOTU sequences
# ============================================================
# Prodigal -p meta is used (metagenome mode) since these are
# environmental viral/plasmid contigs, not isolate genomes.
# Output .faa files are passed directly to the DIAMOND/LCA
# taxonomy script (10_diamond_lca_taxonomy.sh).
#
# TODO: confirm prodigal environment name on your cluster.
# Check with: find .../envs -name "prodigal" -maxdepth 3
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# TODO: update prodigal env name
PRODIGAL_ENV="prodigal_env"   # <── update this

VOTUS_FASTA="${VOTUS_FASTA:-07_viruses/vOTUs/vOTUs.fasta}"
POTUS_FASTA="${POTUS_FASTA:-08_plasmids/pOTUs/dereplicated_genomes}"
OUTDIR="${OUTDIR:-09_proteins}"

mkdir -p "${OUTDIR}"

module load stack/.2024-05-silent  gcc/13.2.0 prodigal/2.6.3

# ── vOTUs ─────────────────────────────────────────────────────
echo "[$(date)] Predicting proteins for vOTUs..."
prodigal \
    -i "${VOTUS_FASTA}" \
    -a "${OUTDIR}/vOTUs.faa" \
    -p meta \
    -f gff \
    -o "${OUTDIR}/vOTUs.gff" \
    -q   # quiet mode

N=$(grep -c "^>" "${OUTDIR}/vOTUs.faa")
echo "[$(date)] vOTU proteins predicted: ${N}"

# ── pOTUs ─────────────────────────────────────────────────────
# pOTUs are split across individual FASTAs in the dRep output —
# concatenate them first then run Prodigal on the merged set
echo "[$(date)] Predicting proteins for pOTUs..."

cat "${POTUS_FASTA}"/*.fa > "${OUTDIR}/pOTUs_merged.fasta"

prodigal \
    -i "${OUTDIR}/pOTUs_merged.fasta" \
    -a "${OUTDIR}/pOTUs.faa" \
    -p meta \
    -f gff \
    -o "${OUTDIR}/pOTUs.gff" \
    -q

N=$(grep -c "^>" "${OUTDIR}/pOTUs.faa")
echo "[$(date)] pOTU proteins predicted: ${N}"

conda deactivate

echo "[$(date)] Protein prediction complete."
echo "[$(date)] Output FAAs:"
echo "    ${OUTDIR}/vOTUs.faa"
echo "    ${OUTDIR}/pOTUs.faa"
echo ""
echo "Submit taxonomy with:"
echo "  bash 10_diamond_lca_taxonomy.sh ${OUTDIR}/vOTUs.faa ${OUTDIR}/pOTUs.faa"
