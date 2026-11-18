#!/bin/bash
#SBATCH --job-name=plasmid_potu
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=6:00:00
#SBATCH --output=plasmid_potu_%j.out
#SBATCH --error=plasmid_potu_%j.err

# ============================================================
# Step 8: Plasmid dereplication → pOTUs
# ============================================================
# 1. Merge plasmid contigs (≥1 kb) from all geNomad outputs
# 2. Dereplicate with dRep at 95% ANI, 50% coverage → pOTUs
#
# dRep is used rather than CD-HIT-EST for plasmids because
# plasmid sequences are typically longer and more structurally
# complex — dRep's two-stage clustering handles this better.
#
# Usage:
#   GENOMAD_DIR=/path/to/06_genomad \
#   sbatch 08_plasmid_potu.sh
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── Config ────────────────────────────────────────────────────
GENOMAD_DIR="${GENOMAD_DIR:-06_genomad}"
OUTDIR="${OUTDIR:-08_plasmids}"
THREADS="${THREADS:-16}"

MIN_CONTIG_LEN=1000   # ≥1 kb as per methods
DREP_ANI=0.95         # 95% ANI for pOTUs
DREP_COV=0.50         # 50% coverage for pOTUs

mkdir -p "${OUTDIR}/merged" \
         "${OUTDIR}/split_fastas" \
         "${OUTDIR}/pOTUs"


# ════════════════════════════════════════════════════════════
# 1. Merge plasmid contigs ≥1 kb from all assemblies
# ════════════════════════════════════════════════════════════
echo "[$(date)] Merging plasmid contigs (≥${MIN_CONTIG_LEN} bp) from all assemblies..."

> "${OUTDIR}/merged/all_plasmid_contigs.fasta"

for SAMPLE_DIR in "${GENOMAD_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    PLASMID_FASTA="${SAMPLE_DIR}/${SAMPLE}_summary/${SAMPLE}_plasmid.fna"

    if [[ ! -f "${PLASMID_FASTA}" ]]; then
        echo "  WARNING: no plasmid FASTA found for ${SAMPLE}, skipping"
        continue
    fi

    # Filter ≥1 kb and prefix contig names with sample name
    awk -v sample="${SAMPLE}" -v min="${MIN_CONTIG_LEN}" \
        'BEGIN { header=""; seq="" }
         /^>/ {
             if (header != "" && length(seq) >= min) {
                 print ">"sample"__"substr(header,2)"\n"seq
             }
             header = $0; seq = ""
             next
         }
         { seq = seq $0 }
         END {
             if (header != "" && length(seq) >= min)
                 print ">"sample"__"substr(header,2)"\n"seq
         }' "${PLASMID_FASTA}" >> "${OUTDIR}/merged/all_plasmid_contigs.fasta"
done

N_MERGED=$(grep -c "^>" "${OUTDIR}/merged/all_plasmid_contigs.fasta")
echo "[$(date)] Total plasmid contigs merged: ${N_MERGED}"


# ════════════════════════════════════════════════════════════
# 2. Split into individual FASTAs for dRep
#    dRep expects one genome per FASTA file
# ════════════════════════════════════════════════════════════
echo "[$(date)] Splitting into per-contig FASTAs for dRep..."

awk -v outdir="${OUTDIR}/split_fastas" \
    '/^>/ {
        if (file != "") close(file)
        # sanitise contig name for use as filename
        name = substr($0, 2)
        gsub(/[^A-Za-z0-9._-]/, "_", name)
        file = outdir"/"name".fa"
        print > file
        next
    }
    { print > file }' \
    "${OUTDIR}/merged/all_plasmid_contigs.fasta"

N_SPLIT=$(ls "${OUTDIR}/split_fastas/"*.fa | wc -l)
echo "[$(date)] Split into ${N_SPLIT} individual FASTA files."


# ════════════════════════════════════════════════════════════
# 3. dRep dereplication → pOTUs
#    No CheckM2 quality filter — plasmids are not chromosomal,
#    so completeness/contamination metrics don't apply.
#    We pass --ignoreGenomeQuality to bypass the CheckM step.
# ════════════════════════════════════════════════════════════
conda activate drep_3.4.5
echo "[$(date)] Running dRep for plasmid dereplication..."
echo "[$(date)] Parameters: ANI=${DREP_ANI}, coverage=${DREP_COV}"

dRep dereplicate \
    "${OUTDIR}/pOTUs" \
    --genomes "${OUTDIR}/split_fastas/"*.fa \
    --ignoreGenomeQuality \
    --P_ani 0.85 \
    --S_ani "${DREP_ANI}" \
    --S_algorithm ANIn \
    --cov_thresh "${DREP_COV}" \
    --coverage_method larger \
    -p "${THREADS}"

N_POTUS=$(ls "${OUTDIR}/pOTUs/dereplicated_genomes/"*.fa 2>/dev/null | wc -l)
echo "[$(date)] pOTUs after dereplication: ${N_POTUS}"
echo "[$(date)] pOTU representatives: ${OUTDIR}/pOTUs/dereplicated_genomes/"
