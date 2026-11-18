#!/bin/bash
#SBATCH --job-name=checkv_votu
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=checkv_votu_%j.out
#SBATCH --error=checkv_votu_%j.err

# ============================================================
# Step 7: CheckV quality filter + vOTU dereplication
# ============================================================
# 1. Merge viral contigs (≥1 kb) from all geNomad outputs
# 2. Run CheckV for quality assessment
# 3. Filter: remove contigs <3 kb or <3 genes; keep
#    medium-quality, high-quality, complete, and proviral contigs
# 4. Dereplicate with CD-HIT-EST at 95% ANI, 85% coverage → vOTUs
#
# TODO: confirm CD-HIT-EST environment name below.
#       Check availability with: conda search -c bioconda cd-hit
#       or ask your group if it is installed.
#
# Usage:
#   GENOMAD_DIR=/path/to/06_genomad \
#   sbatch 07_checkv_votu.sh
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

CDHIT_ENV="cdhit"

# ── Config ────────────────────────────────────────────────────
GENOMAD_DIR="${GENOMAD_DIR:-06_genomad}"
OUTDIR="${OUTDIR:-07_viruses}"
THREADS="${THREADS:-16}"

MIN_CONTIG_LEN=1000   # geNomad output filter (≥1 kb)
MIN_LENGTH_CHECKV=3000 # post-CheckV filter (≥3 kb)
MIN_GENES=3            # post-CheckV filter (≥3 genes)
CDHIT_ANI=0.95         # 95% nucleotide identity for vOTUs
CDHIT_COV=0.85         # 85% alignment coverage for vOTUs

mkdir -p "${OUTDIR}/merged" \
         "${OUTDIR}/checkv" \
         "${OUTDIR}/filtered" \
         "${OUTDIR}/vOTUs"


# ════════════════════════════════════════════════════════════
# 1. Merge viral contigs ≥1 kb from all assemblies
# ════════════════════════════════════════════════════════════
echo "[$(date)] Merging viral contigs (≥${MIN_CONTIG_LEN} bp) from all assemblies..."

> "${OUTDIR}/merged/all_viral_contigs.fasta"

for SAMPLE_DIR in "${GENOMAD_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    VIRAL_FASTA="${SAMPLE_DIR}/${SAMPLE}_summary/${SAMPLE}_virus.fna"

    if [[ ! -f "${VIRAL_FASTA}" ]]; then
        echo "  WARNING: no viral FASTA found for ${SAMPLE}, skipping"
        continue
    fi

    N_TOTAL=$(grep -c "^>" "${VIRAL_FASTA}")

    # Filter ≥1 kb and prefix contig names with sample name to avoid collisions
    # awk accumulates full sequence across multiple lines before checking length
    awk -v sample="${SAMPLE}" -v min="${MIN_CONTIG_LEN}"         'BEGIN { header=""; seq="" }
         /^>/ {
             if (header != "" && length(seq) >= min)
                 print ">"sample"__"substr(header,2)"\n"seq
             header = substr($0,2); seq = ""
             next
         }
         { seq = seq $0 }
         END {
             if (header != "" && length(seq) >= min)
                 print ">"sample"__"substr(header,2)"\n"seq
         }' "${VIRAL_FASTA}" >> "${OUTDIR}/merged/all_viral_contigs.fasta"

    N_KEPT=$(grep -c "^>" "${OUTDIR}/merged/all_viral_contigs.fasta" 2>/dev/null || echo 0)
    echo "  ${SAMPLE}: ${N_TOTAL} contigs in, kept ≥${MIN_CONTIG_LEN} bp"
done

N_MERGED=$(grep -c "^>" "${OUTDIR}/merged/all_viral_contigs.fasta")
echo "[$(date)] Total viral contigs merged: ${N_MERGED}"


# ════════════════════════════════════════════════════════════
# 2. CheckV quality assessment
# ════════════════════════════════════════════════════════════
conda activate checkV_1.0.3
echo "[$(date)] Running CheckV..."

# CHECKVDB should be set in the checkV_1.0.3 env; if not, set it here:
# export CHECKVDB="/path/to/checkv-db"

checkv end_to_end \
    "${OUTDIR}/merged/all_viral_contigs.fasta" \
    "${OUTDIR}/checkv" \
    -t "${THREADS}"

echo "[$(date)] CheckV complete."
conda deactivate


# ════════════════════════════════════════════════════════════
# 3. Filter: ≥3 kb, ≥3 genes, quality ≥ medium-quality
#    Retain: complete, high-quality, medium-quality, and
#    proviral contigs (which have already been trimmed by CheckV)
# ════════════════════════════════════════════════════════════
echo "[$(date)] Filtering viral contigs..."

# CheckV quality_summary.tsv columns:
#   1=contig_id, 3=contig_length, 6=gene_count, 8=checkv_quality, 9=miuvig_quality
# Quality categories: Complete, High-quality, Medium-quality, Low-quality, Not-determined

KEEP_LIST="${OUTDIR}/filtered/keep_contigs.txt"

awk -v min_len="${MIN_LENGTH_CHECKV}" -v min_genes="${MIN_GENES}" \
    'NR==1 { next }  # skip header
     {
         contig=$1; ctg_len=$3; genes=$6; quality=$8
         if (ctg_len >= min_len && genes >= min_genes &&
             (quality == "Complete" || quality == "High-quality" || quality == "Medium-quality")) {
             print contig
         }
     }' \
    "${OUTDIR}/checkv/quality_summary.tsv" > "${KEEP_LIST}"

# Also include proviruses (trimmed by CheckV, separate FASTA)
# Combine CheckV-trimmed proviruses + passing contigs into one FASTA
cat "${OUTDIR}/checkv/proviruses.fna" \
    "${OUTDIR}/checkv/viruses.fna" \
    > "${OUTDIR}/filtered/all_passing.fna.tmp"

# Extract only the kept contig IDs
conda activate checkV_1.0.3   # reuse env for seqtk if available, else use awk
if command -v seqtk &>/dev/null; then
    seqtk subseq \
        "${OUTDIR}/filtered/all_passing.fna.tmp" \
        "${KEEP_LIST}" \
        > "${OUTDIR}/filtered/viral_contigs_filtered.fasta"
else
    # Pure awk fallback
    awk 'NR==FNR { ids[$1]=1; next }
         /^>/ { keep = (substr($1,2) in ids) }
         keep { print }' \
        "${KEEP_LIST}" \
        "${OUTDIR}/filtered/all_passing.fna.tmp" \
        > "${OUTDIR}/filtered/viral_contigs_filtered.fasta"
fi
conda deactivate

rm "${OUTDIR}/filtered/all_passing.fna.tmp"

N_FILTERED=$(grep -c "^>" "${OUTDIR}/filtered/viral_contigs_filtered.fasta")
echo "[$(date)] Viral contigs passing quality filters: ${N_FILTERED}"


# ════════════════════════════════════════════════════════════
# 4. Dereplicate with CD-HIT-EST → vOTUs
# ════════════════════════════════════════════════════════════
if [[ -n "${CDHIT_ENV}" ]]; then
    conda activate "${CDHIT_ENV}"
    echo "[$(date)] Running CD-HIT-EST for vOTU dereplication..."
    echo "[$(date)] Parameters: ANI=${CDHIT_ANI}, coverage=${CDHIT_COV}"

    cd-hit-est \
        -i "${OUTDIR}/filtered/viral_contigs_filtered.fasta" \
        -o "${OUTDIR}/vOTUs/vOTUs.fasta" \
        -c "${CDHIT_ANI}" \
        -aS "${CDHIT_COV}" \
        -aL "${CDHIT_COV}" \
        -n 8 \
        -d 0 \
        -T "${THREADS}" \
        -M 0   # no memory limit (uses available RAM)

    N_VOTUS=$(grep -c "^>" "${OUTDIR}/vOTUs/vOTUs.fasta")
    echo "[$(date)] vOTUs after dereplication: ${N_VOTUS}"
    conda deactivate
else
    echo "[$(date)] CDHIT_ENV not set — skipping vOTU dereplication."
    echo "[$(date)] Run CD-HIT-EST manually on: ${OUTDIR}/filtered/viral_contigs_filtered.fasta"
fi

echo "[$(date)] Virus pipeline complete."
echo "[$(date)] vOTU representatives: ${OUTDIR}/vOTUs/vOTUs.fasta"
