#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=binning_%j.out
#SBATCH --error=binning_%j.err

# ============================================================
# Step 2: Binning with MetaBAT2, CONCOCT, and MaxBin2
# ============================================================
# Each tool lives in its own conda environment. We activate and
# deactivate in sequence: MetaBAT2 → CONCOCT → MaxBin2.
#
# Multi-sample differential coverage (all BAMs passed at once)
# is used wherever possible — the biggest quality lever for
# ONT metagenome binning.
#
# Usage:
#   ASSEMBLY=/path/to/assembly.fasta \
#   MAPPING_DIR=/path/to/01_mapping \
#   sbatch 02_binning.sh
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── Config ────────────────────────────────────────────────────
ASSEMBLY="${ASSEMBLY:-../Wildstrubel_assemblies/WILB_D1.fasta}"
MAPPING_DIR="${MAPPING_DIR:-01_mapping}"
OUTDIR="${OUTDIR:-02_binning}"
THREADS="${THREADS:-16}"
MIN_CONTIG_LEN=2000

ASSEMBLY_NAME=$(basename "${ASSEMBLY}" .fasta)

mkdir -p "${OUTDIR}/metabat2" \
         "${OUTDIR}/concoct/cut_contigs" \
         "${OUTDIR}/concoct/output" \
         "${OUTDIR}/maxbin2"

# Collect all BAMs for this assembly
BAM_FILES=( "${MAPPING_DIR}"/*_vs_"${ASSEMBLY_NAME}".bam )
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${MAPPING_DIR} for assembly ${ASSEMBLY_NAME}"
    exit 1
fi
echo "[$(date)] Found ${#BAM_FILES[@]} BAM file(s): ${BAM_FILES[*]}"


# ════════════════════════════════════════════════════════════
# MetaBAT2
# ════════════════════════════════════════════════════════════
conda activate metaBAT2_2.17
echo "[$(date)] [MetaBAT2] Generating contig depth file..."

jgi_summarize_bam_contig_depths \
    --outputDepth "${OUTDIR}/metabat2/depth.txt" \
    "${BAM_FILES[@]}"

echo "[$(date)] [MetaBAT2] Running binner..."

metabat2 \
    -i "${ASSEMBLY}" \
    -a "${OUTDIR}/metabat2/depth.txt" \
    -o "${OUTDIR}/metabat2/bin" \
    -m "${MIN_CONTIG_LEN}" \
    -t "${THREADS}" \
    --seed 42 \
    --unbinned

echo "[$(date)] [MetaBAT2] Done."
conda deactivate


# ════════════════════════════════════════════════════════════
# CONCOCT
# ════════════════════════════════════════════════════════════
conda activate concoct_1.1.0
echo "[$(date)] [CONCOCT] Cutting contigs to 10 kb chunks..."

CHUNK=10000
cut_up_fasta.py \
    "${ASSEMBLY}" \
    -c "${CHUNK}" \
    -o 0 \
    --merge_last \
    -b "${OUTDIR}/concoct/cut_contigs/contigs_10k.bed" \
    > "${OUTDIR}/concoct/cut_contigs/contigs_10k.fasta"

echo "[$(date)] [CONCOCT] Generating coverage table..."

concoct_coverage_table.py \
    "${OUTDIR}/concoct/cut_contigs/contigs_10k.bed" \
    "${BAM_FILES[@]}" \
    > "${OUTDIR}/concoct/coverage_table.tsv"

echo "[$(date)] [CONCOCT] Running binner..."

concoct \
    --composition_file "${OUTDIR}/concoct/cut_contigs/contigs_10k.fasta" \
    --coverage_file "${OUTDIR}/concoct/coverage_table.tsv" \
    -b "${OUTDIR}/concoct/output/" \
    -t "${THREADS}" \
    -s 42

echo "[$(date)] [CONCOCT] Merging clustering results and extracting bin FASTAs..."

merge_cutup_clustering.py \
    "${OUTDIR}/concoct/output/clustering_gt1000.csv" \
    > "${OUTDIR}/concoct/clustering_merged.csv"

mkdir -p "${OUTDIR}/concoct/bins"
extract_fasta_bins.py \
    "${ASSEMBLY}" \
    "${OUTDIR}/concoct/clustering_merged.csv" \
    --output_path "${OUTDIR}/concoct/bins"

echo "[$(date)] [CONCOCT] Done."
conda deactivate


# ════════════════════════════════════════════════════════════
# MaxBin2
# ════════════════════════════════════════════════════════════
conda activate maxbin2_2.2.7
echo "[$(date)] [MaxBin2] Preparing per-sample abundance files..."

# The MetaBAT2 depth file has the format:
#   contigName  contigLen  totalAvgDepth  sample1_depth  sample1_var  sample2_depth ...
# MaxBin2 wants one file per sample: contig <TAB> depth
ABUND_LIST="${OUTDIR}/maxbin2/abund_list.txt"
> "${ABUND_LIST}"

HEADER_ROW=$(head -1 "${OUTDIR}/metabat2/depth.txt")
N_SAMPLES=$(( ($(echo "${HEADER_ROW}" | awk '{print NF}') - 3) / 2 ))

for i in $(seq 1 "${N_SAMPLES}"); do
    SAMPLE_COL=$(( 3 + (i - 1) * 2 + 1 ))   # columns 4, 6, 8, ...
    ABUND_FILE="${OUTDIR}/maxbin2/abund_sample${i}.txt"
    awk -v col="${SAMPLE_COL}" 'NR>1 {print $1"\t"$col}' \
        "${OUTDIR}/metabat2/depth.txt" > "${ABUND_FILE}"
    echo "${ABUND_FILE}" >> "${ABUND_LIST}"
done

echo "[$(date)] [MaxBin2] Running binner..."

run_MaxBin.pl \
    -contig "${ASSEMBLY}" \
    -abund_list "${ABUND_LIST}" \
    -out "${OUTDIR}/maxbin2/bin" \
    -min_contig_length "${MIN_CONTIG_LEN}" \
    -thread "${THREADS}"

echo "[$(date)] [MaxBin2] Done."
conda deactivate

echo "[$(date)] All binning complete."
