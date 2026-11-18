#!/bin/bash
#SBATCH --job-name=dastool
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=dastool_%j.out
#SBATCH --error=dastool_%j.err

# ============================================================
# Step 3: DAS Tool — select best bins from all three binners
# ============================================================
# DAS Tool scores and dereplicates bins from MetaBAT2, CONCOCT,
# and MaxBin2, outputting a non-redundant refined bin set.
# Diamond is used as the gene-calling search engine (faster than BLAST).
#
# Usage:
#   ASSEMBLY=/path/to/assembly.fasta \
#   BINNING_DIR=/path/to/02_binning \
#   sbatch 03_dastool.sh
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate das_tool_1.1.7

# ── Config ────────────────────────────────────────────────────
ASSEMBLY="${ASSEMBLY:-../Wildstrubel_assemblies/WILB_D1.fasta}"
BINNING_DIR="${BINNING_DIR:-02_binning}"
OUTDIR="${OUTDIR:-03_dastool}"
THREADS="${THREADS:-16}"

ASSEMBLY_NAME=$(basename "${ASSEMBLY}" .fasta)

mkdir -p "${OUTDIR}/scaffolds2bin"

# ── Build scaffolds-to-bin tables ────────────────────────────
# Format required by DAS Tool: contig_id <TAB> bin_id
# We build these directly with awk rather than Fasta_to_Scaffolds2Bin.sh
# to avoid PATH issues and to explicitly exclude non-bin files.
#
# MetaBAT2 non-bin files to exclude:
#   bin.tooShort.fa, bin.lowDepth.fa, bin.unbinned.fa
# MaxBin2 non-bin files to exclude:
#   bin.tooshort, bin.noclass, bin.abundance, bin.marker*, bin.log etc.
#   (these lack the NNN.fasta pattern so the glob handles it)

echo "[$(date)] Generating scaffolds-to-bin tables..."

# Helper function: for each FASTA file, print contig_id <TAB> bin_name
fasta_to_s2b() {
    local FASTA_FILE="$1"
    local BIN_NAME
    BIN_NAME=$(basename "${FASTA_FILE}" | sed 's/\.\(fa\|fasta\)$//')
    grep "^>" "${FASTA_FILE}" | sed "s/^>//; s/ .*//" | awk -v bin="${BIN_NAME}" '{print $1"\t"bin}'
}
export -f fasta_to_s2b

# MetaBAT2: bin.N.fa — exclude tooShort, lowDepth, unbinned
> "${OUTDIR}/scaffolds2bin/metabat2.tsv"
for BIN in "${BINNING_DIR}/metabat2"/bin.[0-9]*.fa; do
    [[ -f "${BIN}" ]] || continue
    fasta_to_s2b "${BIN}" >> "${OUTDIR}/scaffolds2bin/metabat2.tsv"
done

# CONCOCT: N.fa (numeric bins only)
> "${OUTDIR}/scaffolds2bin/concoct.tsv"
for BIN in "${BINNING_DIR}/concoct/bins"/[0-9]*.fa; do
    [[ -f "${BIN}" ]] || continue
    fasta_to_s2b "${BIN}" >> "${OUTDIR}/scaffolds2bin/concoct.tsv"
done

# MaxBin2: bin.NNN.fasta — the numeric glob excludes tooshort/noclass/etc.
> "${OUTDIR}/scaffolds2bin/maxbin2.tsv"
for BIN in "${BINNING_DIR}/maxbin2"/bin.[0-9]*.fasta; do
    [[ -f "${BIN}" ]] || continue
    fasta_to_s2b "${BIN}" >> "${OUTDIR}/scaffolds2bin/maxbin2.tsv"
done

# Report line counts so we can verify all three are non-empty
echo "[$(date)] Scaffolds-to-bin table sizes:"
wc -l "${OUTDIR}/scaffolds2bin/metabat2.tsv" \
       "${OUTDIR}/scaffolds2bin/concoct.tsv" \
       "${OUTDIR}/scaffolds2bin/maxbin2.tsv"

# Sanity check — abort if any table is empty
for TSV in metabat2 concoct maxbin2; do
    if [[ ! -s "${OUTDIR}/scaffolds2bin/${TSV}.tsv" ]]; then
        echo "ERROR: ${TSV}.tsv is empty — check binning output"
        exit 1
    fi
done

# ── Run DAS Tool ─────────────────────────────────────────────
echo "[$(date)] Running DAS Tool..."

DAS_Tool \
    -i "${OUTDIR}/scaffolds2bin/metabat2.tsv","${OUTDIR}/scaffolds2bin/concoct.tsv","${OUTDIR}/scaffolds2bin/maxbin2.tsv" \
    -l MetaBAT2,CONCOCT,MaxBin2 \
    -c "${ASSEMBLY}" \
    -o "${OUTDIR}/${ASSEMBLY_NAME}_dastool" \
    --write_bins \
    --write_unbinned \
    --search_engine diamond \
    --score_threshold 0.4 \
    -t "${THREADS}"

N_BINS=$(ls "${OUTDIR}/${ASSEMBLY_NAME}_dastool_DASTool_bins/"*.fa 2>/dev/null | wc -l)
echo "[$(date)] DAS Tool complete. Refined bins: ${N_BINS}"
echo "[$(date)] Bins written to: ${OUTDIR}/${ASSEMBLY_NAME}_dastool_DASTool_bins/"
