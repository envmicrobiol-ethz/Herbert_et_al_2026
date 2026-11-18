#!/bin/bash
# ============================================================
# Wildstrubel MAG pipeline — master submission script
# ============================================================
# Submits steps 1–3 per assembly, then steps 4–5 across all.
#
# Run this from your MAGs/ working directory:
#   cd [MAG_dir]
#   bash run_pipeline.sh
#
# Each step is submitted as a SLURM job with --dependency
# so the chain runs automatically.
# ============================================================

set -euo pipefail

SCRIPTS_DIR="$(dirname "$0")"
ASSEMBLIES_DIR="../Wildstrubel_assemblies"
SAMPLE_MAP="${SCRIPTS_DIR}/sample_map.tsv"

# All assemblies to process
ASSEMBLIES=(
    WIL1_before
    WIL1_field_before
    WIL2_after
    WIL2_before
    WILB_D0
    WILB_D1
    WILB_D3
    WILB_D10
)

DASTOOL_DIRS=""   # will be built up as jobs are submitted

for ASM in "${ASSEMBLIES[@]}"; do
    echo "=== Submitting pipeline for: ${ASM} ==="

    WORKDIR="${ASM}_MAGs"
    mkdir -p "${WORKDIR}"

    ASSEMBLY_PATH="${ASSEMBLIES_DIR}/${ASM}.fasta"

    # Step 1: Read mapping
    JID1=$(sbatch \
        --job-name="mapping_${ASM}" \
        --output="${WORKDIR}/01_mapping_%j.out" \
        --error="${WORKDIR}/01_mapping_%j.err" \
        --export=ALL,ASSEMBLY="${ASSEMBLY_PATH}",SAMPLE_MAP="${SAMPLE_MAP}",OUTDIR="${WORKDIR}/01_mapping" \
        "${SCRIPTS_DIR}/01_read_mapping.sh" \
        | awk '{print $NF}')
    echo "  Step 1 (mapping):  job ${JID1}"

    # Step 2: Binning — depends on step 1
    JID2=$(sbatch \
        --job-name="binning_${ASM}" \
        --dependency=afterok:"${JID1}" \
        --output="${WORKDIR}/02_binning_%j.out" \
        --error="${WORKDIR}/02_binning_%j.err" \
        --export=ALL,ASSEMBLY="${ASSEMBLY_PATH}",MAPPING_DIR="${WORKDIR}/01_mapping",OUTDIR="${WORKDIR}/02_binning" \
        "${SCRIPTS_DIR}/02_binning.sh" \
        | awk '{print $NF}')
    echo "  Step 2 (binning):  job ${JID2}"

    # Step 3: DAS Tool — depends on step 2
    JID3=$(sbatch \
        --job-name="dastool_${ASM}" \
        --dependency=afterok:"${JID2}" \
        --output="${WORKDIR}/03_dastool_%j.out" \
        --error="${WORKDIR}/03_dastool_%j.err" \
        --export=ALL,ASSEMBLY="${ASSEMBLY_PATH}",BINNING_DIR="${WORKDIR}/02_binning",OUTDIR="${WORKDIR}/03_dastool" \
        "${SCRIPTS_DIR}/03_dastool.sh" \
        | awk '{print $NF}')
    echo "  Step 3 (dastool):  job ${JID3}"

    DASTOOL_DIRS="${DASTOOL_DIRS} ${WORKDIR}/03_dastool"
done

# Collect all DAS Tool job IDs to use as dependency for step 4
ALL_JID3=$(squeue --me --format="%j %i" | grep "^dastool_" | awk '{print $2}' | tr '\n' ':' | sed 's/:$//')
echo ""
echo "=== Submitting dRep (step 4) across all assemblies ==="

JID4=$(sbatch \
    --job-name="drep_all" \
    --dependency=afterok:"${ALL_JID3}" \
    --output="04_drep_%j.out" \
    --error="04_drep_%j.err" \
    --export=ALL,DASTOOL_DIRS="${DASTOOL_DIRS}",OUTDIR="04_drep" \
    "${SCRIPTS_DIR}/04_drep.sh" \
    | awk '{print $NF}')
echo "  Step 4 (drep):     job ${JID4}"

echo ""
echo "=== Submitting GTDB-Tk + CoverM (step 5) ==="

JID5=$(sbatch \
    --job-name="annotation_all" \
    --dependency=afterok:"${JID4}" \
    --output="05_annotation_%j.out" \
    --error="05_annotation_%j.err" \
    --export=ALL,MAG_DIR="04_drep/drep_output/dereplicated_genomes",SAMPLE_MAP="${SAMPLE_MAP}",OUTDIR="05_annotation" \
    "${SCRIPTS_DIR}/05_gtdbtk_coverm.sh" \
    | awk '{print $NF}')
echo "  Step 5 (annotation): job ${JID5}"

echo ""
echo "Pipeline submitted. Monitor with: squeue --me"
