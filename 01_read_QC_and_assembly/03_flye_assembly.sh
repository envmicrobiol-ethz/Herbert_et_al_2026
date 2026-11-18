#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --job-name="EM_LK_EXP001_meta_flye"
#SBATCH --mem-per-cpu=50G
#SBATCH --output="EM_LK_EXP001_meta_flye.out"
#SBATCH --error="EM_LK_EXP001_meta_flye.err"
#SBATCH --open-mode=append

# === CONFIG ===
THREADS=24
READ_DIR="[reads_dir]/fastq_trimmed/trimmed_filtered"
OUTDIR="[reads_dir]/flye_assemblies"
FLYE_LOGS="${OUTDIR}/logs"

# === ENVIRONMENT SETUP ===
source [conda_base]/etc/profile.d/conda.sh
mkdir -p "$OUTDIR" "$FLYE_LOGS"

# === LOOP OVER FASTQ FILES ===
for fq in "$READ_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    sample_out="${OUTDIR}/${sample}"
    echo "=== Processing sample: $sample ==="

    mkdir -p "$sample_out"

    # -------------------------------------------------------
    # 1️⃣  RUN FLYE
    # -------------------------------------------------------
    conda activate flye_V2.9
    echo "[`date`] Running Flye..."
    flye --nano-hq "$fq" \
         --threads "$THREADS" \
         --meta \
         --iterations 3 \
         --out-dir "$sample_out" \
         > "${FLYE_LOGS}/${sample}.flye.log" 2>&1
    echo "Flye finished for $sample"
done

echo "✅ All assemblies processed."

