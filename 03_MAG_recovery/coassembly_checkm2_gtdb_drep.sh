#!/bin/bash
#SBATCH --job-name=coasm_checkm2_drep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=72:00:00
#SBATCH --output=coasm_checkm2_drep_%j.out
#SBATCH --error=coasm_checkm2_drep_%j.err

# ============================================================
# CheckM2 + dRep + GTDB-Tk for co-assembly MAGs
# Only keeps novel MAGs not represented in existing 82
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

COASM_BINS="21_coassembly/MAGs/03_dastool/assembly_dastool_DASTool_bins"
EXISTING_MAGS="04_drep_95_wilbd0/drep_output/dereplicated_genomes"
OUTDIR="21_coassembly/MAGs"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}/checkm2" "${OUTDIR}/all_for_drep" "${OUTDIR}/gtdbtk"

# ════════════════════════════════════════════════════════════
# 1. CheckM2 on co-assembly bins
# ════════════════════════════════════════════════════════════
conda activate checkM2_1.0.2

if [[ -f "${OUTDIR}/checkm2/quality_report.tsv" ]]; then
    echo "[$(date)] [SKIP] CheckM2 already done"
else
    echo "[$(date)] Running CheckM2 on ${COASM_BINS}..."
    N=$(ls "${COASM_BINS}"/*.fa | wc -l)
    echo "  Bins: ${N}"

    checkm2 predict \
        --input "${COASM_BINS}" \
        --output-directory "${OUTDIR}/checkm2" \
        --extension fa \
        --threads "${THREADS}"
fi

N_PASS=$(tail -n +2 "${OUTDIR}/checkm2/quality_report.tsv" | \
    awk '$2 >= 50 && $3 <= 10' | wc -l)
echo "[$(date)] Bins passing QC (compl ≥50%, cont <10%): ${N_PASS}"

conda deactivate

# ════════════════════════════════════════════════════════════
# 2. dRep — combine with existing 82 MAGs and dereplicate
# ════════════════════════════════════════════════════════════
conda activate drep_3.4.5

echo "[$(date)] Collecting all bins for dRep..."

# Copy existing 82 MAGs
cp "${EXISTING_MAGS}"/*.fa "${OUTDIR}/all_for_drep/"
echo "  Existing MAGs: $(ls ${OUTDIR}/all_for_drep/*.fa | wc -l)"

# Copy co-assembly bins with COASM_ prefix
for FA in "${COASM_BINS}"/*.fa; do
    cp "${FA}" "${OUTDIR}/all_for_drep/COASM_$(basename ${FA})"
done
echo "  Total with co-assembly: $(ls ${OUTDIR}/all_for_drep/*.fa | wc -l)"

# Build genome_info — existing MAGs + co-assembly bins
GENOME_INFO="${OUTDIR}/genome_info_combined.csv"
echo "genome,completeness,contamination" > "${GENOME_INFO}"

# Existing MAG scores from previous CheckM2
tail -n +2 04_drep_95_wilbd0/drep_output/data_tables/genomeInfo.csv \
    >> "${GENOME_INFO}" 2>/dev/null || \
tail -n +2 04_drep/genome_info.csv >> "${GENOME_INFO}" 2>/dev/null || true

# Co-assembly bin scores
tail -n +2 "${OUTDIR}/checkm2/quality_report.tsv" | \
    awk -v dir="${OUTDIR}/all_for_drep" \
        'BEGIN{OFS=","} {print dir"/COASM_"$1".fa", $2, $3}' \
    >> "${GENOME_INFO}"

echo "[$(date)] Running dRep 95% ANI..."

if [[ -d "${OUTDIR}/drep_95/dereplicated_genomes" ]]; then
    echo "  [SKIP] dRep already done"
else
    dRep dereplicate \
        "${OUTDIR}/drep_95" \
        --genomes "${OUTDIR}/all_for_drep"/*.fa \
        --genomeInfo "${GENOME_INFO}" \
        --completeness 50 \
        --contamination 10 \
        --P_ani 0.9 \
        --S_ani 0.95 \
        --S_algorithm ANIn \
        --coverage_method larger \
        -p "${THREADS}"
fi

# Identify novel co-assembly MAGs
N_TOTAL=$(ls "${OUTDIR}/drep_95/dereplicated_genomes/"*.fa | wc -l)
N_COASM=$(ls "${OUTDIR}/drep_95/dereplicated_genomes/"COASM_*.fa 2>/dev/null | wc -l)
N_EXISTING=$((N_TOTAL - N_COASM))

echo ""
echo "=== dRep results ==="
echo "  Existing MAGs retained: ${N_EXISTING}"
echo "  Novel co-assembly MAGs: ${N_COASM}"
echo "  Total MAGs: ${N_TOTAL}"

if [[ "${N_COASM}" -gt 0 ]]; then
    echo ""
    echo "Novel co-assembly MAGs:"
    ls "${OUTDIR}/drep_95/dereplicated_genomes/"COASM_*.fa | xargs -I{} basename {}
fi

conda deactivate

# ════════════════════════════════════════════════════════════
# 3. GTDB-Tk on novel co-assembly MAGs only
# ════════════════════════════════════════════════════════════
if [[ "${N_COASM}" -eq 0 ]]; then
    echo "[$(date)] No novel MAGs — skipping GTDB-Tk"
    exit 0
fi

export GTDBTK_DATA_PATH="[gtdbtk_db]"
conda activate gtdb-tk_2.4.0

echo "[$(date)] Running GTDB-Tk on ${N_COASM} novel co-assembly MAGs..."

mkdir -p "${OUTDIR}/gtdbtk_input"
cp "${OUTDIR}/drep_95/dereplicated_genomes/"COASM_*.fa "${OUTDIR}/gtdbtk_input/"

gtdbtk classify_wf \
    --genome_dir "${OUTDIR}/gtdbtk_input" \
    --out_dir "${OUTDIR}/gtdbtk" \
    --extension fa \
    --cpus "${THREADS}" \
    --skip_ani_screen

echo "[$(date)] GTDB-Tk complete."
echo ""
echo "=== Novel MAG taxonomy ==="
cat "${OUTDIR}/gtdbtk/gtdbtk.bac120.summary.tsv" 2>/dev/null | \
    awk -F'\t' 'NR>1{print $1"\t"$2}' | head -20
cat "${OUTDIR}/gtdbtk/gtdbtk.ar53.summary.tsv" 2>/dev/null | \
    awk -F'\t' 'NR>1{print $1"\t"$2}' | head -20

conda deactivate

echo ""
echo "[$(date)] Done."
echo "Novel MAGs:    ${OUTDIR}/drep_95/dereplicated_genomes/COASM_*.fa"
echo "GTDB-Tk:       ${OUTDIR}/gtdbtk/"
echo "CheckM2:       ${OUTDIR}/checkm2/quality_report.tsv"
