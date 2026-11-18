#!/bin/bash
#SBATCH --job-name=coassembly_mags
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=50G
#SBATCH --time=300:00:00
#SBATCH --output=coassembly_mags_%j.out
#SBATCH --error=coassembly_mags_%j.err

# ============================================================
# MAG recovery from Wildstrubel co-assembly
# ============================================================
# Co-assembly: 21_coassembly/assembly/assembly.fasta
# Maps all 8 sample reads for differential coverage binning
# Steps: minimap2 → MetaBAT2+CONCOCT+MaxBin2 → DAS Tool
#        → CheckM2 → dRep 95% ANI
#        → Compare against existing 82 MAGs
# ============================================================

set -euo pipefail

CONDA_BASE="[conda_base]"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

module load stack/.2024-06-silent gcc/12.2.0
module load samtools/1.16.1

COASM="21_coassembly/assembly/assembly.fasta"
OUTDIR="21_coassembly/MAGs"
EXISTING_MAGS="04_drep_95_wilbd0/drep_output/dereplicated_genomes"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUTDIR}/01_mapping" "${OUTDIR}/02_binning/metabat2" \
         "${OUTDIR}/02_binning/concoct/output" "${OUTDIR}/02_binning/maxbin2" \
         "${OUTDIR}/03_dastool/scaffolds2bin" "${OUTDIR}/checkm2" \
         "${OUTDIR}/drep_95" "${OUTDIR}/stats"

# ── Sample reads ──────────────────────────────────────────────
declare -A READS
READS["WIL1_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
READS["WIL2_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz"
READS["WIL2_after"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz"
READS["WILB_D3"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
READS["WILB_D10"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode17_filtered.fastq.gz"
READS["WILB_D0"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
READS["WILB_D1"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"

SAMPLES=(WIL1_before WIL2_before WIL2_after WILB_D0 WILB_D1 WILB_D3 WILB_D10)

# ════════════════════════════════════════════════════════════
# 1. Map all reads to co-assembly
# ════════════════════════════════════════════════════════════
conda activate minimap2_V2.22

echo "[$(date)] === Step 1: Read mapping ==="

# Index co-assembly once
if [[ ! -f "${COASM}.fai" ]]; then
    echo "  Indexing co-assembly..."
    samtools faidx "${COASM}"
fi

BAM_LIST=()
for SAMPLE in "${SAMPLES[@]}"; do
    READS_FILE="${READS[${SAMPLE}]}"
    BAM="${OUTDIR}/01_mapping/${SAMPLE}_vs_coasm.bam"
    BAM_LIST+=("${BAM}")

    if [[ ! -f "${READS_FILE}" ]]; then
        echo "  [SKIP] ${SAMPLE} — reads not found"
        continue
    fi
    if [[ -f "${BAM}" ]]; then
        echo "  [SKIP] ${SAMPLE} — BAM exists"
        continue
    fi

    echo "  [$(date '+%H:%M:%S')] Mapping ${SAMPLE}..."
    minimap2 -ax map-ont -t "${THREADS}" --secondary=no --split-prefix "${OUTDIR}/01_mapping/${SAMPLE}_tmp" \
        "${COASM}" "${READS_FILE}" | \
        samtools sort -@ "${THREADS}" -o "${BAM}" -
    samtools index "${BAM}"
    echo "  ${SAMPLE}: done"
done

conda deactivate

# ════════════════════════════════════════════════════════════
# 2. Binning
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 2: Binning ==="

# ── MetaBAT2 ──────────────────────────────────────────────────
conda activate metaBAT2_2.17

DEPTH="${OUTDIR}/02_binning/depth.txt"
if [[ ! -f "${DEPTH}" ]]; then
    echo "  Generating depth file..."
    jgi_summarize_bam_contig_depths \
        --outputDepth "${DEPTH}" \
        "${OUTDIR}/01_mapping/"*.bam
fi

if [[ ! $(ls "${OUTDIR}/02_binning/metabat2/"*.fa 2>/dev/null) ]]; then
    echo "  Running MetaBAT2..."
    metabat2 \
        -i "${COASM}" \
        -a "${DEPTH}" \
        -o "${OUTDIR}/02_binning/metabat2/bin" \
        -m 1500 \
        -t "${THREADS}" \
        --unbinned 2>/dev/null || true
    echo "  MetaBAT2 bins: $(ls ${OUTDIR}/02_binning/metabat2/bin*.fa 2>/dev/null | wc -l)"
fi

conda deactivate

# ── CONCOCT ───────────────────────────────────────────────────
conda activate concoct_1.1.0

if [[ ! -f "${OUTDIR}/02_binning/concoct/output/clustering_merged.csv" ]]; then
    echo "  Running CONCOCT..."
    cut_up_fasta.py "${COASM}" -c 10000 -o 0 --merge_last \
        -b "${OUTDIR}/02_binning/concoct/contigs_10k.bed" \
        > "${OUTDIR}/02_binning/concoct/contigs_10k.fasta"

    concoct_coverage_table.py \
        "${OUTDIR}/02_binning/concoct/contigs_10k.bed" \
        "${OUTDIR}/01_mapping/"*.bam \
        > "${OUTDIR}/02_binning/concoct/coverage_table.tsv"

    concoct \
        --composition_file "${OUTDIR}/02_binning/concoct/contigs_10k.fasta" \
        --coverage_file "${OUTDIR}/02_binning/concoct/coverage_table.tsv" \
        -b "${OUTDIR}/02_binning/concoct/output/" \
        -t "${THREADS}" 2>/dev/null || true

    merge_cutup_clustering.py \
        "${OUTDIR}/02_binning/concoct/output/clustering_gt1000.csv" \
        > "${OUTDIR}/02_binning/concoct/output/clustering_merged.csv"

    extract_fasta_bins.py "${COASM}" \
        "${OUTDIR}/02_binning/concoct/output/clustering_merged.csv" \
        --output_path "${OUTDIR}/02_binning/concoct/output/bins/" 2>/dev/null || true

    echo "  CONCOCT bins: $(ls ${OUTDIR}/02_binning/concoct/output/bins/*.fa 2>/dev/null | wc -l)"
fi

conda deactivate

# ── MaxBin2 ───────────────────────────────────────────────────
conda activate maxbin2_2.2.7

if [[ ! $(ls "${OUTDIR}/02_binning/maxbin2/"*.fasta 2>/dev/null) ]]; then
    echo "  Running MaxBin2..."
    awk '{print $1"\t"$4}' "${DEPTH}" | tail -n +2 \
        > "${OUTDIR}/02_binning/maxbin2/depth_maxbin.txt"

    run_MaxBin.pl \
        -contig "${COASM}" \
        -abund "${OUTDIR}/02_binning/maxbin2/depth_maxbin.txt" \
        -out "${OUTDIR}/02_binning/maxbin2/bin" \
        -thread "${THREADS}" 2>/dev/null || true

    echo "  MaxBin2 bins: $(ls ${OUTDIR}/02_binning/maxbin2/bin*.fasta 2>/dev/null | wc -l)"
fi

conda deactivate

# ════════════════════════════════════════════════════════════
# 3. DAS Tool
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 3: DAS Tool ==="
conda activate das_tool_1.1.7

DAS_OUT="${OUTDIR}/03_dastool/coasm_dastool"

if [[ -d "${DAS_OUT}_DASTool_bins" ]]; then
    echo "  [SKIP] DAS Tool already done"
else
    # Build scaffolds2bin tables
    Fasta_to_Scaffolds2Bin.sh \
        -i "${OUTDIR}/02_binning/metabat2" -e fa \
        > "${OUTDIR}/03_dastool/scaffolds2bin/metabat2.tsv" 2>/dev/null || true

    Fasta_to_Scaffolds2Bin.sh \
        -i "${OUTDIR}/02_binning/concoct/output/bins" -e fa \
        > "${OUTDIR}/03_dastool/scaffolds2bin/concoct.tsv" 2>/dev/null || true

    Fasta_to_Scaffolds2Bin.sh \
        -i "${OUTDIR}/02_binning/maxbin2" -e fasta \
        > "${OUTDIR}/03_dastool/scaffolds2bin/maxbin2.tsv" 2>/dev/null || true

    # Build input lists (skip empty files)
    BINS_INPUT=""
    LABELS_INPUT=""
    for BINNER in metabat2 concoct maxbin2; do
        TSV="${OUTDIR}/03_dastool/scaffolds2bin/${BINNER}.tsv"
        if [[ -s "${TSV}" ]]; then
            BINS_INPUT="${BINS_INPUT:+${BINS_INPUT},}${TSV}"
            LABELS_INPUT="${LABELS_INPUT:+${LABELS_INPUT},}${BINNER}"
        fi
    done

    DAS_Tool \
        -i "${BINS_INPUT}" \
        -l "${LABELS_INPUT}" \
        -c "${COASM}" \
        -o "${DAS_OUT}" \
        --write_bins \
        --threads "${THREADS}"

    echo "  DAS Tool bins: $(ls ${DAS_OUT}_DASTool_bins/*.fa 2>/dev/null | wc -l)"
fi

conda deactivate

# ════════════════════════════════════════════════════════════
# 4. CheckM2
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 4: CheckM2 ==="
conda activate checkM2_1.0.2

if [[ -f "${OUTDIR}/checkm2/quality_report.tsv" ]]; then
    echo "  [SKIP] CheckM2 already done"
else
    checkm2 predict \
        --input "${DAS_OUT}_DASTool_bins" \
        --output-directory "${OUTDIR}/checkm2" \
        --extension fa \
        --threads "${THREADS}"
fi

N_PASS=$(tail -n +2 "${OUTDIR}/checkm2/quality_report.tsv" | \
    awk '$2 >= 50 && $3 <= 10' | wc -l)
echo "  Bins passing QC (compl ≥50%, cont <10%): ${N_PASS}"

conda deactivate

# ════════════════════════════════════════════════════════════
# 5. dRep 95% ANI — co-assembly bins only
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 5: dRep 95% ANI ==="
conda activate drep_3.4.5

# Build genome_info
GENOME_INFO="${OUTDIR}/genome_info.csv"
echo "genome,completeness,contamination" > "${GENOME_INFO}"
tail -n +2 "${OUTDIR}/checkm2/quality_report.tsv" | \
    awk -v dir="${DAS_OUT}_DASTool_bins" \
        'BEGIN{OFS=","} {print dir"/"$1".fa", $2, $3}' \
    >> "${GENOME_INFO}"

if [[ -d "${OUTDIR}/drep_95/dereplicated_genomes" ]]; then
    echo "  [SKIP] dRep already done"
else
    dRep dereplicate \
        "${OUTDIR}/drep_95" \
        --genomes "${DAS_OUT}_DASTool_bins"/*.fa \
        --genomeInfo "${GENOME_INFO}" \
        --completeness 50 \
        --contamination 10 \
        --P_ani 0.9 \
        --S_ani 0.95 \
        --S_algorithm ANIn \
        --coverage_method larger \
        -p "${THREADS}"
fi

N_DREP=$(ls "${OUTDIR}/drep_95/dereplicated_genomes/"*.fa 2>/dev/null | wc -l)
echo "  Co-assembly MAGs (95% ANI): ${N_DREP}"

conda deactivate

# ════════════════════════════════════════════════════════════
# 6. Compare with existing 82 MAGs
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 6: Compare with existing MAGs ==="
conda activate drep_3.4.5

# Combine existing 82 MAGs + new co-assembly MAGs and dereplicate
COMBINED_DIR="${OUTDIR}/all_for_drep"
mkdir -p "${COMBINED_DIR}"

# Copy existing MAGs
cp "${EXISTING_MAGS}"/*.fa "${COMBINED_DIR}/"
echo "  Existing MAGs: $(ls ${COMBINED_DIR}/*.fa | wc -l)"

# Copy co-assembly MAGs with prefix
for FA in "${OUTDIR}/drep_95/dereplicated_genomes/"*.fa; do
    cp "${FA}" "${COMBINED_DIR}/COASM_$(basename ${FA})"
done
echo "  After adding co-assembly MAGs: $(ls ${COMBINED_DIR}/*.fa | wc -l)"

# Build combined genome_info
COMBINED_INFO="${OUTDIR}/genome_info_combined.csv"
echo "genome,completeness,contamination" > "${COMBINED_INFO}"

# Existing MAG genome info
tail -n +2 04_drep_95_wilbd0/drep_output/data_tables/genomeInfo.csv \
    >> "${COMBINED_INFO}" 2>/dev/null || \
tail -n +2 04_drep/genome_info.csv >> "${COMBINED_INFO}" 2>/dev/null || true

# Co-assembly MAG genome info
tail -n +2 "${OUTDIR}/checkm2/quality_report.tsv" | \
    awk -v dir="${COMBINED_DIR}" \
        'BEGIN{OFS=","} {print dir"/COASM_"$1".fa", $2, $3}' \
    >> "${COMBINED_INFO}"

if [[ -d "${OUTDIR}/drep_95_combined/dereplicated_genomes" ]]; then
    echo "  [SKIP] Combined dRep already done"
else
    dRep dereplicate \
        "${OUTDIR}/drep_95_combined" \
        --genomes "${COMBINED_DIR}"/*.fa \
        --genomeInfo "${COMBINED_INFO}" \
        --completeness 50 \
        --contamination 10 \
        --P_ani 0.9 \
        --S_ani 0.95 \
        --S_algorithm ANIn \
        --coverage_method larger \
        -p "${THREADS}"
fi

N_COMBINED=$(ls "${OUTDIR}/drep_95_combined/dereplicated_genomes/"*.fa 2>/dev/null | wc -l)
N_NEW_COASM=$(ls "${OUTDIR}/drep_95_combined/dereplicated_genomes/"COASM_*.fa 2>/dev/null | wc -l)

echo ""
echo "=== Summary ==="
echo "  Co-assembly stats: 211,503 contigs, 4.1 Gb"
echo "  Co-assembly MAGs (95% ANI): ${N_DREP}"
echo "  Existing MAGs: 82"
echo "  Combined dRep 95% ANI: ${N_COMBINED}"
echo "  Novel co-assembly MAGs not in existing 82: ${N_NEW_COASM}"

echo "[$(date)] Done."
echo "Co-assembly MAGs:   ${OUTDIR}/drep_95/dereplicated_genomes/"
echo "Novel MAGs:         ${OUTDIR}/drep_95_combined/dereplicated_genomes/COASM_*.fa"
