#!/bin/bash
#SBATCH --job-name=mapping_stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --output=mapping_stats_%j.out
#SBATCH --error=mapping_stats_%j.err

# ============================================================
# Assembly mapping statistics
# ============================================================
# 1. % contigs from each assembly that ended up in dRep 95% MAGs
# 2. % reads (per sample) that map to dRep 95% MAGs
# ============================================================

set -euo pipefail

module load stack/.2024-06-silent gcc/12.2.0
module load samtools/1.16.1

MAG_DIR="04_drep_95_wilbd0/drep_output/dereplicated_genomes"
ASSEMBLIES_DIR="../Wildstrubel_assemblies"
OUTDIR="assembly_stats_wilbd0"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTDIR}"

ASSEMBLIES=(WIL1_before WIL1_field_before WIL2_before WIL2_after
            WILD_D0 WILB_D1 WILB_D3 WILB_D10)

# Raw read paths — used to calculate total reads before assembly
declare -A RAW_READS
RAW_READS["WIL1_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
RAW_READS["WIL1_field_before"]="[reads_dir_run3]/barcode17_filtered.fastq.gz"
RAW_READS["WIL2_before"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz"
RAW_READS["WIL2_after"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz"
RAW_READS["WILD_D0"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz"
RAW_READS["WILB_D1"]="[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
RAW_READS["WILB_D3"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz"
RAW_READS["WILB_D10"]="[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode17_filtered.fastq.gz"

# ════════════════════════════════════════════════════════════
# 0. Count total reads per sample from raw fastqs
# ════════════════════════════════════════════════════════════
echo "[$(date)] Counting total reads per sample..."

TOTAL_READS_OUT="${OUTDIR}/total_reads_per_sample.tsv"
echo -e "sample	total_reads" > "${TOTAL_READS_OUT}"

declare -A TOTAL_READ_COUNTS
for ASM in "${ASSEMBLIES[@]}"; do
    FASTQ="${RAW_READS[${ASM}]}"
    if [[ -f "${FASTQ}" ]]; then
        N=$(zcat "${FASTQ}" | awk 'NR%4==1' | wc -l)
        TOTAL_READ_COUNTS["${ASM}"]="${N}"
        echo -e "${ASM}	${N}"
        echo "[$(date)] ${ASM}: ${N} total reads"
    else
        TOTAL_READ_COUNTS["${ASM}"]="NA"
        echo -e "${ASM}	NA"
        echo "[$(date)] WARNING: raw reads not found for ${ASM}: ${FASTQ}"
    fi
done >> "${TOTAL_READS_OUT}"

# ════════════════════════════════════════════════════════════
# 1. Build contig-to-MAG lookup from 95% dereplicated MAGs
# ════════════════════════════════════════════════════════════
echo "[$(date)] Building contig-to-MAG lookup from dRep 95% MAGs..."

LOOKUP="${OUTDIR}/contig_to_mag_95.tsv"
> "${LOOKUP}"
for MAG in "${MAG_DIR}"/*.fa; do
    [[ -f "${MAG}" ]] || continue
    MAG_NAME=$(basename "${MAG}" .fa)
    grep "^>" "${MAG}" | sed 's/^>//' | awk -v mag="${MAG_NAME}" '{print $1"\t"mag}'
done > "${LOOKUP}"
echo "[$(date)] Total binned contigs: $(wc -l < ${LOOKUP})"

# ════════════════════════════════════════════════════════════
# 2. % contigs per assembly in dRep 95% MAGs
# ════════════════════════════════════════════════════════════
echo "[$(date)] Calculating contig binning rates..."

BINNING_OUT="${OUTDIR}/contig_binning_rates_95.tsv"
echo -e "assembly\ttotal_contigs\tbinned_contigs\tunbinned_contigs\tpct_contigs_binned\ttotal_bases\tbinned_bases\tpct_bases_binned" \
    > "${BINNING_OUT}"

for ASM in "${ASSEMBLIES[@]}"; do
    ASSEMBLY="${ASSEMBLIES_DIR}/${ASM}.fasta"
    [[ -f "${ASSEMBLY}" ]] || continue

    TOTAL_CONTIGS=0
    BINNED_CONTIGS=0
    TOTAL_BASES=0
    BINNED_BASES=0

    while IFS=$'\t' read -r CONTIG LEN; do
        CONTIG_CLEAN=$(echo "${CONTIG}" | sed 's/^>//' | awk '{print $1}')
        TOTAL_CONTIGS=$((TOTAL_CONTIGS + 1))
        TOTAL_BASES=$((TOTAL_BASES + LEN))
        IN_MAG=$(grep -m1 "^${CONTIG_CLEAN}"$'\t' "${LOOKUP}" 2>/dev/null || true)
        if [[ -n "${IN_MAG}" ]]; then
            BINNED_CONTIGS=$((BINNED_CONTIGS + 1))
            BINNED_BASES=$((BINNED_BASES + LEN))
        fi
    done < <(awk '/^>/{if(name) print name"\t"len; name=substr($0,2); len=0; next}{len+=length($0)}END{if(name) print name"\t"len}' "${ASSEMBLY}")

    UNBINNED_CONTIGS=$((TOTAL_CONTIGS - BINNED_CONTIGS))
    PCT_CONTIGS=$(awk "BEGIN{printf \"%.1f\", ${BINNED_CONTIGS}/${TOTAL_CONTIGS}*100}")
    PCT_BASES=$(awk "BEGIN{printf \"%.1f\", ${BINNED_BASES}/${TOTAL_BASES}*100}")

    echo -e "${ASM}\t${TOTAL_CONTIGS}\t${BINNED_CONTIGS}\t${UNBINNED_CONTIGS}\t${PCT_CONTIGS}\t${TOTAL_BASES}\t${BINNED_BASES}\t${PCT_BASES}"
    echo "[$(date)] ${ASM}: ${BINNED_CONTIGS}/${TOTAL_CONTIGS} contigs binned (${PCT_CONTIGS}%), ${PCT_BASES}% bases"

done >> "${BINNING_OUT}"

# ════════════════════════════════════════════════════════════
# 3. % reads mapping to dRep 95% MAGs per sample
# ════════════════════════════════════════════════════════════
echo "[$(date)] Calculating read mapping rates to dRep 95% MAGs..."

MAPPING_OUT="${OUTDIR}/read_mapping_rates_95.tsv"
echo -e "assembly\tsample\ttotal_raw_reads\treads_mapped_to_assembly\tpct_assembled\treads_to_MAGs\tpct_reads_to_MAGs\tpct_unassembled" \
    > "${MAPPING_OUT}"

for ASM in "${ASSEMBLIES[@]}"; do
    BAM_NAME="${ASM}"; [[ "${ASM}" == "WILD_D0" ]] && BAM_NAME="WILB_D0"; BAM_DIR="${BAM_NAME}_MAGs/01_mapping"
    ASSEMBLY="${ASSEMBLIES_DIR}/${ASM}.fasta"
    [[ -d "${BAM_DIR}" && -f "${ASSEMBLY}" ]] || continue

    # Get binned contig names for this assembly
    BINNED_LIST="${OUTDIR}/${ASM}_binned_contigs_95.txt"
    awk '/^>/{name=substr($0,2); gsub(/ .*/,"",name); print name}' "${ASSEMBLY}" | \
        while read -r CONTIG; do
            grep -m1 "^${CONTIG}"$'\t' "${LOOKUP}" > /dev/null 2>&1 && echo "${CONTIG}"
        done > "${BINNED_LIST}"

    N_BINNED=$(wc -l < "${BINNED_LIST}")
    echo "[$(date)] ${ASM}: ${N_BINNED} binned contigs to check"

    for BAM in "${BAM_DIR}"/*.bam; do
        [[ -f "${BAM}" ]] || continue
        SAMPLE=$(basename "${BAM}" .bam | sed "s/_vs_${ASM}//")

        FLAGSTAT=$(samtools flagstat -@ "${THREADS}" "${BAM}")
        TOTAL=$(echo    "${FLAGSTAT}" | awk 'NR==1{print $1}')
        MAPPED=$(echo   "${FLAGSTAT}" | awk '/^[0-9]+ \+ [0-9]+ mapped \(/{print $1}' | head -1)
        PCT_OVERALL=$(echo "${FLAGSTAT}" | awk '/^[0-9]+ \+ [0-9]+ mapped \(/{gsub(/[()%:]/,"",$5); print $5}' | head -1)

        # Count reads mapping specifically to binned contigs
        if [[ -s "${BINNED_LIST}" ]]; then
            READS_TO_MAGS=$(samtools view -@ "${THREADS}" -c "${BAM}" \
                $(head -100 "${BINNED_LIST}" | tr '\n' ' ') 2>/dev/null || echo 0)
        else
            READS_TO_MAGS=0
        fi

        PCT_TO_MAGS=$(awk "BEGIN{if(${MAPPED}>0) printf \"%.1f\", ${READS_TO_MAGS}/${MAPPED}*100; else print 0}")

        # Total raw reads for this sample (same assembly = same sample name)
        TOTAL_RAW="${TOTAL_READ_COUNTS[${ASM}]:-NA}"
        if [[ "${TOTAL_RAW}" != "NA" && "${TOTAL_RAW}" -gt 0 ]]; then
            PCT_ASSEMBLED=$(awk "BEGIN{printf \"%.1f\", ${MAPPED}/${TOTAL_RAW}*100}")
            PCT_UNASSEMBLED=$(awk "BEGIN{printf \"%.1f\", (${TOTAL_RAW}-${MAPPED})/${TOTAL_RAW}*100}")
        else
            PCT_ASSEMBLED="NA"
            PCT_UNASSEMBLED="NA"
        fi

        echo -e "${ASM}\t${SAMPLE}\t${TOTAL_RAW}\t${MAPPED}\t${PCT_ASSEMBLED}\t${READS_TO_MAGS}\t${PCT_TO_MAGS}\t${PCT_UNASSEMBLED}"
    done >> "${MAPPING_OUT}"
done

echo ""
echo "=== Complete ==="
echo "Contig binning rates: ${OUTDIR}/contig_binning_rates_95.tsv"
echo "Read mapping rates:   ${OUTDIR}/read_mapping_rates_95.tsv"
