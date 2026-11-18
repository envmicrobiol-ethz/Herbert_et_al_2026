#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --time=48:00:00
#SBATCH --job-name="EM_LK_EXP002_ont_qc_trimming"
#SBATCH --mem-per-cpu=50G # Depending on the size of your data this may need to be increased
#SBATCH --output="EM_LK_EXP002_ont_qc_trimming.out"
#SBATCH --error="EM_LK_EXP002_ont_qc_trimming.err"
#SBATCH --open-mode=append

##################### Nanopore long-read QC and trimming pipeline #####################

 ## This script performs QC on your raw data using MinIONQC and Nanostat then trims adapters and remaining barcodes using Porechop, filters using Chopper for q score and length. The script thwen uses Minimap2 which will remove any DNA mapping to the lambda control used in library preparation, and then these files go through QC again using Nanostat to check the quality has improved. In the end it will remove the intermediate files and recompress the raw files and the new trimmed and filtered files for downstream use. The nanostat_combine script can be run to create a table to see how quality has improved and how much data was lost in filtered - it is good to check this to make sure it is not too much.

##################### CONFIGURABLE VARIABLES #####################

# This section contains the variables you need to change according to the location of your data and what filtering parameters you want to use.

basedir="[reads_dir]"
#summary_file="${basedir}/sequencing_summary_*.txt"
contaminant_db="[contamination_db]/ONT_lambda_DCS.fasta"
NS_OUT="${basedir}/Nanostat"

MINLEN=100
MAXLEN=""
QSCORE=15
THREADS=16

RAW_DIR="${basedir}/fastq_cat"
TRIMMED_DIR="${basedir}/fastq_trimmed"
FILTERED_DIR="${TRIMMED_DIR}/trimmed_filtered"

mkdir -p "$TRIMMED_DIR" "$FILTERED_DIR" "$NS_OUT/raw" "$NS_OUT/trimmed"

##################### PIPELINE START #####################

cd "$basedir"

# Decompress raw .fastq.gz files
echo "Decompressing .fastq.gz files..."
gunzip -f ${RAW_DIR}/*.fastq.gz

# 1. MinIONQC - this can be useful to problem shoot your run if there were issues but is not 100% necessary
#module load stack/2024-06 r/4.4.0
#Rscript [software_dir]/MinIONQC.R -i "$summary_file" -o MinIONQC

# 2. NanoStat on raw reads
source [conda_base]/etc/profile.d/conda.sh
conda activate nanostat_V1.4.0

for fq in ${RAW_DIR}/*.fastq; do
    sample=$(basename "$fq" .fastq)
    NanoStat -t "$THREADS" --fastq "$fq" > "$NS_OUT/raw/${sample}_raw.nanostat.txt"
done

# 3. Adapter trimming with Porechop
module load stack/2024-06 gcc/12.2.0 python/3.11.6
cd [software_dir]/venvs
source porechop_V0.2.4/bin/activate
cd "$basedir"

for fq in ${RAW_DIR}/*.fastq; do
    sample=$(basename "$fq" .fastq)
    porechop -i "$fq" -o "${TRIMMED_DIR}/${sample}_trimmed.fastq"
done

# 4. Filtering with Chopper
conda activate chopper_V0.10.0b
maxlen_arg=""
[ -n "$MAXLEN" ] && maxlen_arg="--maxlength $MAXLEN"

for file in ${TRIMMED_DIR}/*_trimmed.fastq; do
    sample=$(basename "$file" _trimmed.fastq)
    chopper --input "$file" \
            --minlength "$MINLEN" \
            $maxlen_arg \
            -q "$QSCORE" \
            --threads "$THREADS" \
            > "${FILTERED_DIR}/${sample}_filtered.fastq"
done

# 5. Filter contaminants using minimap2
conda activate minimap2_V2.22
module load samtools/1.16.1

for fq in ${FILTERED_DIR}/*.fastq; do
    sample=$(basename "$fq" .fastq)
    bam_unmapped="${FILTERED_DIR}/${sample}_unmapped.bam"
    clean_fq="${FILTERED_DIR}/${sample}.fastq.gz"

    # Map to contaminants and retain only unmapped
    minimap2 -t "$THREADS" -ax map-ont "$contaminant_db" "$fq" | \
        samtools view -@ "$THREADS" -b -f 4 - > "$bam_unmapped"

    samtools fastq -@ "$THREADS" "$bam_unmapped" | gzip > "$clean_fq"
    rm -f "$fq" "$bam_unmapped"  # Remove original .fastq and intermediate BAM
done

# 6. NanoStat on final filtered reads
conda activate nanostat_V1.4.0
for fq in ${FILTERED_DIR}/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    NanoStat -t "$THREADS" --fastq "$fq" > "$NS_OUT/trimmed/${sample}_trimmed.nanostat.txt"
done

# 7. Cleanup
echo "Removing *_trimmed.fastq..."
rm -f ${TRIMMED_DIR}/*_trimmed.fastq

echo "Recompressing original raw fastqs..."
gzip -f ${RAW_DIR}/*.fastq
