#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=500000  # 500 GB per CPU
#SBATCH --time=48:00:00
#SBATCH --job-name="kraken2_gtdb"
#SBATCH --output="kraken2_gtdb.out"
#SBATCH --error="kraken2_gtdb.err"

##################### Kraken2 read taxonomic classification and abundance #####################

## This script runs Kraken2 to taxonomically classify the reads and then Bracken to estimate relative abundance.

# ---- SETTINGS ----
THREADS=4
DB_PATH=[kraken2_gtdb_db]
BRACKEN_DB=$DB_PATH
READS_DIR="[reads_dir]/fastq_trimmed/trimmed_filtered"
OUTPUT_DIR="[reads_dir]/Kraken_Bracken_gtdb"
READ_LENGTH=300   # set this to your mean read length after QC
KMER_LEN=300        # should match your Kraken2 DB

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ---- RUN KRAKEN2 ----

source [conda_base]/etc/profile.d/conda.sh
conda activate kraken2_bracken

for fq in "$READS_DIR"/*.fastq.gz; do
    [ -e "$fq" ] || continue  # skip if no match
    sample=$(basename "$fq" | sed 's/.fastq.gz$//;s/.fastq$//')

    echo "Running Kraken2 on $sample..."
    kraken2 --db "$DB_PATH" \
            --threads "$THREADS" \
            --gzip-compressed \
            --report "$OUTPUT_DIR/${sample}_kraken2.report" \
            --output "$OUTPUT_DIR/${sample}_kraken2.output" \
            --use-names \
            "$fq"
done

# ---- RUN BRACKEN ----

# ---- CHECK FOR BRACKEN DB ----
KMER_FILE="$BRACKEN_DB/database${READ_LENGTH}mers.kmer_distrib"

if [ ! -f "$KMER_FILE" ]; then
    echo "Bracken DB for ${READ_LENGTH} bp not found. Building..."
    bracken-build -d "$BRACKEN_DB" -t "$THREADS" -k "$KMER_LEN" -l "$READ_LENGTH"
else
    echo "Bracken DB already exists for ${READ_LENGTH} bp."
fi

for report in "$OUTPUT_DIR"/*_kraken2.report; do
    sample=$(basename "$report" _kraken2.report)

    echo "Running Bracken on $sample..."
    bracken -d "$BRACKEN_DB" \
            -i "$report" \
            -o "$OUTPUT_DIR/${sample}_bracken.tsv" \
            -r "$READ_LENGTH" \
            -l S \
            -t "$THREADS"
done
