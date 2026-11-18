#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G  # 5X0 GB per CPU
#SBATCH --time=48:00:00
#SBATCH --job-name="bracken_gtdb"
#SBATCH --output="bracken_gtdb.out"
#SBATCH --error="bracken_gtdb.err"

##################### Bracken read taxonomic classification and abundance #####################

# ---- SETTINGS ----
THREADS=24
DB_PATH=[kraken2_gtdb_db]
BRACKEN_DB="$DB_PATH"
READS_DIR="[reads_dir]/fastq_trimmed/trimmed_filtered"
OUTPUT_DIR="Kraken_Bracken_gtdb"
READ_LENGTH=300   # Set this to your mean read length
KMER_LEN=300    # Should match Kraken2 build

# ---- ENVIRONMENT ----
source [conda_base]/etc/profile.d/conda.sh
conda activate kraken2_bracken

# ---- CHECK FOR BRACKEN DB ----
KMER_FILE="$BRACKEN_DB/database${READ_LENGTH}mers.kmer_distrib"

if [ ! -f "$KMER_FILE" ]; then
    echo "Bracken DB for ${READ_LENGTH} bp not found. Building..."
    bracken-build -d "$BRACKEN_DB" -t "$THREADS" -k "$KMER_LEN" -l "$READ_LENGTH"
else
    echo "Bracken DB already exists for ${READ_LENGTH} bp."
fi

# ---- RUN BRACKEN ----
LEVELS=(D P C O F G S)

for report in "$OUTPUT_DIR"/*_kraken2.report; do
    sample=$(basename "$report" _kraken2.report)

    for level in "${LEVELS[@]}"; do
        echo "Running Bracken on $sample at level $level..."
        bracken -d "$BRACKEN_DB" \
                -i "$report" \
                -o "$OUTPUT_DIR/${sample}_bracken_${level}.tsv" \
                -r "$READ_LENGTH" \
                -l "$level" \
                -t "$THREADS"
    done
done
