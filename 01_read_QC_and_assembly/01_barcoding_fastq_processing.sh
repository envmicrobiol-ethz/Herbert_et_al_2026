#!/bin/bash

# Set the input and output directories
input_dir="fastq_pass"
output_dir="fastq_cat"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each barcode folder
for barcode_dir in "$input_dir"/barcode*/; do
    # Get folder name without path (e.g., barcode01)
    barcode_name=$(basename "$barcode_dir")
    
    # Define output file path
    output_file="${output_dir}/${barcode_name}.fastq.gz"
    
    # Concatenate all fastq files in the barcode folder
    cat "${barcode_dir}"/*.fastq.gz > "$output_file"
    
    echo "Created: $output_file"
done
