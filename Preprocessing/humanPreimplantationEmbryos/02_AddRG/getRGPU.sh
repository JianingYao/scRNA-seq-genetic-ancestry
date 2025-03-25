#!/bin/bash

input_file="humanPreimplantationEmbryos_RG.tsv"
output_file="RGPU.tsv"
echo -e "Subdirectory\tRGPU" > "$output_file"

# Function to extract RGPU from a FASTQ file
extract_rgpu() {
    local fastq_file=$1
    # Extract the first read identifier and get the RGPU part
    zcat "$fastq_file" | head -n 1 | awk -F'_' '{print $1 "_" $2 "_" $3}'
}

while IFS=$'\t' read -r line; do
    subdir=$(echo "$line" | cut -f1)
    # Check if the directory exists
    if [[ -d "/scratch1/yaojiani/humanPreimplantationEmbryos/$subdir" ]]; then
        # Loop through each FASTQ file in the directory
        for fastq in "/scratch1/yaojiani/humanPreimplantationEmbryos/$subdir"/*.fastq.gz; do
            # Check if the file exists
            if [[ -f "$fastq" ]]; then
                rgpu=$(extract_rgpu "$fastq")
                echo -e "${subdir}\t${rgpu}" >> "$output_file"
            fi
        done
    else
        echo "Directory $subdir not found."
    fi
done < "$input_file"

echo "Processing complete. Output saved to $output_file"
