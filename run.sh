#!/bin/bash
# This script is used to run the pipeline with all the sra_ids in the samples.txt file
# Usage: bash run.sh

# Check if the samples.txt file exists
if [ ! -f samples.txt ]; then
    echo "samples.txt file not found!"
    exit 1
fi

# Check if fasterq-dump is installed
if ! command -v fasterq-dump &> /dev/null; then
    echo "fasterq-dump could not be found. Please install SRA Toolkit."
    exit 1
fi

rm -rf ./data
mkdir -p ./data

# Read sample IDs
sra_ids=()
while IFS= read -r line || [ -n "$line" ]; do
    sra_id=$(echo "$line" | awk '{print $1}')
    sra_ids+=("$sra_id")
done < samples.txt

total=${#sra_ids[@]}
echo "Total samples: $total"

# Process each ID
count=0
for sra_id in "${sra_ids[@]}"; do
    count=$((count + 1))
    echo "checking if $sra_id is already processed..."
    # check if files exist in the data directory
    if [ -f "./data/$sra_id/anticp.csv" ]; then
        echo "$sra_id is already processed. Skipping..."
        continue
    fi

    echo "$count/$total - Downloading SRA file for $sra_id..."
    mkdir -p ./data/$sra_id
    fasterq-dump "$sra_id" --outdir ./data/$sra_id
    echo "Downloaded $sra_id to ./data/$sra_id"

    peptominer -i ./data/$sra_id -o ./data/$sra_id --kraken --threads 8 --map

    rm -rf ./data/$sra_id/*.fastq
    rm -rf ./work/*
    rm -rf ./.nextflow.log.*
done

rm -rf ./work
echo "$total/$total sequences processed."
