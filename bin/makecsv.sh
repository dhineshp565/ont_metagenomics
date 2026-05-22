#!/usr/bin/env bash

# This script makes a csv file with sample name and sample path with headers (SampleName,SamplePath).
# $1 = input path of fastq directories



input_dir="$1"
output_csv="samplelist.csv"

if [[ ! -d "$input_dir" ]]; then
    echo "Error: input directory not found: $input_dir" >&2
    exit 1
fi

# Write header (overwrite existing file for a fresh run).
echo "SampleName,SamplePath" > "$output_csv"

# Iterate over subdirectories in the input path.
for dir in "$input_dir"/*; do
    if [[ -d "$dir" ]]; then
        sample="$(basename "$dir")"
        samplename="${sample//_/-}"
        path="$(realpath "$dir")"
        echo "${samplename},${path}" >> "$output_csv"
    fi
done