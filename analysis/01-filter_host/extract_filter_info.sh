#!/bin/bash

# Directory containing the files
input_dir="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/silvereyes/metrics/"
# Output CSV file
output_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/silvereyes/metrics/extracted_metrics.tsv"

# Write the header to the output file (tab-separated)
echo -e "SAMPLE_NAME\tPRIMARY_READS\tREADS_AFTER_PREALIGNED_HOST_FILTER\tREADS_AFTER_QUALITY_AND_COMPLEXITY_FILTER\tREADS_AFTER_HOST_FILTER\tREADS_AFTER_DEDUPLICATION\tFINAL_PAIRED_READS\tFINAL_UNPAIRED_READS\tFINAL_TOTAL_READS\tLOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED\tHOST_READS_FILTERED\tDUPLICATE_READS_FILTERED" > "$output_file"

# Loop through each file in the directory
for filepath in "$input_dir"/*.bam; do
    # Extract the sample name from the file name
    filename=$(basename "$filepath")
    sample_name="${filename#metrics_}"
    sample_name="${sample_name%.sorted.bam}"

    # Extract the line containing the metrics
    metrics_line=$(grep -A 1 "PRIMARY_READS" "$filepath" | tail -n 1)

    # Check if metrics_line is not empty
    if [ -n "$metrics_line" ]; then
        # Append the sample name and the metrics line to the output file, using tabs
        echo -e "$sample_name\t$metrics_line" >> "$output_file"
    else
        echo "Metrics not found in $filename"
    fi
done

echo "Extraction complete! Data saved in $output_file."
