#!/bin/bash

# Directory containing the kraken report files
input_dir="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/silvereyes/kraken_report_paired"
# Output file
output_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/silvereyes/extracted_unassembled_kraken_reads.tsv"

# Write the header to the output file
echo -e "SAMPLE_NAME\tUNIDENTIFIED_READS\tIDENTIFIED_READS" > "$output_file"

# Loop through each file in the directory
for filepath in "$input_dir"/*.report; do
    # Extract the sample name from the file name
    filename=$(basename "$filepath")
    sample_name="${filename%.sorted.kraken_report}"

    # Extract the unidentified and identified reads from the file
    unidentified_reads=$(awk 'NR==1 {print $2}' "$filepath")
    identified_reads=$(awk 'NR==2 {print $2}' "$filepath")

    # Append the sample name and the extracted values to the output file
    echo -e "$sample_name\t$unidentified_reads\t$identified_reads" >> "$output_file"
done

echo "Extraction complete! Data saved in $output_file."
