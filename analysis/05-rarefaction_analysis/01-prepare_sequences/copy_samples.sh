#!/bin/bash

mkdir -p /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/rarefaction/own_samples/raw

# Path to the text file containing the list of .bam files
BAM_LIST="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/sample-info/novogene_sample_paths.txt"

# Destination directory for the files
DEST_DIR="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/rarefaction/own_samples/raw"

# Iterate over each line in the text file
while read -r BAM_PATH; do
  # Copy the .bam file
  cp "$BAM_PATH" "$DEST_DIR/"
  
  # Copy the corresponding .bam.bai index file
  BAI_PATH="${BAM_PATH%.bam}.bam.bai"
  if [ -f "$BAI_PATH" ]; then
    cp "$BAI_PATH" "$DEST_DIR/"
  else
    echo "Index file $BAI_PATH not found!"
  fi
done < "$BAM_LIST"
