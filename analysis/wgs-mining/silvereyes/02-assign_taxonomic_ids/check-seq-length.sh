#!/bin/bash

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

module load SAMtools/1.16.1-GCC-11.3.0

<<'COMMENT'
# Variable to hold total length and count of all reads
total_length=0
total_count=0

# Loop through each BAM file in the directory
for bam_file in "$BAM_ZOST_PATH"/paired/*.bam; do
  # Calculate the average length of reads in the current BAM file
  read_lengths=$(samtools view "$bam_file" | awk '{total += length($10); count++} END {print total, count}')
  
  # Extract total length and count from the output
  file_total_length=$(echo $read_lengths | cut -d' ' -f1)
  file_total_count=$(echo $read_lengths | cut -d' ' -f2)
  
  # Add to overall total length and count
  total_length=$((total_length + file_total_length))
  total_count=$((total_count + file_total_count))
  
  # Optionally, print average length for this file
  echo "Average length for $(basename "$bam_file"): $((file_total_length / file_total_count))"
done

# Calculate and print the overall average length
echo "Overall average length paired bam: $((total_length / total_count))"

# Loop through each BAM file in the directory
for bam_file in "$BAM_ZOST_PATH"/paired/*.bam; do
  # Count the reads in the current BAM file
  read_count=$(samtools view -c "$bam_file")
  
  # Print the file name and read count
  echo "$(basename "$bam_file"): $read_count reads"
done

# Variable to hold total length and count of all reads
total_length=0
total_count=0

# Loop through each BAM file in the directory
for bam_file in "$BAM_ZOST_PATH"/unpaired/*.bam; do
  # Calculate the average length of reads in the current BAM file
  read_lengths=$(samtools view "$bam_file" | awk '{total += length($10); count++} END {print total, count}')
  
  # Extract total length and count from the output
  file_total_length=$(echo $read_lengths | cut -d' ' -f1)
  file_total_count=$(echo $read_lengths | cut -d' ' -f2)
  
  # Add to overall total length and count
  total_length=$((total_length + file_total_length))
  total_count=$((total_count + file_total_count))
  
  # Optionally, print average length for this file
  echo "Average length for $(basename "$bam_file"): $((file_total_length / file_total_count))"
done

# Calculate and print the overall average length
echo "Overall average length unpaired bam: $((total_length / total_count))"

for bam_file in "$BAM_ZOST_PATH"/unpaired/*.bam; do
  # Count the reads in the current BAM file
  read_count=$(samtools view -c "$bam_file")
  
  # Print the file name and read count
  echo "$(basename "$bam_file"): $read_count reads"
done



# Initialize variables for overall totals
overall_total_length=0
overall_total_count=0

# Loop through each FASTQ file in the directory
for fastq_file in "$PAIRED_F_FASTQS_ZOST"/*.fastq; do
  # Process the FASTQ file with awk
  read_stats=$(awk '{if(NR%4==2) {total_length+=length($0); count++}} END {print total_length, count}' "$fastq_file")
  
  # Extract total length and count from the output
  file_total_length=$(echo $read_stats | cut -d' ' -f1)
  file_total_count=$(echo $read_stats | cut -d' ' -f2)
  
  # Add to overall totals
  overall_total_length=$((overall_total_length + file_total_length))
  overall_total_count=$((overall_total_count + file_total_count))
  
  # Print mean read length and count for this file
  echo "$(basename "$fastq_file"): Mean read length = $((file_total_length / file_total_count)), Count = $file_total_count"
done

# Calculate and print the overall mean read length and total count
echo "Overall: Mean read length forward = $((overall_total_length / overall_total_count)), Total count = $overall_total_count"

# Initialize variables for overall totals
overall_total_length=0
overall_total_count=0

# Loop through each FASTQ file in the directory
for fastq_file in "$PAIRED_R_FASTQS_ZOST"/*.fastq; do
  # Process the FASTQ file with awk
  read_stats=$(awk '{if(NR%4==2) {total_length+=length($0); count++}} END {print total_length, count}' "$fastq_file")
  
  # Extract total length and count from the output
  file_total_length=$(echo $read_stats | cut -d' ' -f1)
  file_total_count=$(echo $read_stats | cut -d' ' -f2)
  
  # Add to overall totals
  overall_total_length=$((overall_total_length + file_total_length))
  overall_total_count=$((overall_total_count + file_total_count))
  
  # Print mean read length and count for this file
  echo "$(basename "$fastq_file"): Mean read length = $((file_total_length / file_total_count)), Count = $file_total_count"
done

# Calculate and print the overall mean read length and total count
echo "Overall: Mean read length reverse = $((overall_total_length / overall_total_count)), Total count = $overall_total_count"

COMMENT
set +e

overall_total_count=0
overall_total_length=0

# Assuming this loop is for a specific set of FASTQ files
overall_total_count=0
overall_total_length=0

# Assuming this loop is for a specific set of FASTQ files
for fastq_file in "$ASSEMBLED_FASTQS_ZOST"/*.fastq; do
  if [ -f "$fastq_file" ]; then
    # Process the FASTQ file with awk, adding error handling
    read_stats=$(awk '{if(NR%4==2) {total_length+=length($0); count++}} END {print total_length, count}' "$fastq_file" 2>/dev/null)
    
    if [ -z "$read_stats" ]; then
      echo "Error processing $fastq_file. Skipping."
      continue
    fi

    file_total_length=$(echo $read_stats | cut -d' ' -f1)
    file_total_count=$(echo $read_stats | cut -d' ' -f2)
    
    if [ "$file_total_count" -gt 0 ]; then
      overall_total_length=$((overall_total_length + file_total_length))
      overall_total_count=$((overall_total_count + file_total_count))
      
      echo "$(basename "$fastq_file"): Mean read length = $((file_total_length / file_total_count)), Count = $file_total_count"
    else
      echo "$(basename "$fastq_file") has no reads. Skipping."
    fi
  else
    echo "$fastq_file does not exist. Skipping."
  fi
done

if [ "$overall_total_count" -gt 0 ]; then
  echo "Overall: Mean read length = $((overall_total_length / overall_total_count)), Total count = $overall_total_count"
else
  echo "No reads found in any files."
fi