#!/bin/sh

#SBATCH --job-name=conifer_conf
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

cd "$CONIFER"

# Create an empty file to store the final output
> $CONIFER_OUTPUT_ZOST

# Print the header line
echo -e "sample_id\ttaxon_name\ttaxid\treads\tP25\tP50\tP75\tC25\tC50\tC75" > $CONIFER_OUTPUT_ZOST

# Loop through all files in the directory
for file in $KRAKEN_OUTPUT_ZOST_PAIRED/*.report
do
  # Run conifer on the file
  ./conifer --both_scores -s -i "$file" -d $KRAKEN_CUSTOM_TAXID > output.txt
  ./conifer --both_scores -i test_files/example.out.txt -d test_files/taxo.k2d


  # Get the filename without the path
  filename=$(basename "$file")

  # Extract the sample ID from the filename
  sample_id=$(echo $filename | awk '{match($0, /^[A-Za-z]{3}[0-9]{1,2}/, arr); print arr[0]}')

  echo "Filename: $filename"
  echo "Sample ID: $sample_id"
  # Add a column for the sample ID to the output file, skipping the header line
  awk -v id="$sample_id" 'NR>1 {printf("%s\t%s\n", id, $0)}' output.txt > output_with_id.txt

  # Append the output to the final output file
  cat output_with_id.txt >> $CONIFER_OUTPUT_ZOST
done

