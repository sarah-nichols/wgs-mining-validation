#!/bin/sh

#SBATCH --job-name=kraken_output
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi

cd "$CONIFER"

# Create an empty file to store the final output
> $KRAKEN_OUTPUT_PAIRED


# Print the header line
echo -e "sample_id\ttaxon_name\ttaxid\treads\tP25\tP50\tP75" > $KRAKEN_OUTPUT_PAIRED

# Loop through all files in the directory
for file in $KRAKEN_PAIRED_STATS/*.report
do
  # Run conifer on the file
  ./conifer --rtl -s -i "$file" -d $KRAKEN_CUSTOM_TAXID > output.txt

  # Get the filename without the path
  filename=$(basename "$file")

  # Extract the sample ID from the filename
  sample_id=$(echo $filename | awk -F'_' '{print $2}' | awk -F'.' '{print $1}')

  # Add a column for the sample ID to the output file, skipping the header line
  awk -v id="$sample_id" 'NR>1 {printf("%s\t%s\n", id, $0)}' output.txt > output_with_id.txt

  # Append the output to the final output file
  cat output_with_id.txt >> $KRAKEN_OUTPUT_PAIRED
done