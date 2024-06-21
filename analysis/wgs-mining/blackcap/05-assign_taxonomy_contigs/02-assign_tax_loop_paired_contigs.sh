#!/bin/sh

#SBATCH --job-name=kraken_assign_paired_contigs
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk
#SBATCH --array=1-388%10

set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi

module purge
module load Kraken2/2.1.1-gompi-2020b

# Get an array of all the fasta files
fasta_files=($OUTPUT_MEGA_PAIRED/*.fa)

# Get the fasta file for this array task
fasta_file=${fasta_files[$SLURM_ARRAY_TASK_ID-1]}

# Extract the sample ID from the filename
sample=$(basename "$fasta_file" .fa)

# Run kraken2 on the file
kraken2 --db "$KRAKEN_CAPPED_CUSTOM" \
    --threads 16 \
    --classified-out \
    --use-names \
    --output "$KRAKEN_PAIRED_CONTIGS_MEGA_STATS"/"$sample".stats.report \
    --report "$KRAKEN_PAIRED_CONTIGS_MEGA_REPORT"/"$sample".report \
    "$fasta_file"


