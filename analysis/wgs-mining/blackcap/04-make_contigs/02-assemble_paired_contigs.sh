#!/bin/sh

# assemble contigs from non-host reads

#SBATCH --job-name=assemble_contigs_paired
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=30G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load MEGAHIT/1.2.9-GCCcore-9.3.0

set -eu
set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi


for file in "$NONHOST_PAIRED_FASTAS"/paired_*.fasta; do
    filename=${file##/*/}
    output="$OUTPUT_MEGA_PAIRED/${filename/%.fasta/_assembled}"
    megahit -r "$file" -t 16 -o "$output"
    
done

for subdir in "$OUTPUT_MEGA_PAIRED"/*; do
    mv "$subdir"/final.contigs.fa "$subdir".fa; done;
