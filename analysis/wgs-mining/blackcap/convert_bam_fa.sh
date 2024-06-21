#!/bin/sh

#SBATCH --job-name=convert_bam
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/zool-zost/sann7416/wgs-mining/.env

module purge
module load SAMtools/1.14-GCC-11.2.0
module load seqtk/1.3-foss-2018b

# loop through each interleaved paired bam file in the input directory
for bam_file in "${BAM_INPUT}"/paired_*.bam; do
  # extract the filename without extension
  filename=$(basename -- "${bam_file}")
  filename="${filename%.*}"

  # convert the interleaved paired bam file to separate R1 and R2 fasta files using samtools
  samtools fastq -1 "${FA_OUTPUT}/${filename}_R1.fastq" -2 "${FA_OUTPUT}/${filename}_R2.fastq" -0 /dev/null -s /dev/null -n "${bam_file}"

  # convert the R1 and R2 fastq files to fasta format using seqtk
  seqtk seq -A "${FA_OUTPUT}/${filename}_R1.fastq" > "${FA_OUTPUT}/${filename}_R1.fasta"
  seqtk seq -A "${FA_OUTPUT}/${filename}_R2.fastq" > "${FA_OUTPUT}/${filename}_R2.fasta"
done

