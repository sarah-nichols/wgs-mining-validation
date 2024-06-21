#!/bin/sh

# convert filtered output files to fastq
# visualise quality on fastq


#SBATCH --job-name=quality_check
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load SAMtools/1.14-GCC-11.2.0
module load FastQC/0.11.9-Java-11

set -eu
source /data/zool-zost/sann7416/wgs-mining/.env

for bam in "$OUTPATH"/filter_test/paired_*.sorted.bam;
do OUTPUT="${bam/%.sorted.bam/.fastq}"; INPUT=$bam; samtools fastq $INPUT > $OUTPUT;
fastqc $OUTPUT -o "$OUTPATH"/filter_test/fastqc_output/;
done
