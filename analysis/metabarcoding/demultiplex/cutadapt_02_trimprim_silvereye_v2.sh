#!/bin/bash

#SBATCH --job-name=trimprim
#SBATCH --output=trimprim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o logs/trimprim_%a.out # Standard output
#SBATCH -e logs/trimprim_%a.err # Standard error
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=e.harney@sheffield.ac.uk

## load profile
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat file_list.txt))
File_fastq=$(echo "$data" | cut -f1 )

## The final cutadapt is based on the 3NDF2 and UNonMet_DB primers, provided as a mandatory linked primer (but not anchored). 
## Also important to note is that 3NDF2 is in its original sense but UNonMet_DB is provided as a reverse complement.
cutadapt -e 0.1 -O 10 -g ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNGCAAGTCTGGTGCCAG...CGCAAGNNTGAAANTTAAAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--untrimmed-o demux_not_trimmed.fastq.gz -j 8 \
-o trimmed/trim${File_fastq} \
demux/${File_fastq}
