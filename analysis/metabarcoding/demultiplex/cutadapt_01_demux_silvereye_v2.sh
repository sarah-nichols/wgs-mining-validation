#!/bin/bash

#SBATCH --job-name=demux
#SBATCH --output=demux
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o demux_%a.out # Standard output
#SBATCH -e demux_%a.err # Standard error
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=e.harney@sheffield.ac.uk

## load profile
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

## Stage 1: Search and trims away the generic illumina library adapters. Sequences are not anchored. 
## Reads pass if either primer is found (not necessary for both to be found)
cutadapt -e 0.0 -O 10 -g AATGATACGGCGACCACCGAGATCTACAC -a ATCTCGTATGCCGTCTTCTGCTTG \
--revcomp -j 32 \
--untrimmed-o no_illumina.fastq.gz \
-o silvereye_illumina_hifi.fastq.gz \
silvereye_hifi.fastq.gz

## Stage 2: Demultiplex based on the i5 and i7 indexes. Sequences are anchored to the start and end of the reads.
## Reads pass if either primer is found (not necessary for both), the Forward adapters are searched first, then the Reverse
## Thus further trimming must be performed afterwards
cutadapt -e 0.0 -O 10 -g file:adapter_sarah_pb_indexF.fasta -a file:adapter_sarah_pb_indexR.fasta \
--revcomp -j 32 \
--untrimmed-o not_demux_silvereye_illumina_hifi.fastq.gz \
-o demux/demux_{name}.fastq.gz \
silvereye_illumina_hifi.fastq.gz

