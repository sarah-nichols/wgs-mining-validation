#!/bin/bash

#SBATCH --job-name=filter_host
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem=40G
#SBATCH --clusters=all
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL

module purge
module load SAMtools/1.14-GCC-11.2.0
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8


    gatk --java-options "-Xmx80G" PathSeqFilterSpark  \
    --input "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/rarefaction/raw_data/ERS16373033.sorted.bam" \
    --paired-output "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/rarefaction/test/paired_output.sorted.bam" \
    --unpaired-output "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/rarefaction/test/unpaired_output.sorted.bam" \
    --min-adapter-length 1 \
    --min-clipped-read-length 30 \
    --kmer-file /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/host_reference/blackcap_ref_reg.hss \
    --is-host-aligned TRUE \
    --filter-bwa-image /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/wgs-mining/host_reference/blackcap_ref_reg.img \
    --filter-metrics "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/rarefaction/test/filter_metrics.txt" \
    --bam-partition-size 4000000 \
    --verbosity DEBUG