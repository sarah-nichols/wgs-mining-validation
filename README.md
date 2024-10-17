# WGS Mining for Parasite Sequences

This repository contains scripts and configuration files for identifying parasite sequences from Whole Genome Sequence (WGS) data of the host. Firstly, the host genome sequences (in BAM file format) are filtered for host DNA (using a reference in FASTA format); then there is a step to build a reference database for parasites of interest in Kraken2 (from FASTA files); finally, Kraken2 is used to assign taxonomy to filtered DNA sequences, and Conifer is used to generate confidence scores for assignments. It is also possible to run this analysis by first assembling longer reads from the filtered data and then assigning taxonomies using Kraken2. The scripts for mining WGS sequences are written in bash. The following sections describe the purpose of each script and how to use them.

## Table of Contents

- [Requirements](#requirements)
- [Script Overview and Usage](#script-overview-and-usage)
- [Filter Host](#filter-host)
  - [01-make_reg_exp.sh](#01-make_reg_expsh)
  - [02-build_host_ref.sh](#02-build_host_refsh)
  - [03-filter_bam.sh](#03-filter_bamsh)
- [Build Reference Database](#build-reference-database)
  - [01-sequence_downloads.sh](#01-sequence-downloadssh)
  - [02-build_kraken.sh](#02-build_krakensh)
- [Assign Taxonomy](#assign-taxonomy)
  - [01-assign_tax_paired.sh](#01-assign_tax_pairedsh)
  - [02-assign_confidence_conifer.sh](#02-assign_confidence_conifersh)
- [Assemble contigs and assign taxonomy](#assign-taxonomy)
  - [01-assemble_reads.sh](#01-assemble_readssh)
  - [02-kraken_contigs.sh](#02-kraken_contigssh)
  - [03-assign_confidence_conifer.sh](#03-assign_contig_confidencesh)


## [Requirements](#requirements)

- **SLURM**: A workload manager to submit jobs on a cluster. (Version 20.11.7 or later)
- **GATK**: Genome Analysis Toolkit for processing genome data. (Version 4.1.5.0)
- **SAMtools**: A suite of programs for interacting with high-throughput sequencing data. (Version 1.16.1)
- **SeqKit**: A toolset for working with FASTA/Q files. (Version 2.2.0)
- **Kraken2**: A tool for assigning taxonomic labels to short DNA sequences. (Version 2.1.2)
- **Conifer**: A tool for generating confidence scores for taxonomic assignments. (available here: https://github.com/Ivarz/Conifer.git)
- **Anaconda/Miniconda**: For managing Python environments used in sequence downloads. (Version 4.10.3 or later)
- **TaxonKit**: For processing taxonomy data.
- **wget** 
- **gzip**: For downloading and decompressing data.
- **datasets tool (NCBI)**: For downloading genomic data from NCBI.

Ensure the relevant software is available on your system.

## [Script Overview and Usage](#script-overview-and-usage)

Each script is designed to be submitted via a workload manager such as SLURM. For each script, an email address can be optionally included. Submission either requires submission with input file paths or using a configuration file to specify the necessary parameters. It may also be necessary to define the batch size/number of jobs to be submitted to SLURM for some scripts. The scripts are intended to be run in the order below.

## [Filter Host](#filter-host)

### [01-make_reg_exp.sh](#01-make_reg_expsh)

In order to make the host reference files necessary for the filtering step, we need a FASTA file for the host in a standard format. This script is designed to clean up the host FASTA file by removing commas from the headers.

**Usage:**
Run this script by submitting it to SLURM with the input file, output file, and optionally, an email address for notifications.

```bash
sbatch 01-make_reg_exp.sh <input_file> <output_file> [--mail-user=youremail@mail.com]
```

**Output:**
A FASTA file without commas in the headers to be input in `02-build_host_ref.sh`.

--------

### [02-build_host_ref.sh](#02-build_host_refsh)

The next script uses the host reference in FASTA format to generate the `.hss` and `.img` files that are required for the filtering step. A configuration file is required at submission with `input_reg`, `output_img`, and `output_hss` defined. These specify the pathway to the host FASTA file, corrected to have regular expressions, and the paths to the two output files, respectively.

**Configuration**

This file contains the paths needed for building the host reference.

- `input_reg`: Path to the input FASTA file.
- `output_img`: Path to save the BWA index image.
- `output_hss`: Path to save the PathSeq Kmers file.

**Usage:**
Submit this script to SLURM with the path to the configuration file and optionally, an email address for notifications.

```bash
sbatch 02-build_host_ref.sh <config_file.txt> [--mail-user=youremail@mail.com]
```

**Output:**
Two files (names are specified in the config file): a `.hss` file and an `.img` file to be used in `03-filter_bam.sh`.

-------

### [03-filter_bam.sh](#03-filter_bamsh)

This script filters host BAM files using the `.hss` and `.img` files generated in the previous step. It requires a configuration file that specifies the paths to the sample IDs file, host reference files, and output directory.

**Configuration**

This file contains the paths needed for filtering BAM files.

- `sample_ids_file`: Path to the file containing the list of sample IDs.
- `host_hss`: Path to the PathSeq Kmers file.
- `host_img`: Path to the BWA index image.
- `output_directory`: Directory to save the filtered BAM files and metrics.
- 
**Usage:**
Submit this script to SLURM with the path to the configuration file and optionally, an email address for notifications. N = Total number of samples to process, M = batch size. 

```bash
sbatch --array=1-N%M 03-filter_bam.sh 03-filter_bam_config.txt --mail-user=youremail@mail.com

```

**Output:**
Filtered BAM files and associated metrics, which will be saved in the specified output directory.

-------
## [Build Reference Database](#build-reference-database)


### [01-sequence_downloads.sh](#01-sequence-downloadssh)

The sequence download pipeline facilitates the retrieval of genomic sequences from multiple databases: 
1. **EuPathDB**: Includes various parasite datasets such as PlasmoDB and GiardiaDB. 
2. **NCBI Genomes**: Downloads genomic data for taxonomic groups like Apicomplexa and Euglenozoa. 
3. **PR2 Database**: Provides ribosomal RNA sequences with filtering based on taxonomy. 
4. **NCBI Entrez**: Downloads sequences on the NCBI nucleotide database using taxonomic identifiers. 
5. **WormBase ParaSite**: Downloads genome data from parasitic worm database.

These sources provide a comprehensive dataset for building a custom reference database using Kraken2. Just a warning, this is a bit of a work in progress atm!! I have built the whole database in my directory though so it might make sense just to use the file path to that for now. 

**Configuration**

The primary configuration file is `01-download_all_config.txt`, which specifies paths and settings required for downloading sequences: 
- `download_dir`: Directory where all downloaded sequences will be stored.
- `ncbi_datasets_conda`: Path to the conda environment for downloading NCBI datasets.
- `taxonkit_conda`: Path to the conda environment for TaxonKit (used in NCBI Entrez downloads).
- `entrez_conda`: Path to the conda environment for Entrez-based sequence downloads.
- `taxid_file`: Path to the file containing taxonomic IDs to include in the NCBI Entrez download.
- `taxdump_dir`: Directory for storing NCBI taxonomy data (used in NCBI Entrez downloads).
- `entrez_py_script`: Path to the Python script used for downloading sequences from NCBI Entrez. 

An example configuration file is included in this directory. 

**Usage **

Edit `01-download_all_config.txt` to specify the paths and parameters for your system. Run the main script to submit all download jobs: 

``` sbatch 01-download_all.sh 01-download_all_config.txt```

**Output**

A directory of downloaded sequences ready to be incorporated into a kraken2 database and a .txt file with the pathnames to the downloaded sequences listed. **Warning: this bit is a work in progress!**

### Scripts Overview 
#### `01-download_all.sh` 

The main script that submits jobs for downloading sequences from various sources via SLURM. It reads the configuration file and runs individual download scripts for each data source. You can decide which datasets you want to include in the final database by adding/ removing them here. 

 1. `download_eupath.sh` Downloads and processes datasets from EuPathDB, including sources like PlasmoDB, CryptoDB, and ToxoDB. The script concatenates all sequences into a single FASTA file with cleaned headers. 

2. `download_ncbi_genomes.sh` Downloads genome data for specified taxonomic groups (e.g., Apicomplexa, Euglenozoa) from NCBI. The script consolidates genome data into categorized FASTA files. 

3. `download_pr2.sh` Downloads the PR2 database, filters sequences based on taxonomy to exclude certain groups (e.g. fungus and uncultured), and outputs a cleaned FASTA file. 

4. `download_ncbi_entrez.sh` Uses a Python script to download sequences from NCBI based on taxonomic IDs. It first generates the list of taxonomic IDs using TaxonKit and then downloads sequences using Entrez. 

5. `download_wormbase.sh` Downloads genome data from WormBase ParaSite and processes the files by decompressing, modifying headers, and combining them into a single FASTA file. 


### [02-build_kraken.sh](#02-build_krakensh)

This script builds a custom Kraken2 database by adding sequences to it and then constructing the database. The script requires a configuration file to specify the sequence files and output paths.

**Configuration**

This file contains the paths required for building a custom Kraken2 database.

- `seqs_to_add`: Path to a file listing sequences to add to the Kraken2 library.
- `kraken_custom`: Directory to save the custom Kraken2 database.

**Usage:**
Submit this script to SLURM with the path to the configuration file and optionally, an email address for notifications. This takes about 72 hours to run.

```bash
sbatch 02-build_kraken.sh <config_file.txt> [--mail-user=youremail@mail.com]
```

**Output:**
A custom Kraken2 database built from the specified sequences, ready for use in the taxonomy assignment step.

---------

## Assign Taxonomy

### [01-assign_tax_paired.sh](#01-assign_tax_pairedsh)

This script processes BAM files, converts them to FASTQ format, and classifies sequences using Kraken2. The SLURM array settings should be provided at the time of submission. The script requires a configuration file to specify the input paths and output directory.

**Configuration**

This file contains the paths required for assigning taxonomy using Kraken2.

- `bam_path`: Path to the directory containing BAM files to be processed.
- `kraken_custom`: Path to the custom Kraken2 database directory.
- `output_dir`: Directory to save the Kraken2 output files. 

**Usage:**
Submit this script to SLURM with the appropriate array settings, the path to the configuration file, and optionally, an email address for notifications. N = Total number of samples to process, M = batch size. 

```bash
sbatch --array=1-N%M 01-assign_tax_paired.sh 01-assign_tax_paired_config.txt --mail-user=youremail@mail.com

```
**Output:**
Taxonomic classification reports and associated statistics, saved in the specified output directory.

-------

### [02-assign_confidence_conifer.sh](#02-assign_confidence_conifersh)

This script runs Conifer on Kraken2 output files to generate confidence scores and a final report. The script requires a configuration file to specify the Kraken2 results, the custom taxid database, and the output directory.

**Configuration**
This file contains the paths required for running Conifer to assign confidence scores.

- `kraken_results`: Path to the directory containing Kraken2 output reports.
- `kraken_custom_taxid`: Path to the Kraken2 custom taxid database.
- `output_directory`: Directory to save the Conifer output files.


**Usage:**
Submit this script to SLURM with the path to the configuration file and optionally, an email address for notifications.

```
sbatch 02-assign_confidence_conifer.sh <config_file.txt> [--mail-user=youremail@mail.com]
```

**Output:**

A final report with confidence scores for the taxonomic assignments, saved in the specified output directory.

-------

## [Assemble contigs and assign taxonomy](#assign-taxonomy)

### [01-assemble_reads.sh](#01-assemble_readssh)

This script converts BAM files to FASTQ format and then assembles contigs using metaSPAdes. It also aligns reads back to the assembled contigs and extracts unmapped reads for further processing. The final output combines assembled contigs with unmapped reads, which are saved as a FASTA file.

#### Configuration

The configuration file (`01-assemble_reads_config.txt`) specifies the following:

-   `bam_path`: Directory containing the filtered BAM files.
-   `output_dir`: Directory where all output files will be stored.

#### Usage

Submit the script to SLURM, specifying the configuration file and optionally, an email address for job notifications. The script can be run as a batch job using SLURM's array settings to parallelize processing of multiple BAM files. N = Total number of samples to process, M = batch size.

`sbatch --array=1-N%M 01-assemble_reads.sh 01-assemble_reads_config.txt --mail-user=youremail@mail.com` 

#### Output

The script generates the following output:

-   Paired-end FASTQ files converted from BAM.
-   Contigs assembled with metaSPAdes.
-   Unassembled reads in both FASTQ and FASTA formats.
-   A combined FASTA file containing both assembled contigs and unassembled reads.

----------

### [02-kraken_contigs.sh](#02-kraken_contigssh)

This script assigns taxonomy to contigs using Kraken2. It requires a configuration file to specify the paths to the contig files and the custom Kraken2 database.

#### Configuration

The configuration file (`02-kraken_contigs_config.txt`) specifies the following:

-   `contig_files`: Path to a file listing the combined sequences files to be processed.
-   `kraken_custom`: Path to the custom Kraken2 database directory.
-   `output_dir`: Directory where Kraken2 results will be saved.

#### Usage

Submit the script to SLURM, specifying the configuration file and optionally, an email address for job notifications. The script can be run as a batch job using SLURM's array settings to parallelize processing of multiple contig files.

`sbatch --array=1-N%M 02-kraken_contigs.sh 02-kraken_contigs_config.txt --mail-user=youremail@mail.com` 

#### Output

Kraken2 reports and classification results, saved in the specified output directory.

----------

### [03-assign_confidence_conifer.sh](#03-assign_contig_confidencesh)

This script processes Kraken2 results to generate confidence scores using Conifer, producing a final report.

#### Configuration

The configuration file (`03-assign_confidence_conifer_config.txt`) specifies:

-   `kraken_results`: Path to the directory containing Kraken2 output files.
-   `kraken_custom_taxid`: Path to the Kraken2 custom taxid database.
-   `output_directory`: Directory for storing Conifer's output files.

#### Usage

Submit the script to SLURM with the path to the configuration file and optionally, an email address for notifications.

`sbatch 03-assign_confidence_conifer.sh 03-assign_confidence_conifer_config.txt --mail-user=youremail@mail.com` 

#### Output

Conifer generates a final report containing confidence scores for each taxonomic assignment, which is saved in the specified output directory.

___
Now download the conifer output files and start analysing in R! 
