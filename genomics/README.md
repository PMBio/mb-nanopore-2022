# Genomics analysis code

This repository contains the genomic analysis scripts.

# Installation

To create a conda environment with the required tools just use

`make all`

# Downloading and preparing reference files

The genome files are not included in the repository but can be downloaded using

`cd genome/ && ./prepare_genome.sh`

The input FASTQ files are assumed to be in the `./fastq` subdirectory.

# Alignment

Alignments for the ONT and Illumina data can be carried out using

`cd alignment/ont/ && ./align.sh`

`cd alignment/illumina/ && ./altalign.sh`

# Alignment statistics

To calculate the alignment error rate, genome coverage and other statistics.

`cd qc/ont/ && ./qc.sh`
