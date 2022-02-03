# Genomics analysis code

This repository contains the genomic analysis scripts.

## Installation

To create a conda environment with the required tools just use

`make all`

## Downloading and preparing reference files

The genome files are not included in the repository but can be downloaded using

`cd genome/ && ./prepare_genome.sh`

The input FASTQ files are assumed to be in the `./fastq` subdirectory.

## Alignment

Alignments for the ONT and Illumina data can be carried out using

`cd alignment/ont/ && ./align.sh`

`cd alignment/illumina/ && ./altalign.sh`

## Alignment statistics

To calculate the alignment error rate, genome coverage and other statistics.

`cd qc/ont/ && ./qc.sh`

`cd qc/illumina/ && ./qc.sh`

## Read-depth profiles and copy-number variants

Genome-wide read-depth profiles

`cd coverage/ont/ && ./cov.sh`

`cd coverage/illumina/ && ./cov.sh`

## Single-nucleotide variant (SNVs) calling and small insertions and deletions (InDels)

SNVs and InDel calling using short-read data.

`cd snv/illumina/ && ./snv.sh`

## Structural variant calling

For the short-read data we applied [delly](https://github.com/dellytools/delly) to call structural variants and the pre-defined somatic workflow.

`cd sv/illumina/ && ./delly.sh`

For the long read data we applied multiple methods

`cd sv/ont/ && ./sv.sh`

## Telomere sequences associated with SVs

For the long-read data, [lorax](https://github.com/tobiasrausch/lorax) can be used to call SV to telomere junctions

`cd telomere/ont/ && ./telomere.sh`

For the short-read data, [alfred](https://github.com/tobiasrausch/alfred) computes binned counts of telomeric motifs that can be post-processed to filter tumor-only candidate windows.

`cd telomere/illumina/ && ./telomere.sh`
