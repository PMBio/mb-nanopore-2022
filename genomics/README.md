# Genome analysis code

This repository contains the genome analysis scripts. Associated tools are available in the [lorax](https://github.com/tobiasrausch/lorax) and [rayas](https://github.com/tobiasrausch/rayas) GitHub repositories.

## Installation

To create a conda environment with the required tools just use

`git clone https://github.com/PMBio/mb-nanopore-2022.git`

`cd mb-nanopore-2022/genomics/`

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

## Phasing of single-nucleotide variants

Assuming low to moderate long-read coverage, haplotyping is implemented as a two-step process using (1) [whatshap](https://whatshap.readthedocs.io/) to compute initial haplotype blocks based on long-reads and (2) [shapeit](https://odelaneau.github.io/shapeit4/) to scaffold haplotype blocks using the 1000 Genomes reference panel. For regions in the tumor genome that deviate from the expected 1:1 haplotype ratio we implemented a somatic copy-number aware switch-error correction algorithm to correct these switch errors. The pipeline can be run using

`cd phasing/ && ./phase.sh`

and requires the SNV calls, the genetic maps, the reference panel and the long-read alignments.

## Structural variant calling

For the short-read data we applied [delly](https://github.com/dellytools/delly) to call structural variants and the pre-defined somatic workflow.

`cd sv/illumina/ && ./delly.sh`

For the long read data we applied multiple methods

`cd sv/ont/ && ./sv.sh`

## Figures and plots

Read-depth plots with phased heterozygous allele frequencies and structural variants were generated using

`cd plotting/baf_rd_sv/ && ./rd.sh`

Alignment plots for SVs were generated with [IGV](https://software.broadinstitute.org/software/igv/) and [wally](https://github.com/tobiasrausch/wally) using

`cd plotting/sv_alignments/ && ./wally.sh`

## Telomere sequences associated with SVs

For the long-read data, [lorax](https://github.com/tobiasrausch/lorax) can be used to call SV to telomere junctions

`cd telomere/ont/ && ./telomere.sh`

For the short-read data, [alfred](https://github.com/tobiasrausch/alfred) computes binned counts of telomeric motifs that can be post-processed to filter tumor-only candidate windows.

`cd telomere/illumina/ && ./telomere.sh`

## Templated insertion threads

For the long-read data, [lorax](https://github.com/tobiasrausch/lorax) can be used to call templated insertion threads

`cd tithreads/ont/ && ./tithreads.sh`

For the short-read data, [rayas](https://github.com/tobiasrausch/rayas) infers templated insertion threads using split-reads and depth of coverage

`cd tithreads/illumina/ && ./tithreads.sh`

To compute the overlap between long-read and short-read predictions

`cd tithreads/analysis/overlap_ont_illumina/ && ./overlap.sh`

Various plotting scripts for self-alignments and alignments against templated insertion source sequences for selected long reads

`cd tithreads/analysis/plotting/ && ./run.sh`

Enrichment/depletion analysis with respect to annotated repeats

`cd tithreads/analysis/repeats/ && ./rmsk.sh`

Internal repeat motif length using subsampling

`cd tithreads/analysis/repetitions/ && ./subsample.sh`

PCAWG related analysis scripts are available in a separate folder

`cd tithreads/analysis/pcawg/ && ./run.sh`

## Targeted assembly

Using phased germline variants and amplicon regions from a somatic copy-number alteration analysis (SCNAs) [lorax](https://github.com/tobiasrausch/lorax) can be used to select reads for a targeted amplicon assembly. These reads can then be used with any long-read assembler, our choice in this project was [wtdbg2](https://github.com/ruanjue/wtdbg2). Lastly, we project the alignments back to the reference to infer breakpoints.

`cd targeted_assembly/ && ./assemble.sh`

