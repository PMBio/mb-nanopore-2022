#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# GRCh38
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz > hg38.fa
samtools faidx hg38.fa
rm GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz

# Telomere-to-telomere assembly (T2T)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
zcat chm13v2.0_maskedY_rCRS.fa.gz > t2t.fa
samtools faidx t2t.fa
rm chm13v2.0_maskedY_rCRS.fa.gz

# Circle-seq reference genome
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip -c hs37d5.fa.gz > hs37d5_mm10.fa
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.p6.genome.fa.gz
gunzip -c GRCm38.p6.genome.fa.gz > mm10.fa &&  sed -i -e 's/>/>m\./g' mm10.fa && cat mm10.fa >> hs37d5_mm10.fa
rm hs37d5.fa.gz GRCm38.p6.genome.fa.gz mm10.fa

# Index
bwa index hg38.fa
bwa index t2t.fa
bwa index hs37d5_mm10.fa

