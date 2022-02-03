#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# GRCh38 primary assembly used for ONT mapping
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/^>/>chr/' > Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
rm Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# GRCh38 with decoy sequences for Illumina mapping
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
bwa index GRCh38_full_analysis_set_plus_decoy_hla.fa

# Telomere-to-telomere assembly (T2T)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
gunzip chm13.draft_v1.1.fasta.gz
samtools faidx chm13.draft_v1.1.fasta

