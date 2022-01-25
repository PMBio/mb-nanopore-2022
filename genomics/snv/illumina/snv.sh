#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# FreeBayes
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities ${BASEDIR}/../../alignment/illumina/*.bam -v freebayes.vcf
bgzip freebayes.vcf
tabix freebayes.vcf.gz
