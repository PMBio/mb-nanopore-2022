#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate phasing

# Read-based phasing statistics
whatshap stats --sample blood --chromosome chr3 --tsv stats.phasing.tsv whatshap.vcf.gz

# ONT, primary tumor
samtools depth -Q 10 ${BASEDIR}/split_final/Primary.h1.bam split_final/Primary.h2.bam | awk '$3>=1 || $4>=1' | wc -l

# ONT, relapse
samtools depth -Q 10 ${BASEDIR}/split_final/Relapse.h1.bam split_final/Relapse.h2.bam | awk '$3>=1 || $4>=1' | wc -l

# Non-N masked, total=2933638686
faCount ${BASEDIR}/../genome/hg38.fa
