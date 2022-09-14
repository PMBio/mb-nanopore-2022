#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate variants

# ONT, primary tumor
samtools depth -Q 10 ont/Primary.bam | awk '$3>=5' | wc -l

# Illumina, primary tumor
samtools depth -Q 10 illumina/tumor.bam | awk '$3>=5' | wc -l

# Non-N masked
faCount ../genome/hg38.fa
