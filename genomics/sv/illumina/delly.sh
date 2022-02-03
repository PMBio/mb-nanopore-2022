#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Delly
delly call -g ${HG} -x human.hg38.excl.tsv -o delly.bcf ${BASEDIR}/../../alignment/illumina/*.bam