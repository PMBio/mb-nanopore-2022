#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
THREADS=8
INFA=${1}
if [ -f ${INFA} ]
then
    minimap2 -ax map-ont -t ${THREADS} ${HG} ${INFA} | samtools sort -o assembly.bam -
    samtools index assembly.bam
    alfred bam2match -r ${HG} assembly.bam
    rm assembly.bam assembly.bam.bai
fi
