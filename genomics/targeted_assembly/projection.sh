#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly
HG=${BASEDIR}/../genome/hg38.fa
      
if [ $# -ne 1 ]
then
    echo ""
    echo "Usage: $0 <assembly.fasta>"
    echo ""
    exit -1
fi

INFA=${1}
if [ -f ${INFA} ]
then
    minimap2 -ax map-ont ${HG} ${INFA} | samtools sort -o assembly.bam -
    samtools index assembly.bam
    alfred bam2match -r ${HG} assembly.bam
    #rm assembly.bam assembly.bam.bai
fi
