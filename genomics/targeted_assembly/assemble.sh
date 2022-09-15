#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/hg38.fa
THREADS=8


# Install Shasta
if [ ! -f shasta-Linux-0.10.0 ]
then
    curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/0.10.0/shasta-Linux-0.10.0
    chmod ugo+x shasta-Linux-0.10.0
fi

if [ -f ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf ]
then
    # Collect reads
    /opt/dev/lorax/src/lorax amplicon -o tumor -g ${HG} -s blood -v ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf -b amplicons.bed ${BASEDIR}/../alignment/ont/Primary.bam
    /opt/dev/lorax/src/lorax amplicon -o relapse -g ${HG} -s blood -v ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf -b amplicons.bed ${BASEDIR}/../alignment/ont/Relapse.bam
    /opt/dev/lorax/src/lorax extract -a -o tumor.match.gz -f tumor.fa.gz -g ${HG} -r tumor.reads ${BASEDIR}/../alignment/ont/Primary.bam
    /opt/dev/lorax/src/lorax extract -a -o relapse.match.gz -f relapse.fa.gz -g ${HG} -r relapse.reads ${BASEDIR}/../alignment/ont/Relapse.bam
    zcat tumor.fa.gz relapse.fa.gz > in.fa
    rm relapse.* tumor.*

    # Assemble
    rm -rf ShastaRun
    ./shasta-Linux-0.10.0 --Reads.minReadLength 3000 --input in.fa --config Nanopore-May2022
    
    # Align to reference
    rm -rf assembly.bam assembly.bam.bai match.gz
    minimap2 -ax map-ont ${HG} ShastaRun/Assembly.fasta | samtools sort -o assembly.bam -
    samtools index assembly.bam
    /opt/dev/alfred/bin/alfred bam2match -r ${HG} assembly.bam
    rm assembly.bam assembly.bam.bai
fi
