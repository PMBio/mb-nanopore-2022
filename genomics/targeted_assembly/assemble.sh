#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/hg38.fa
THREADS=8

if [ -f ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf ]
then
    # Collect reads
    /opt/dev/lorax/src/lorax amplicon -o tumor -g ${HG} -s blood -v ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf -b amplicons.bed ${BASEDIR}/../alignment/ont/Primary.bam
    /opt/dev/lorax/src/lorax amplicon -o relapse -g ${HG} -s blood -v ${BASEDIR}/../phasing/split_rephase/blood.phased.bcf -b amplicons.bed ${BASEDIR}/../alignment/ont/Relapse.bam
    /opt/dev/lorax/src/lorax extract -a -o tumor.match.gz -f tumor.fa.gz -g ${HG} -r tumor.reads ${BASEDIR}/../alignment/ont/Primary.bam
    /opt/dev/lorax/src/lorax extract -a -o relapse.match.gz -f relapse.fa.gz -g ${HG} -r relapse.reads ${BASEDIR}/../alignment/ont/Relapse.bam
    zcat tumor.fa.gz relapse.fa.gz > in.fa
    #rm relapse.* tumor.*

    # Assemble
    wtdbg2 -x ont -g 2m -i in.fa -t ${THREADS} -fo dbg
    wtpoa-cns -t ${THREADS} -i dbg.ctg.lay.gz -fo dbg.raw.fa
    minimap2 -t ${THREADS} -ax map-ont -r 2000 dbg.raw.fa in.fa | samtools sort -@ ${THREADS} > dbg.bam
    samtools view -F0x900 dbg.bam | wtpoa-cns -t ${THREADS} -d dbg.raw.fa -i - -fo dbg.cns.fa
    mv dbg.cns.fa assembly.wtdbg2.fa
    rm -f dbg.* in.fa

    # Align to reference
    minimap2 -ax map-ont ${HG} assembly.wtdbg2.fa | samtools sort -o assembly.bam -
    samtools index assembly.bam
    /opt/dev/alfred/bin/alfred bam2match -r ${HG} assembly.bam
    rm assembly.bam assembly.bam.bai
fi
