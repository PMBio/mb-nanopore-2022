#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/hg38.fa

if [ -f ${HG} ]
then
    cd ${BASEDIR}
    bwa mem -t 12 -R"@RG\tID:blood\tSM:blood" ${HG} ${BASEDIR}/../../fastq/illumina/blood.1.fq.gz ${BASEDIR}/../../fastq/illumina/blood.2.fq.gz | samtools sort -o blood.bam -
    samtools index blood.bam
    
    bwa mem -t 12 -R"@RG\tID:tumor\tSM:tumor" ${HG} ${BASEDIR}/../../fastq/illumina/tumor.1.fq.gz ${BASEDIR}/../../fastq/illumina/tumor.2.fq.gz | samtools sort -o tumor.bam -
    samtools index tumor.bam

    bwa mem -t 12 -R"@RG\tID:relapse\tSM:relapse" ${HG} ${BASEDIR}/../../fastq/illumina/relapse01.1.fq.gz ${BASEDIR}/../../fastq/illumina/relapse01.2.fq.gz | samtools sort -o relapse.bam -
    samtools index relapse.bam
fi
