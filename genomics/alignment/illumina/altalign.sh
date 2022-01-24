#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

if [ -f ${HG} ]
then
    if [ ! -d ${BASEDIR}/bwa.kit ]
    then
	cd ${BASEDIR}
	wget -O- http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download > bwakit.tar
	tar -xf bwakit.tar
	rm bwakit.tar
    fi
    ${BASEDIR}/bwa.kit/run-bwamem -o blood -t 8 -HR"@RG\tID:blood\tSM:blood" ${HG} ${BASEDIR}/../../fastq/illumina/blood.1.fq.gz ${BASEDIR}/../../fastq/illumina/blood.2.fq.gz | sh
    samtools sort -o blood.bam blood.aln.bam
    samtools index blood.bam
    rm blood.aln.bam
    
    ${BASEDIR}/bwa.kit/run-bwamem -o tumor -t 8 -HR"@RG\tID:tumor\tSM:tumor" ${HG} ${BASEDIR}/../../fastq/illumina/tumor.1.fq.gz ${BASEDIR}/../../fastq/illumina/tumor.2.fq.gz | sh
    samtools sort -o tumor.bam tumor.aln.bam
    samtools index tumor.bam
    rm tumor.aln.bam

    ${BASEDIR}/bwa.kit/run-bwamem -o relapse -t 8 -HR"@RG\tID:relapse\tSM:relapse" ${HG} ${BASEDIR}/../../fastq/illumina/relapse01.1.fq.gz ${BASEDIR}/../../fastq/illumina/relapse01.2.fq.gz | sh
    samtools sort -o relapse.bam relapse.aln.bam
    samtools index relapse.bam
    rm relapse.aln.bam
fi
