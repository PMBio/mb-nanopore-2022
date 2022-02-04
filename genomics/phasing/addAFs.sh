#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate phasing

HG=${BASEDIR}/../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

ISFILE=${BASEDIR}/../snv/illumina/freebayes.vcf.gz
PHASEFILE=${1}

if [ -f ${ISFILE} ]
then
    if [ -f ${PHASEFILE} ]
    then
	bcftools view -p -s blood ${PHASEFILE} | bcftools annotate -x INFO,^FORMAT/GT | bgzip > partB.vcf.gz
	tabix partB.vcf.gz
	bcftools view -s tumor,relapse ${ISFILE} | grep "^#" | bgzip > partA.vcf.gz
	bcftools view -s tumor,relapse ${ISFILE} | grep -v "^#" | grep -w -Ff <(zcat partB.vcf.gz | grep -v "^#" | cut -f 1-3) | bgzip >> partA.vcf.gz
	tabix partA.vcf.gz
	bcftools merge -O b -o merged.bcf partA.vcf.gz partB.vcf.gz
	bcftools index merged.bcf
	rm partA.vcf.gz* partB.vcf.gz*
    fi
fi
