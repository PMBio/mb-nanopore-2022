#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# FreeBayes
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities ${BASEDIR}/../../alignment/illumina/*.bam -v freebayes.vcf
bgzip freebayes.vcf
tabix freebayes.vcf.gz

# Strelka
cd ${BASEDIR}
if [ ! -d strelka-2.9.2.centos6_x86_64 ]
then
    wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
    tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
    rm strelka-2.9.2.centos6_x86_64.tar.bz2 
fi

# Run strelka
./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam ${BASEDIR}/../../alignment/illumina/blood.bam  --bam ${BASEDIR}/../../alignment/illumina/relapse.bam  --bam ${BASEDIR}/../../alignment/illumina/tumor.bam --referenceFasta ${HG} --runDir strelka_germline
python strelka_germline/runWorkflow.py -m local -j 4
mv strelka_germline/results/variants/variants.vcf.gz strelka.vcf.gz
mv strelka_germline/results/variants/variants.vcf.gz.tbi strelka.vcf.gz.tbi
