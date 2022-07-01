#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# Bi-allelic variants
bcftools view -s blood -m 2 -M 2 --min-ac 1 freebayes.vcf.gz | grep -v -P "\t\./1:" | bcftools norm -f ${HG} -m -both - | bgzip > fb.vcf.gz
tabix fb.vcf.gz
bcftools view -s blood -m 2 -M 2 --min-ac 1 strelka.vcf.gz | grep -v -P "\t\./1:" | bcftools norm -f ${HG} -m -both - | bgzip > st.vcf.gz
tabix st.vcf.gz

# Consensus
rm -rf cons/
bcftools isec fb.vcf.gz st.vcf.gz -p cons -n =2 -w 1
bgzip cons/0000.vcf 
mv cons/0000.vcf.gz cons.vcf.gz
tabix cons.vcf.gz 
rm -rf cons/ fb.vcf.gz* st.vcf.gz*

# Counts
bcftools view -m2 -M2 -v snps cons.vcf.gz | grep -v "^#" | wc -l
bcftools view -m2 -M2 -v indels cons.vcf.gz | grep -v "^#" | wc -l
