#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate truvari

# Delly
bcftools view -s Germline --min-ac 1 delly.bcf | bgzip > delly.germ.vcf.gz
tabix delly.germ.vcf.gz

# Sniffles
bcftools view -s Germline --min-ac 1 sniffles.vcf.gz | bgzip > sniffles.germ.vcf.gz
tabix sniffles.germ.vcf.gz

# Comparison of calls
rm -rf germ_stats/
truvari bench -p 0 --gtcom --passonly --no-ref a -r 1000 -C 1000 -b delly.germ.vcf.gz -c sniffles.germ.vcf.gz -f ${BASEDIR}/../../genome/hg38.fa -o germ_stats
rm delly.germ.vcf.gz* sniffles.germ.vcf.gz*

# Fetch consensus variants
bcftools sort germ_stats/tp-base.vcf | grep -v "SVTYPE=INV" | grep -v "SVTYPE=DUP" | grep -v "SVTYPE=BND" | bgzip > germline.vcf.gz
bcftools index germline.vcf.gz
bcftools view -O b -o germline.bcf germline.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
bcftools index germline.bcf
rm germline.vcf.gz*
rm -rf germ_stats/
