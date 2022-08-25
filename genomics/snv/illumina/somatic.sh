#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# Somatic variants
python filterSomaticSNPs.py -v freebayes.vcf.gz -g blood -t tumor,relapse -q 20 -c 10 -n 1 -o somatic.snv.vcf
bgzip somatic.snv.vcf
tabix somatic.snv.vcf.gz

# Normalize
bcftools norm -O b -o norm.bcf -f ${HG} -m -both somatic.snv.vcf.gz
bcftools index norm.bcf
rm somatic.snv.vcf.gz*

# Filter
bcftools filter -O b -o filtered.bcf -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2 || RPR<=2 || RPL<=2' norm.bcf
bcftools index filtered.bcf
rm norm.bcf norm.bcf.csi

# Subset 
bcftools view filtered.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX | bgzip > somatic.snv.vcf.gz
tabix somatic.snv.vcf.gz
rm filtered.bcf*

# Annotate variants online with VEP
python somatic.py -v somatic.snv.vep.vcf.gz > somatic.tsv
Rscript somatic.R somatic.tsv

# Stats
TOTAL=`tail -n +2 somatic.tsv | wc -l`
SHARED=`cat somatic.tsv  | cut -f 5,6 | awk '$1!=0 && $2!=0' | wc -l`
TUMOR=`cat somatic.tsv  | cut -f 5,6 | awk '$1!=0 && $2==0' | wc -l`
RELAPSE=`cat somatic.tsv  | cut -f 5,6 | awk '$1==0 && $2!=0' | wc -l`
echo ${TOTAL} ${SHARED} ${TUMOR} ${RELAPSE}


