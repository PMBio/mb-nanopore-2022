#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate truvari

# Delly, somatic filtering
delly filter -f somatic -o somatic.delly.bcf -s samples.tsv delly.bcf

# Sniffles, require absence in Germline
bcftools query -f "%ID\n" sniffles.vcf.gz | sort | uniq > all.ids
bcftools query -f "[%GT\t%DV\t]%ID\n" sniffles.vcf.gz  | grep -P "^0/0\t0"  | sed 's/^.*Sniffles/Sniffles/' | sort | uniq > selected.ids
sort selected.ids all.ids | uniq -u > remove.ids
bcftools view sniffles.vcf.gz | grep -v -w -Ff remove.ids | bcftools view -O b -o somatic.sniffles.bcf -
bcftools index somatic.sniffles.bcf
rm all.ids selected.ids remove.ids

# Comparison
for SAMPLE in Primary_tumor Relapse
do
    echo ${SAMPLE}
    bcftools view -s ${SAMPLE} somatic.delly.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | bgzip > somatic.delly.vcf.gz
    tabix somatic.delly.vcf.gz
    bcftools view -s ${SAMPLE} somatic.sniffles.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | bgzip > somatic.sniffles.vcf.gz
    tabix somatic.sniffles.vcf.gz
    rm -rf soma_stats/
    truvari bench -s 0 -S 0 --sizemax 250000000 --multimatch -t -p 0 -b somatic.delly.vcf.gz -c somatic.sniffles.vcf.gz -f ${BASEDIR}/../../genome/hg38.fa -o soma_stats
    bcftools sort soma_stats/tp-base.vcf | bgzip > ${SAMPLE}.vcf.gz
    bcftools index ${SAMPLE}.vcf.gz
    rm -rf soma_stats/
done
rm somatic.delly* somatic.sniffles*
