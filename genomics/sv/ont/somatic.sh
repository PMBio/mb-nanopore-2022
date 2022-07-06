#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate truvari

# Delly, somatic filtering (tumor and relapse)
delly filter -f somatic -o somatic.delly.bcf -s samples.tsv delly.bcf
bcftools view -O b -o tmp.bcf somatic.delly.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
mv tmp.bcf somatic.delly.bcf
bcftools index -f somatic.delly.bcf

# Sniffles, somatic filtering (tumor and relapse)
bcftools query -f "%ID\n" sniffles.vcf.gz | sort | uniq > all.ids
bcftools query -f "[%GT\t%DV\t]%ID\n" sniffles.vcf.gz  | grep -P "^0/0\t0"  | sed 's/^.*Sniffles/Sniffles/' | sort | uniq > selected.ids
sort selected.ids all.ids | uniq -u > remove.ids
bcftools view sniffles.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | grep -v -w -Ff remove.ids | bcftools view -O b -o somatic.sniffles.bcf -
bcftools index somatic.sniffles.bcf
rm all.ids selected.ids remove.ids

# Primary_tumor and relapse
for SAMPLE in Primary_tumor Relapse
do
    echo ${SAMPLE}

    # Delly
    bcftools view -s ${SAMPLE} somatic.delly.bcf | bcftools query -f "%ID[\t%RV]\n" -  | awk '$2<2'  | cut -f 1 | sort | uniq > remove.ids
    bcftools view -s ${SAMPLE} somatic.delly.bcf | grep -v -w -Ff remove.ids > somatic.delly.vcf
    sniffles --input ${BASEDIR}/../../alignment/ont/Primary_tumor.bam --genotype-vcf somatic.delly.vcf --vcf genotypes.tumor.vcf
    sniffles --input ${BASEDIR}/../../alignment/ont/Germline.bam --genotype-vcf somatic.delly.vcf --vcf genotypes.germline.vcf
    bcftools view genotypes.tumor.vcf | bcftools query -f "%ID[\t%DV]\n" - | awk '$2<2'  | cut -f 1 | sort | uniq > remove.ids
    bcftools view genotypes.germline.vcf | bcftools query -f "%ID[\t%DV]\n" - | awk '$2>0'  | cut -f 1 | sort | uniq >> remove.ids
    sort remove.ids | uniq > remove.ids.tmp
    mv remove.ids.tmp remove.ids
    cat somatic.delly.vcf | grep -v -w -Ff remove.ids | bgzip > ${SAMPLE}.vcf.gz
    tabix ${SAMPLE}.vcf.gz
    rm remove.ids genotypes.tumor.vcf genotypes.germline.vcf somatic.delly.vcf 

    # Clean-up
    #rm somatic.delly.vcf.gz* somatic.sniffles.vcf.gz*
done
rm somatic.delly.bcf* somatic.sniffles.bcf*
