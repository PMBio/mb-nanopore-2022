#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate phasing

HG=${BASEDIR}/../genome/hg38.fa

# Input snp set
if [ ! -f input.snps.vcf.gz ]
then
    # Strelka & FreeBayes consensus
    bcftools view -s blood -m 2 -M 2 --min-ac 1 --max-ac 1 ${BASEDIR}/../snv/illumina/freebayes.vcf.gz | grep -v -P "\t\./1:" | bgzip > freebayes.vcf.gz
    tabix freebayes.vcf.gz
    bcftools view -s blood -m 2 -M 2 --min-ac 1 --max-ac 1 ${BASEDIR}/../snv/illumina/strelka.vcf.gz | grep -v -P "\t\./1:" | grep -v -P "\t1:" | bgzip > strelka.vcf.gz
    tabix strelka.vcf.gz
    zcat strelka.vcf.gz | grep -v "^#" | cut -f 1-3 > strelka.bed
    zcat freebayes.vcf.gz | grep -v "^#" | cut -f 1-3 > freebayes.bed
    sort -k1,1V -k2,2n freebayes.bed strelka.bed | uniq -u > removed.bed
    bcftools view freebayes.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | grep -v -w -Ff removed.bed | bcftools norm -f ${HG} -m -both - | bgzip > input.snps.vcf.gz
    bcftools annotate -x INFO,^FORMAT/GT input.snps.vcf.gz | bgzip > in.vcf.gz
    mv in.vcf.gz input.snps.vcf.gz
    tabix input.snps.vcf.gz
    rm strelka.vcf.gz* freebayes.vcf.gz*
    rm removed.bed freebayes.bed strelka.bed
fi

# Long-read phasing
if [ ! -f whatshap.vcf.gz ]
then
    whatshap phase --indels --ignore-read-groups --reference ${HG} input.snps.vcf.gz ${BASEDIR}/../alignment/ont/*.bam -o whatshap.vcf
    bgzip whatshap.vcf 
    tabix whatshap.vcf.gz
    ./addAFs.sh whatshap.vcf.gz
    mv merged.bcf whatshap.af.bcf
    mv merged.bcf.csi whatshap.af.bcf.csi
fi

# Statistical phasing
if [ ! -f shapeit.af.bcf ]
then
    FILES=""
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
    do
	if [ -f panel/chr${CHR}.bcf ]
	then
	    echo "Prepare chr${CHR}"
	    bcftools view -O b -o chr${CHR}.in.bcf whatshap.vcf.gz chr${CHR}
	    bcftools index chr${CHR}.in.bcf

	    echo "Shapeit4 phasing chr${CHR}"
	    ${BASEDIR}/shapeit/shapeit4-4.2.2/bin/shapeit4.2 --input chr${CHR}.in.bcf --use-PS 0.0001 --thread 4 --map ${BASEDIR}/shapeit/shapeit4-4.2.2/maps/chr${CHR}.b38.gmap.gz --region chr${CHR} --reference panel/chr${CHR}.bcf --output chr${CHR}.shapeit.bcf
	    bcftools index chr${CHR}.shapeit.bcf
	    rm chr${CHR}.in.bcf chr${CHR}.in.bcf.csi
	    FILES=${FILES}" chr${CHR}.shapeit.bcf"
	fi
    done
    # Concatenate chromosomes
    bcftools concat -O b -o shapeit.bcf ${FILES}
    bcftools index shapeit.bcf
    ./addAFs.sh shapeit.bcf
    mv merged.bcf shapeit.af.bcf
    mv merged.bcf.csi shapeit.af.bcf.csi
    rm chr*.shapeit.bcf*
fi

# Adjust for h1:h2 not equal 1:1
if [ ! -f rephase.af.bcf ]
then
    ./rd.adjust.phase.sh
fi

# Phase remaining variants (non-panel variants) onto scaffold
if [ ! -f final.af.bcf ]
then
    bcftools view -s blood rephase.af.bcf | bcftools annotate -x INFO,^FORMAT/GT - | bgzip > scaffold.vcf.gz
    tabix scaffold.vcf.gz
    whatshap phase --indels --ignore-read-groups --reference ${HG} input.snps.vcf.gz scaffold.vcf.gz ${BASEDIR}/../alignment/ont/Germline.bam -o final.vcf
    bgzip final.vcf
    tabix final.vcf.gz
    bcftools view -O b -o final.bcf final.vcf.gz
    bcftools index final.bcf
    ./addAFs.sh final.bcf
    mv merged.bcf final.af.bcf
    mv merged.bcf.csi final.af.bcf.csi
    rm final.vcf.gz final.vcf.gz.tbi scaffold.vcf.gz scaffold.vcf.gz.tbi
fi

# Split by haplotype
if [ ! -d split_rephase ]
then
    ./split.sh
fi

# Germline statistics
for CS in rephase final
do
    echo ${CS} "biallelic SNPs"
    bcftools view -m2 -M2 -v snps split_${CS}/blood.phased.bcf | grep -v "^#" | wc -l
    echo ${CS} "biallelic InDel"
    bcftools view -m2 -M2 -v indels split_${CS}/blood.phased.bcf | grep -v "^#" | wc -l
done
