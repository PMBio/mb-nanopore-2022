#!/bin/bash

for CHR in `seq 1 22` X
do
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
    FILE=`ls *chr${CHR}.*.vcf.gz`
    if [ -f ${FILE} ]
    then
	bcftools view --no-version -O b -o chr${CHR}.bcf ${FILE}
	bcftools index chr${CHR}.bcf
    fi
done
