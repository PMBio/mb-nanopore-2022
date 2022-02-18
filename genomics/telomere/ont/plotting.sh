#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate telomere

# Plotting (Primary tumor)
Rscript match.R Primary_tumor.match.gz 11df23a3-af83-4770-ba05-daecde6a2481 "Tumor, chr5p telomere"
/opt/dev/wally/src/wally region -cu -x 512 -y 512 -g ../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -r chr7:7802367-7804367 -b /opt/dev/wally/bed/gencode.hg38.bed.gz ../../alignment/ont/Primary_tumor.bam

Rscript match.R Primary_tumor.match.gz f47d2994-7116-4741-abd7-ec9527bf5615 "Tumor, chr22q telomere"
/opt/dev/wally/src/wally region -cu -x 512 -y 512 -g ../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -r chr2:43463192-43465192 -b /opt/dev/wally/bed/gencode.hg38.bed.gz ../../alignment/ont/Primary_tumor.bam

Rscript match.R Primary_tumor.chm13.match.gz 501ef0df-09a4-4769-b9f0-941d7dab7f97 "Tumor, chr16q telomere (T2T)"
/opt/dev/wally/src/wally region -cu -x 512 -y 512 -g ../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -r chr19:8056491-8058491 -b /opt/dev/wally/bed/gencode.hg38.bed.gz ../../alignment/ont/Primary_tumor.bam

# Plotting (Relapse)
Rscript match.R Relapse.match.gz 4df1499d-79b5-42f5-b1bf-ef2263484d2f "Relapse, chr5p telomere"
/opt/dev/wally/src/wally region -cu -x 512 -y 512 -g ../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -r chr14:25169796-25171796  -b /opt/dev/wally/bed/gencode.hg38.bed.gz ../../alignment/ont/Relapse.bam

Rscript match.R Relapse.chm13.match.gz df5064a5-e792-4d94-bdee-f987e7d32a2e "Relapse, chr7p telomere (T2T)"
/opt/dev/wally/src/wally region -cu -x 512 -y 512 -g ../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -r chr11:24853194-24855194  -b /opt/dev/wally/bed/gencode.hg38.bed.gz ../../alignment/ont/Relapse.bam
