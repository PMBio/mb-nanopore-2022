#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate visualization

HG=${BASEDIR}/../../genome/hg38.fa

# Plot some random germline deletions and insertions
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" germline.bcf | grep "DEL$" | bgzip > indel.bed.gz
tabix -p bed indel.bed.gz
zcat indel.bed.gz | awk 'NR%100==1 && $2>1000 {print $1"\t"($2-1000)"\t"($3+1000)"\t"$4;}' > regions.bed
# Illumina data
# wally region -pcu -g ${HG} -b indel.bed.gz -R regions.bed -x 2048 ${BASEDIR}/../../alignment/illumina/*.bam
# or ONT data
wally region -pcu -g ${HG} -b indel.bed.gz -R regions.bed -x 2048 ${BASEDIR}/../../alignment/ont/*.bam
rm indel.bed.gz* regions.bed

# ToDo
#bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SVTYPE\n" ../tumor.somatic.bcf | grep "BND$"| grep -w -Ff fetch | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$5"Left\n"$3"\t"($4-1000)"\t"($4+1000)"\t"$5"Right";}' > regions.bed
#bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ../tumor.somatic.bcf | grep "DUP$"| grep -w -Ff fetch | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$4"Left\n"$1"\t"($3-1000)"\t"($3+1000)"\t"$4"Right";}' > regions.bed
#/opt/dev/wally/bin/wally region -p -g ${HG} -b /opt/dev/wally/bed/gencode.hg38.bed.gz -R regions.bed -s 2 -u -c -x 2048 -y 2048 ${BASEDIR}/../../alignment/illumina/*.bam
