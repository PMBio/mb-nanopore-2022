#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate visualization

HG=${BASEDIR}/../../genome/hg38.fa

# Plot somatic SVs
for SAMPLE in Primary_tumor Relapse
do
    mkdir -p ${SAMPLE}
    bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ${SAMPLE}.vcf.gz | egrep "DEL$|INS$|INV$|DUP$" > sv.bed
    bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SVTYPE\n" ${SAMPLE}.vcf.gz | grep "BND$" >> sv.bed
    sort -k1,1V -k2,2n sv.bed | cut -f 1-5 | awk '$2>1000 && $4>1000 {print $1"\t"($2-1000)"\t"($2+1000)"\t"$5"Left\n"$3"\t"($4-1000)"\t"($4+1000)"\t"$5"Right";}' > regions.bed
    rm sv.bed
    # Illumina data
    #/opt/dev/wally/bin/wally region -puc -g ${HG} -b /opt/dev/wally/bed/gencode.hg38.bed.gz -R regions.bed -s 2 -x 2048 ${BASEDIR}/../../alignment/illumina/*.bam
    # or ONT data
    /opt/dev/wally/bin/wally region -puc -g ${HG} -b /opt/dev/wally/bed/gencode.hg38.bed.gz -R regions.bed -s 2 -x 2048 ${BASEDIR}/../../alignment/ont/*.bam
    rm regions.bed
    mv *.png ${SAMPLE}
done
exit;

# Plot some random germline deletions and insertions
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" germline.bcf | grep "DEL$" | bgzip > indel.bed.gz
tabix -p bed indel.bed.gz
zcat indel.bed.gz | awk 'NR%100==1 && $2>1000 {print $1"\t"($2-1000)"\t"($3+1000)"\t"$4;}' > regions.bed
# Illumina data
# wally region -pcu -g ${HG} -b indel.bed.gz -R regions.bed -x 2048 ${BASEDIR}/../../alignment/illumina/*.bam
# or ONT data
wally region -pcu -g ${HG} -b indel.bed.gz -R regions.bed -x 2048 ${BASEDIR}/../../alignment/ont/*.bam
rm indel.bed.gz* regions.bed

