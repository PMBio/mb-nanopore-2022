#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# Joint calling with delly
/opt/dev/delly/bin/delly lr -y ont -g ${HG} -o combined.bcf ${BASEDIR}/../../alignment/ont/*.bam
/opt/dev/delly/src/delly filter -f somatic -o filtered.bcf -s samples.tsv combined.bcf

# Plot insertion
#bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SVTYPE\n" ../tumor.somatic.bcf | grep "BND$"| grep -w -Ff fetch | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$5"Left\n"$3"\t"($4-1000)"\t"($4+1000)"\t"$5"Right";}' > regions.bed
#bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ../tumor.somatic.bcf | grep "DUP$"| grep -w -Ff fetch | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$4"Left\n"$1"\t"($3-1000)"\t"($3+1000)"\t"$4"Right";}' > regions.bed
#/opt/dev/wally/bin/wally region -p -g ${HG} -b /opt/dev/wally/bed/gencode.hg38.bed.gz -R regions.bed -s 2 -u -c -x 2048 -y 2048 ${BASEDIR}/../../alignment/illumina/*.bam
exit;


for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then

	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	# Nanovar
	#nanovar -x ont -t 8 -f hg38 ${BAM} ${HG} ${ID}_dir
	# Sniffles
	#samtools calmd -b ${BAM} ${HG} > tmp.bam 2> /dev/null
	#samtools index tmp.bam
	#sniffles -m tmp.bam -v ${ID}.vcf
	#bgzip ${ID}.vcf
	#rm tmp.bam*
    fi
done
