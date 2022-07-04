#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# Joint calling with delly
/opt/dev/delly/bin/delly lr -y ont -g ${HG} -o combined.bcf ${BASEDIR}/../../alignment/ont/*.bam
/opt/dev/delly/src/delly filter -f somatic -o filtered.bcf -s samples.tsv combined.bcf

# By sample with sniffles
for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then

	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	# Sniffles
	samtools calmd -b ${BAM} ${HG} > tmp.bam 2> /dev/null
	samtools index tmp.bam
	sniffles -m tmp.bam -v ${ID}.vcf
	cat ${ID}.vcf | grep -v "STRANDBIAS" | bcftools sort -O b -o ${ID}.bcf -
	bcftools index ${ID}.bcf
	rm ${ID}.vcf tmp.bam*
    fi
done
