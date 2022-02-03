#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then

	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	# Delly
	delly lr -y ont -p 20 -g ${HG} -o ${ID}.bcf ${BAM} > ${ID}.log 2> ${ID}.err &
	# Nanovar
	nanovar -x ont -t 8 -f hg38 ${BAM} ${HG} ${ID}_dir
	# Sniffles
	samtools calmd -b ${BAM} ${HG} > tmp.bam 2> /dev/null
	samtools index tmp.bam
	sniffles -m tmp.bam -v ${ID}.vcf
	bgzip ${ID}.vcf
	rm tmp.bam*
    fi
done
