#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate truvari

HG=${BASEDIR}/../../genome/hg38.fa

# Joint calling with delly
delly lr -y ont -g ${HG} -o delly.bcf ${BASEDIR}/../../alignment/ont/*.bam
delly filter -f somatic -o somatic.delly.bcf -s samples.tsv delly.bcf

# By sample with sniffles
for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then

	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	# Sniffles
	sniffles --reference ${HG} --input ${BAM} --snf ${ID}.snf --reference ${HG}
    fi
done

# Combine
sniffles --input Germline.snf Primary_tumor.snf Relapse.snf --reference ${HG} --vcf sniffles.vcf
bgzip sniffles.vcf
tabix sniffles.vcf.gz
rm Germline.snf Primary_tumor.snf Relapse.snf
