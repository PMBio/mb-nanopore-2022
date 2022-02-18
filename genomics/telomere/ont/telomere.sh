#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate telomere

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
CHM13=${BASEDIR}/../../genome/chm13.draft_v1.1.fasta

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ${ID} == "Germline" ]; then continue; fi
	echo ${ID}
	# Identify telomere associated SVs
	/opt/dev/lorax/src/lorax telomere -o ${ID}.bed.gz -g ${HG} -m ${BASEDIR}/../../alignment/ont/Germline.bam ${BAM}
	zcat ${ID}.bed.gz  | tail -n +2 | cut -f 4 | sort | uniq > ${ID}.reads
	# Extract reads
	/opt/dev/lorax/src/lorax extract -g ${HG} -r ${ID}.reads -o ${ID}.match.gz -f ${ID}.fa.gz ${BAM}

	# Map to CHM13
	minimap2 -ax map-ont ${CHM13} ${ID}.fa.gz | samtools sort -o assembly.bam -
	samtools index assembly.bam
	/opt/dev/lorax/src/lorax extract -g ${CHM13} -r ${ID}.reads -o ${ID}.chm13.match.gz -f ${ID}.chm13.fa.gz assembly.bam
	rm assembly.bam assembly.bam.bai
    fi
done

