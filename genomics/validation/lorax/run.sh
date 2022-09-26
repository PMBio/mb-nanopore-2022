#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

for BAM in ${BASEDIR}/../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ${ID} == "Germline" ]; then continue; fi
	echo ${ID}
	ls ${BASEDIR}/../../alignment/ont/Germline.bam ${BAM}
	lorax tithreads -o ${ID} -g ${HG} -m ${BASEDIR}/../../alignment/ont/Germline.bam ${BAM}
	lorax extract -a -g ${HG} -o ${ID}.match.gz -f ${ID}.fa.gz -r ${ID}.reads ${BAM}
	for CLUSTER in `tail -n +2 ${ID}.reads | cut -f 3 | sort | uniq`
	do
	    echo ${ID} ${CLUSTER}
	    zcat ${ID}.fa.gz | grep --no-group-separator -A 1 -w -Ff <(awk '$3=="'${CLUSTER}'"' ${ID}.reads | cut -f 1 | sort | uniq) > input.fasta
	    ${BASEDIR}/../../assembly/shasta-Linux-0.10.0 --input input.fasta --config Nanopore-May2022
	    cp ShastaRun/Assembly.fasta ${ID}.assembly.${CLUSTER}.fasta
	    rm -rf ShastaRun/ input.fasta
	    if [ -f ${ID}.assembly.${CLUSTER}.fasta ]
	    then
		minimap2 -t 6 -a -x map-ont -L ${HG} ${ID}.assembly.${CLUSTER}.fasta | samtools sort -o ${ID}.assembly.${CLUSTER}.bam
		samtools index ${ID}.assembly.${CLUSTER}.bam
		
		# Visualize contigs
		#samtools view ${ID}.assembly.${CLUSTER}.bam | cut -f 1 | sort | uniq > reads
		#wally dotplot -a -g ${HG} -R reads ${ID}.assembly.${CLUSTER}.bam
	    fi
	done
    fi
done
