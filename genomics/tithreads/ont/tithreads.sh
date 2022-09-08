#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ${ID} == "Germline" ]; then continue; fi
	echo ${ID}
	/usr/bin/time -v /opt/dev/lorax/src/lorax tithreads -o ${ID} -g ${HG} -m ${BASEDIR}/../../alignment/ont/Germline.bam ${BAM} 
	/usr/bin/time -v /opt/dev/lorax/src/lorax extract -a -g ${HG} -o ${ID}.match.gz -f ${ID}.fa.gz -r ${ID}.reads ${BAM}
	for CLUSTER in `tail -n +2 ${ID}.reads | cut -f 3 | sort | uniq`
	do
	    echo ${ID} ${CLUSTER}
	    zcat ${ID}.fa.gz | grep --no-group-separator -A 1 -w -Ff <(awk '$3=="'${CLUSTER}'"' ${ID}.reads | cut -f 1 | sort | uniq) > input.fasta
	    /usr/bin/time -v ${BASEDIR}/../../assembly/shasta-Linux-0.10.0 --input input.fasta --config Nanopore-May2022
	    cp ShastaRun/Assembly.fasta ${ID}.assembly.${CLUSTER}.fasta
	    rm -rf ShastaRun/ input.fasta
	    if [ -f ${ID}.assembly.${CLUSTER}.fasta ]
	    then
		minimap2 -t 6 -a -x map-ont -L ${HG} ${ID}.assembly.${CLUSTER}.fasta | samtools sort -o ${ID}.assembly.${CLUSTER}.bam
		samtools index ${ID}.assembly.${CLUSTER}.bam

		# Visualize contigs
		#samtools view ${ID}.assembly.${CLUSTER}.bam | cut -f 1 | sort | uniq > reads
		#/opt/dev/wally/src/wally dotplot -a -g ${HG} -R reads ${ID}.assembly.${CLUSTER}.bam
	    fi
	done
    fi
done

# Filter for self-concatenation and high copy-number
echo -e "sample\tinstances\tloci" > stats.tsv
for TYPE in Primary Relapse
do
    INST=`tail -n +2 ${TYPE}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq  | wc -l | cut -f 1`
    NUM=0
    if [ ${INST} -gt 0 ]
    then
	for CLUSTID in `tail -n +2 ${TYPE}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq`
	do
	    SEGMENTS=`cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 7 | awk '{SUM+=$1; BASE+=2;} END {print (SUM-BASE);}'`
	    echo ${SEGMENTS} | sed "s/^/${TYPE}\t${CLUSTID}\t/" >> stats.ti.seglen
	    cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | sed "s/^/${TYPE}\t/" >> stats.ti.tsv
	    NUMCLUST=`cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 1-3 | sort | uniq  | wc -l | cut -f 1`
	    NUM=`expr ${NUM} + ${NUMCLUST}`
	done
    fi		
    echo -e "${TYPE}\t${INST}\t${NUM}" >> stats.tsv
done
