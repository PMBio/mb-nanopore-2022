#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../../genome/hg38.fa
ID=contig
if [ -f ${HG} ]
then
    # Align to GRCh38
    minimap2 -t 6 -a -x map-ont -L ${HG} shasta_contig_9046.fa | samtools sort -o ${ID}.bam
    samtools index ${ID}.bam

    # Get matches
    alfred bam2match -r ${HG} -o ${ID}.hg38.match.gz ${ID}.bam

    # Get segment counts
    bedtools merge -i <(zcat contig.hg38.match.gz | tail -n +2) | awk '{print $0"\tSegmentID"NR;}' > segments.bed
    bedtools intersect -a segments.bed -b <(zcat contig.hg38.match.gz | tail -n +2) -wao | cut -f 1-4 | sort | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1;}' | awk '$5>=3' > segment.counts.bed
    rm segments.bed

    # Get segment CN estimate
    echo -e "chr\tstart\tend\tid\tcount\till\tont" > segment.counts.cov.bed
    
    while read CHR START END ID COUNT
    do
	echo ${ID}
	if [ ${START} == 7804862 ]
	then
	    # This segment got enlarged because of a larger overlapping segment at low copy-number
	    START=7805260
	    END=7805460
	fi
	ILL=`samtools depth -aa -r ${CHR}:${START}-${END} ../../../alignment/illumina/tumor.bam ../../../alignment/illumina/blood.bam | awk '{T+=$3; N+=$4;} END {print T/N * 2;}'`
	#ILL=`samtools depth -aa -r ${CHR}:${START}-${END} ../../../alignment/illumina/relapse.bam ../../../alignment/illumina/blood.bam | awk '{T+=$3; N+=$4;} END {print T/N * 2;}'`
	ONT=`samtools depth -aa -r ${CHR}:${START}-${END} ../../../alignment/ont/Primary.bam ../../../alignment/ont/Germline.bam | awk '{T+=$3; N+=$4;} END {print T/N * 2;}'`
	#ONT=`samtools depth -aa -r ${CHR}:${START}-${END} ../../../alignment/ont/Relapse.bam ../../../alignment/ont/Germline.bam | awk '{T+=$3; N+=$4;} END {print T/N * 2;}'`
	echo -e "${CHR}\t${START}\t${END}\t${ID}\t${COUNT}\t${ILL}\t${ONT}" >> segment.counts.cov.bed
    done < segment.counts.bed
    rm segment.counts.bed
fi

# Plot
Rscript  occ_count.R
