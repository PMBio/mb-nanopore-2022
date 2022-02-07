#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Execute rayas for each tumor-normal pair (Embassy cloud)
#rayas call -g ${HG} -o rayas/<project>/<id>.bed -m <control.bam> <tumor.bam> > rayas/<project>/<id>.log

# Summarize templated insertion threads in PCAWG
echo -e "tumour_specimen_aliquot_id\tproject\taffected\tinstances\tloci" > stats.tsv
head -n 1 ${BASEDIR}/../../illumina/tumor.bed | sed 's/^/tumour_specimen_aliquot_id\t/' > stats.ti.tsv
echo -e "tumour_specimen_aliquot_id\tclusterid\tseglen" > stats.ti.seglen
for F in ${BASEDIR}/rayas/*/*.log
do
    ID=`echo ${F} | sed 's/.log$//' | sed 's/^.*\///'`
    if [ `cut -f 2 ${BASEDIR}/table.tsv | grep -c -w "${ID}"` -eq 0 ]; then continue; fi
    PROJECT=`echo ${F} | sed 's/\/[^\/]*$//' | sed 's/^.*\///'`
    AFFECTED="no"
    INST=0
    NUM=0
    if [ -f rayas/${PROJECT}/${ID}.bed ]
    then
	if [ -f rayas/${PROJECT}/${ID}.control.bam ]
	then
	    if [ -f rayas/${PROJECT}/${ID}.tumor.bam ]
	    then
		# min. degree 100
		INST=`tail -n +2 rayas/${PROJECT}/${ID}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq  | wc -l | cut -f 1`
		NUM=0
		if [ ${INST} -gt 0 ]
		then
		    for CLUSTID in `tail -n +2 rayas/${PROJECT}/${ID}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq`
		    do
			SEGMENTS=`cat rayas/${PROJECT}/${ID}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 7 | awk '{SUM+=$1; BASE+=2;} END {print (SUM-BASE);}'`
			echo ${SEGMENTS} | sed "s/^/${ID}\t${CLUSTID}\t/" >> stats.ti.seglen
			cat rayas/${PROJECT}/${ID}.bed | awk '$8=="'${CLUSTID}'"' | sed "s/^/${ID}\t/" >> stats.ti.tsv
			NUMCLUST=`cat rayas/${PROJECT}/${ID}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 1-3 | sort | uniq  | wc -l | cut -f 1`
			NUM=`expr ${NUM} + ${NUMCLUST}`
		    done
		    AFFECTED="yes"
		fi		
	    fi
	fi
    fi
    echo -e "${ID}\t${PROJECT}\t${AFFECTED}\t${INST}\t${NUM}" >> stats.tsv
done

# Plot
Rscript stats.R stats.tsv

# Chromothripsis overlap
Rscript chromo_tithread.R table.tsv stats.tsv

# Liposarcoma
cat stats.ti.tsv | grep -w -Ff <(grep "SARC-US" table.tsv | grep "SoftTissue-Liposarc" | cut -f 2) > lipo.sarc-us.tmpl_ins.bed
Rscript lipo.chr.R
rm lipo.sarc-us.tmpl_ins.bed

# Statistics
echo "Number of samples with TI theads"
cat stats.tsv | cut -f 3 | grep "yes" | wc -l
echo "Number of instances"
cat stats.tsv | awk '$4!=0' | tail -n +2 | awk '{SUM+=$4;} END {print SUM}'
echo "Source segments per instance"
cat stats.tsv | awk '$4!=0' | tail -n +2 | awk '{SUM+=$5/$4;} END {print SUM/NR}'
echo "Median source segment size"
cat stats.ti.tsv | cut -f 1-4 | awk '{print $0"\t"$4-$3;}' | sort -k5,5n | awk '{count[NR]=$5;} END {print count[int(NR/2)];}'
echo "Median number of concatenated segments"
cat stats.ti.seglen | sort -k3,3n | awk '{count[NR]=$3;} END {print count[int(NR/2)];}'
