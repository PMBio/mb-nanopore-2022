#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate assembly

READ=${1}
SIZE=${2}
LEN=`tail -n 1 ${READ}.fa | awk '{print length($1);}'`
OFFSET=`echo "0.1 * ${SIZE}" | bc -l | awk '{print int($1);}'`
MATCH=5
MISMATCH=-4
GOP=-4
GEX=-4

MAXHITS=0
MAXBEG=0
rm -f hit.${SIZE}.collector
for BEG in `seq 1 ${OFFSET} ${LEN}`
do
    END=`expr ${BEG} + ${SIZE}`
    if [ ${END} -gt ${LEN} ]
    then
	break
    fi
    samtools faidx ${READ}.fa ${READ}:${BEG}-${END} | sed 's/^>.*$/>sub/' > sub.${SIZE}.fwd.fa
    samtools faidx -i ${READ}.fa ${READ}:${BEG}-${END} | sed 's/^>.*$/>sub/' > sub.${SIZE}.rev.fa
    cp ${READ}.fa ${READ}.${SIZE}.fa
    
    # Find all hits
    RUN=1
    NUMHITS=0
    while [ ${RUN} -eq 1 ]
    do
	alfred pwalign -q -m ${MATCH} -n ${MISMATCH} -g ${GOP} -e ${GEX} -a al.${SIZE}.fwd.fa.gz ${READ}.${SIZE}.fa sub.${SIZE}.fwd.fa > /dev/null && gunzip -f al.${SIZE}.fwd.fa.gz
	python score.py -a al.${SIZE}.fwd.fa -m ${MATCH} -n ${MISMATCH} -g ${GOP} -e ${GEX} > hit.${SIZE}.fwd.txt
	SCOREFWD=`cut -f 4 hit.${SIZE}.fwd.txt`
	alfred pwalign -q -m ${MATCH} -n ${MISMATCH} -g ${GOP} -e ${GEX} -a al.${SIZE}.rev.fa.gz ${READ}.${SIZE}.fa sub.${SIZE}.rev.fa > /dev/null && gunzip -f al.${SIZE}.rev.fa.gz
	python score.py -a al.${SIZE}.rev.fa -m ${MATCH} -n ${MISMATCH} -g ${GOP} -e ${GEX} > hit.${SIZE}.rev.txt
	SCOREREV=`cut -f 4 hit.${SIZE}.rev.txt`
	if [ ${SCOREFWD} -gt ${SCOREREV} ]
	then
	    cat hit.${SIZE}.fwd.txt | sed 's/$/\tfwd/' > hit.${SIZE}.txt
	    mv al.${SIZE}.fwd.fa al.${SIZE}.fa
	    SCORE=${SCOREFWD}
	else
	    cat hit.${SIZE}.rev.txt | sed 's/$/\trev/' > hit.${SIZE}.txt
	    mv al.${SIZE}.rev.fa al.${SIZE}.fa
	    SCORE=${SCOREREV}
	fi
	rm -f hit.${SIZE}.fwd.txt hit.${SIZE}.rev.txt al.${SIZE}.fwd.fa al.${SIZE}.rev.fa
	if [ ${SCORE} -gt 0 ]
	then
	    # Debug
	    #cat hit.${SIZE}.txt >> hit.${SIZE}.collector
	    # Mask hit
	    bedtools maskfasta -fi ${READ}.${SIZE}.fa -fo ${READ}.${SIZE}.fa.tmp -bed hit.${SIZE}.txt
	    mv ${READ}.${SIZE}.fa.tmp ${READ}.${SIZE}.fa
	    NUMHITS=`expr ${NUMHITS} + 1`
	else
	    RUN=0
	fi
	# Debug
	#cat hit.${SIZE}.txt
	#echo -e "${READ}\t${LEN}\t${SIZE}\t${BEG}\t${NUMHITS}"
	rm hit.${SIZE}.txt al.${SIZE}.fa
    done
    if [ ${NUMHITS} -gt ${MAXHITS} ]
    then
	MAXHITS=${NUMHITS}
	MAXBEG=${BEG}
    fi
    rm ${READ}.${SIZE}.fa
    rm sub.${SIZE}.fwd.fa sub.${SIZE}.rev.fa
done
echo -e "${READ}\t${LEN}\t${SIZE}\t${MAXBEG}\t${MAXHITS}"
