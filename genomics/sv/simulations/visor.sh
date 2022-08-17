#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate sv

HG=${BASEDIR}/../../genome/hg38.fa
NTITHREADS=10
NSOURCES=3
NSEGMENTS=100

# Simulate using chr18
samtools faidx ../../genome/hg38.fa chr18 > chr18.fa
samtools faidx chr18.fa

# Benchmark SV calling
echo -e "svtype\tmode\tcoverage\treadlen\trecall\tprecision\tf1" > summary.stats.tsv
SVT=TITHREAD
for COV in 5 10 20 30 40
do
    #pb ont
    for MODE in ill pb ont
    do
	LEN=150
	if [ ${MODE} == "ont" ]
	then
	    MODEL="nanopore2018"
	    RT="nanopore"
	    LEN=15000
	fi
	if [ ${MODE} == "pb" ]
	then
	    MODEL="pacbio2016"
	    RT="pacbio"
	    LEN=15000
	fi
	# Simulate anew
	if [ ! -d sim_svt${SVT}_${MODE}_cov${COV}_len${LEN} ]
	then
	    mkdir sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}

	    # Control
	    cat chr18.fa.fai  | awk '{print $1"\t1\t"$2"\t100.0\t100.0";}' > sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/simulate.bed
	    touch sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap1.bed sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap2.bed
	    VISOR HACk -g chr18.fa -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap1.bed sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap2.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_control
	    
	    # Simulate haplotypes
	    for J in `seq 1 1 ${NTITHREADS}`
	    do
		rm -f templates.fa
		for I in `seq 1 1 ${NSOURCES}` 
		do
		    ST1=`echo $(($RANDOM%30000)) | awk '{print ($1*1000 + 1000000);}'`
		    SIZE=`echo $(($RANDOM%500)) | awk '{print ($1 + 250);}'`
		    ED1=`expr ${ST1} + ${SIZE}`
		    SEQUENCE=`samtools faidx chr18.fa chr18:${ST1}-${ED1} | tail -n +2 | tr '\n' '#' | sed 's/#//g'`
		    echo -e "chr18\t${ST1}\t${ED1}\tSegment${I}\tThread${J}" >> sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.bed
		    echo ${SEQUENCE} >> templates.fa
		done
		rm -f tithread.fa
		for I in `seq 1 1 ${NSEGMENTS}`
		do
		    CHOICE=`echo $(($RANDOM%${NSOURCES})) | awk '{print ($1 + 1);}'`
		    FLIP=`echo $(($RANDOM%2))`
		    if [ ${FLIP} -eq 1 ]
		    then
			sed "${CHOICE}q;d" templates.fa | tr 'acgtACGT' 'tgcaTGCA' | rev >> tithread.fa
		    else
			sed "${CHOICE}q;d" templates.fa >> tithread.fa
		    fi
		done
		ST=`echo $(($RANDOM%30000)) | awk '{print ($1*1000 + 50000000);}'`
		ED=`expr ${ST} + 1`
		SEQUENCE=`cat tithread.fa | tr '\n' '#' | sed 's/#//g'`
		rm -f templates.fa tithread.fa
		echo -e "chr18\t${ST}\t${ED}\tinsertion\t${SEQUENCE}\t0" >> sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap1.bed
	    done
	    VISOR HACk -g chr18.fa -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap1.bed sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.hap2.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_tumor
	    # Draw reads
	    if [ ${MODE} == "ill" ]
	    then
		VISOR SHORtS --coverage ${COV} --length ${LEN} -g chr18.fa -s sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_tumor/ -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/simulate.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/tumor
		VISOR SHORtS --coverage ${COV} --length ${LEN} -g chr18.fa -s sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_control/ -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/simulate.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/control
	    else
		VISOR LASeR --coverage ${COV} --length_mean ${LEN} --error_model ${MODEL} --qscore_model ${MODEL} --read_type ${RT} -g chr18.fa -s sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_tumor/ -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/simulate.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/tumor
		VISOR LASeR --coverage ${COV} --length_mean ${LEN} --error_model ${MODEL} --qscore_model ${MODEL} --read_type ${RT} -g chr18.fa -s sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/hack_control/ -b sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/simulate.bed -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/control  
	    fi
	fi
	
	# Call SVs
	if [ ! -f sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed ]
	then
	    if [ ${MODE} == "ill" ]
	    then
		/opt/dev/rayas/bin/rayas call -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed -g chr18.fa -m sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/control/sim.srt.bam sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/tumor/sim.srt.bam
		cut -f 4,9 sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.dot
		dot -Tpdf sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.dot -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.pdf
	    else
		echo ${MODE}
		#/opt/dev/delly/bin/delly lr -y ${MODE} -o sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.bcf -g data/chr18.fa sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/out/sim.srt.bam
	    fi
	fi

	# Evaluation
	if [ -f sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed ]
	then
	    bedtools intersect -a sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.bed -b <(tail -n +2 sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed) -wao | awk '$6!="."' | cut -f 6- | sort -k1,1V -k2,2n | uniq > sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/true_positives.bed
	    TOTTRUTH=`cat sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/${SVT}.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
	    TOTCALLED=`tail -n +2 sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/results.bed | cut -f 1-3 | sort | uniq | wc -l | cut -f 1`
	    TP=`cut -f 1-3 sim_svt${SVT}_${MODE}_cov${COV}_len${LEN}/true_positives.bed | sort | uniq | wc -l | cut -f 1`
	    RECALL=`echo "${TP} / ${TOTTRUTH}" | bc -l`
	    PREC=0
	    if [ ${TOTCALLED} -gt 0 ]
	    then
		PREC=`echo "${TP} / ${TOTCALLED}" | bc -l`
	    fi
	    F1=0
	    if [ ${RECALL} != "0" ]
	    then
		if [ ${PREC} != "0" ]
		then
		    F1=`echo "2 * (${RECALL} * ${PREC}) / (${RECALL} + ${PREC})" | bc -l`
		fi
	    fi
	    echo -e "${SVT}\t${MODE}\t${COV}\t${LEN}\t${RECALL}\t${PREC}\t${F1}" >> summary.stats.tsv
	fi
    done
done
rm chr18.fa*

# PLot
Rscript plot.R
