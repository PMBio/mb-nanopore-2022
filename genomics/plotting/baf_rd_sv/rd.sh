#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Whole chromosome plots
for TUMOR in relapse tumor
do
    if [ -f ${BASEDIR}/../../coverage/illumina/${TUMOR}.cov.gz ]
    then
	COV=${BASEDIR}/../../coverage/illumina/${TUMOR}.cov.gz
	if [ -f ${BASEDIR}/../../sv/${TUMOR}.somatic.bcf ]
	then
	    SVFILE=${BASEDIR}/../../sv/${TUMOR}.somatic.bcf
	    if [ -f ${BASEDIR}/../../phasing/split_rephase/blood.phased.bcf ]
	    then
		PHASED=${BASEDIR}/../../phasing/split_rephase/blood.phased.bcf
		VCF=${BASEDIR}/../../snv/illumina/freebayes.vcf.gz
		
		# Fetch delly SVs
		if [ ! -f ${TUMOR}.sv ]
		then
		    bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%INFO/SVTYPE\n" ${SVFILE} | grep -v "BND" > ${TUMOR}.sv 
		    bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%INFO/SVTYPE\n" ${SVFILE} | grep "BND" >> ${TUMOR}.sv
		fi
		
		# Fetch SNVs
		if [ ! -f ${TUMOR}.vaf.tsv ]
		then
		    bcftools view -v snps -m2 -M2 -s blood ${PHASED} | cut -f 1,2,4,5,10 | sed 's/:.*$//' | grep "1|0$" | cut -f 1-4 > ${TUMOR}.snp.sites
		    bcftools view -v snps -m2 -M2 -s ${TUMOR} ${VCF} | cut -f 1-9,10 | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n" - | sed 's/,/\t/' | grep -w -Ff ${TUMOR}.snp.sites | awk '($5+$6)>0 {print $1"\t"$2"\t"$6/($5+$6);}' > ${TUMOR}.alt.vaf.tsv
		    rm ${TUMOR}.snp.sites
		    bcftools view -v snps -m2 -M2 -s blood ${PHASED} | cut -f 1,2,4,5,10 | sed 's/:.*$//' | grep "0|1$" | cut -f 1-4 > ${TUMOR}.snp.sites
		    bcftools view -v snps -m2 -M2 -s ${TUMOR} ${VCF} | cut -f 1-9,10 | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n" - | sed 's/,/\t/' | grep -w -Ff ${TUMOR}.snp.sites | awk '($5+$6)>0 {print $1"\t"$2"\t"$5/($5+$6);}' > ${TUMOR}.ref.vaf.tsv
		    rm ${TUMOR}.snp.sites
		    cat ${TUMOR}.alt.vaf.tsv ${TUMOR}.ref.vaf.tsv | sort -k1,1V -k2,2n | uniq > ${TUMOR}.vaf.tsv
		    wc -l ${TUMOR}.alt.vaf.tsv ${TUMOR}.ref.vaf.tsv ${TUMOR}.vaf.tsv
		    rm ${TUMOR}.alt.vaf.tsv ${TUMOR}.ref.vaf.tsv
		fi
		
		# Plot
		for CHR in `seq 1 22`
		do
		    echo "Plotting" chr${CHR}
		    Rscript covWithSV.R ${COV} ${TUMOR}.vaf.tsv ${TUMOR}.sv chr${CHR}:0-250000000
		done
	    fi
	fi
    fi
done

# Customized plots
TUMOR=tumor

# Templated insertion thread, chr7:20507358
Rscript covWithoutSV.R ${BASEDIR}/../../coverage/illumina/${TUMOR}.cov.gz ${TUMOR}.vaf.tsv chr7:20000000-21000000

# Templated insertion thread, chr5:17208573
Rscript covWithoutSV.R ${BASEDIR}/../../coverage/illumina/${TUMOR}.cov.gz ${TUMOR}.vaf.tsv chr5:16700000-17700000

# Chromothripsis example (chr5p)
Rscript covWithoutSV.R ${BASEDIR}/../../coverage/illumina/${TUMOR}.cov.gz ${TUMOR}.vaf.tsv chr5:0-30000000
