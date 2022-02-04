#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

#HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# SV plots for tumor or relapse
for TUMOR in tumor 
do
    if [ -f ${BASEDIR}/../../sv/${TUMOR}.somatic.bcf ]
    then
	SVFILE=${BASEDIR}/../../sv/${TUMOR}.somatic.bcf
	# Intra-chromosomal SVs <50kbp
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ${SVFILE} | grep -v "BND" | awk '$3-$2<50000' | sed "s/$/\t${SAMPLE}\tblood/" > ${TUMOR}.bed
	# Intra-chromosomal SVs >= 50kbp
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ${SVFILE} | grep -v "BND" | awk '$3-$2>=50000' | awk '{print $1"\t"($2-10)"\t"($2+10)"\t"$4"Left\n"$1"\t"($3-10)"\t"($3+10)"\t"$4"Right";}' | sed "s/$/\t${SAMPLE}\tblood/" >> ${TUMOR}.bed
	# Inter-chromosomal translocations (BND)
	bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SVTYPE\n" ${SVFILE} | grep "BND"  | awk '{print $1"\t"($2-10)"\t"($2+10)"\t"$5"Left\n"$3"\t"($4-10)"\t"($4+10)"\t"$5"Right";}' | sed "s/$/\t${SAMPLE}\tblood/" >> ${TUMOR}.bed

	# Sort and index
	cat ${TUMOR}.bed | sort -k1,1V -k2,2n > ${TUMOR}.bed.tmp
	mv ${TUMOR}.bed.tmp ${TUMOR}.bed
	bgzip ${TUMOR}.bed
	tabix -p bed ${TUMOR}.bed.gz

	# Add 500bp buffer for viewing
	zcat ${TUMOR}.bed.gz | awk '{print $1"\t"($2-500)"\t"($3+500)"\t"$4;}' | awk '$2>=0' > ${TUMOR}.bed

	# Wally (Illumina data)
	#/opt/dev/wally/src/wally region -y 2048 -b ${TUMOR}.bed.gz -g ${HG} -R ${TUMOR}.bed -p -u -c ${BASEDIR}/../../alignment/illumina/*.bam
	# or using ONT data
	/opt/dev/wally/src/wally region -y 2048 -b ${TUMOR}.bed.gz -g ${HG} -R ${TUMOR}.bed -p -u -c ${BASEDIR}/../../alignment/ont/*.bam
    fi
done
