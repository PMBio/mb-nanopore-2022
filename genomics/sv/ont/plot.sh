#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

# Germline SVs
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" germline.bcf | grep "DEL$" | awk '{print $0"\t"(($3-$2) * -1);}' | cut -f 1,2,4,5,6 > size.tsv
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/SVTYPE\n" germline.bcf | grep "INS$" | awk '{print $0"\t"(length($4) - length($3));}' | cut -f 1,2,5,6,7 >> size.tsv

# Plotting
Rscript plot.R
