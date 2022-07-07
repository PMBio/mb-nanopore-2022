#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate sv

HG=${BASEDIR}/../../genome/hg38.fa

# Illumina
bcftools view ${BASEDIR}/../illumina/germline.bcf | grep -v "SVTYPE=INS" | bgzip > illumina.vcf.gz
tabix illumina.vcf.gz

# ONT
bcftools view ${BASEDIR}/../ont/germline.bcf | grep -v "SVTYPE=INS" | bgzip > ont.vcf.gz
tabix ont.vcf.gz

# Comparison
rm -rf germ_stats/
truvari bench -p 0 --gtcom --passonly --no-ref a -r 1000 -C 1000 -b ont.vcf.gz -c illumina.vcf.gz -f ${BASEDIR}/../../genome/hg38.fa -o germ_stats
cat germ_stats/tp-base.vcf  | grep -v "^#" | awk '{print "Common"NR;}' > ont.calls
cat germ_stats/tp-base.vcf  | grep -v "^#" | awk '{print "Common"NR;}' > illumina.calls
cat germ_stats/fp.vcf | grep -v "^#" | awk '{print "IlluminaOnly"NR;}' >> illumina.calls
cat germ_stats/fn.vcf | grep -v "^#" | awk '{print "OntOnly"NR;}' >> ont.calls
#truvari bench -p 0 --gtcom --passonly --no-ref a -r 1000 -C 1000 -b illumina.vcf.gz -c ont.vcf.gz -f ${BASEDIR}/../../genome/hg38.fa -o germ_stats

# Size distribution
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" illumina.vcf.gz | awk '{print $0"\t"($3-$2)"\tillumina";}' > size.tsv
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ont.vcf.gz | awk '{print $0"\t"($3-$2)"\tont";}' >> size.tsv

# Plotting
Rscript plot.R
Rscript venn.R
rm -rf illumina.vcf.gz* ont.vcf.gz* size.tsv ont.calls illumina.calls germ_stats/
