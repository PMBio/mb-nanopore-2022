#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate visualization

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Templated insertion thread (1st instance, example reads)
echo -e "0645e9e8-3cc3-4f9f-b6f1-382412b524b6\n5cce2684-dd22-40c5-8e2b-97164edfb753\nb6678d21-0115-4a5a-81b8-eb5e8e269890" > reads.lst

# Templated insertion thread (2nd instance, example reads)
echo -e "f60d5507-ab9b-4617-a08f-a93b5639e80a\nfd080231-152f-4557-b45c-b3ebfd941791\nfea65b28-fe0e-40f7-acca-56ba21839675" >> reads.lst

# Plot
wally matches -s -R reads.lst -g ${HG} ${BASEDIR}/../../../alignment/ont/Primary_tumor.bam

# Reduce to high-copy regions
#wally matches -s -n 1000 -x 800 -y 10 -m 10 -R reads.lst -g ${HG} ${BASEDIR}/../../../alignment/ont/Primary_tumor.bam
