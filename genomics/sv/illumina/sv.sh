#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate truvari

HG=${BASEDIR}/../../genome/hg38.fa

# Manta
if [ ! -d manta-1.6.0.centos6_x86_64 ]
then
    wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
    tar -xf manta-1.6.0.centos6_x86_64.tar.bz2
    rm manta-1.6.0.centos6_x86_64.tar.bz2
fi
./manta-1.6.0.centos6_x86_64/bin/configManta.py --normalBam ${BASEDIR}/../../alignment/illumina/blood.bam --tumorBam ${BASEDIR}/../../alignment/illumina/tumor.bam --referenceFasta ${HG} --runDir manta.tumor
./manta.tumor/runWorkflow.py -j 8

exit;

# Delly
delly call -g ${HG} -x human.hg38.excl.tsv -o delly.bcf ${BASEDIR}/../../alignment/illumina/*.bam

# Install svaba
if [ ! -d svaba ]
then
    git clone --recursive https://github.com/walaj/svaba
    cd svaba
    ./configure
    make
    make install
    ./bin/svaba 
    cd ..
fi

# SVaba
./svaba/bin/svaba run -t ${BASEDIR}/../../alignment/illumina/tumor.bam -t ${BASEDIR}/../../alignment/illumina/relapse.bam -n ${BASEDIR}/../../alignment/illumina/blood.bam -G ${HG} -a tumor.svaba
for F in *.vcf
do
    sed -i 's/\/[^\t]*\/..\/..\/alignment\/illumina\/\([a-z]*\).bam/\1/g' ${F}
    grep -m 1 "^#CHROM" ${F} | cut -f 10-
    bgzip ${F}
    tabix ${F}.gz
done
