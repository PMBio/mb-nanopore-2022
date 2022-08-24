#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa


if [ ! -d strelka_somatic ]
then
    ./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${BASEDIR}/../../alignment/illumina/blood.bam --tumorBam ${BASEDIR}/../../alignment/illumina/tumor.bam --referenceFasta ${HG} --runDir strelka_somatic
    python strelka_somatic/runWorkflow.py -m local -j 4
fi

if [ ! -d strelka_relapse_somatic ]
then
    ./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${BASEDIR}/../../alignment/illumina/blood.bam --tumorBam ${BASEDIR}/../../alignment/illumina/relapse.bam --referenceFasta ${HG} --runDir strelka_relapse_somatic
    python strelka_relapse_somatic/runWorkflow.py -m local -j 4
fi
