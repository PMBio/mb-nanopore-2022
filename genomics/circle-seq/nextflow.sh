#!/bin/bash

if [ -f nextflow ]
then
    wget -qO- https://get.nextflow.io | bash
fi

# Run circdna workflow
./nextflow run nf-core/circdna --input sample_sheet.csv --outdir ecDNA --genome GRCh38 -profile singularity --circle_identifier circle_map_realign,circle_map_repeats,circle_finder,circexplorer2
