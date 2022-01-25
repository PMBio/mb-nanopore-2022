import pysam
import pandas as pd
import matplotlib.pyplot as plt
from nanoepitools.plotting.general_plotting import PlotArchiver
from mb_analysis.config import module_config
from nanoepitools.annotations.annotations import GFFAnnotationsReader
import numpy as np


"""
The question was raised whether long insertions imagined by the gene fusion explanation via SV tracing are
realistic. Here we plot the distribution of intron lengths.
 
It shows that they are indeed realistic.
 """


def compute_intron_length_distribution(gff):
    all_distances = []
    for chrom in gff.chromosomes:
        chrom = gff.chromosomes[chrom]
        for gene in chrom.children:
            gene = chrom.children[gene]
            # To be more conservatives we even remove long non-coding rna
            if not gene.name.startswith("LINC"):
                continue
            for transcript in gene.children:
                transcript = gene.children[transcript]
                exons = transcript.leafs
                exons = sorted(exons, key=lambda x: x.start)
                if len(exons) == 1:
                    continue
                exons1 = exons[:-1]
                exons2 = exons[1:]
                
                assert len(exons1) == len(exons2)
                dist = [b.start - a.end for a, b in zip(exons1, exons2)]
                all_distances += dist
    return np.array(all_distances)


if __name__ == "__main__":
    pa = PlotArchiver("summary", config={"plot_archive_dir": "/homes/snajder/data1/plots_medulloblastoma"})
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    intron_len_dist = compute_intron_length_distribution(gff)
    
    with pa.open_multipage_pdf("intron_length_distribution"):
        pa.figure()
        plt.hist(intron_len_dist, bins=100)
        pa.savefig()
        plt.hist(np.log10(intron_len_dist), bins=100)
        pa.savefig()
