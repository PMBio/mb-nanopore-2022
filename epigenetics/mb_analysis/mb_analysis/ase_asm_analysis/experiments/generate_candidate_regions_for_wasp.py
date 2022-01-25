import pandas as pd

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from mb_analysis.config import module_config

"""
Generates a list of regions to hook ase calls on. Used by WASP. This needed because otherwise WASP will report
based on the haplotype of some variant far far away, and we can no longer trust that the HP assignment of the
ASE calls and ASM calls are consistent.
"""

if __name__ == "__main__":
    out_file = module_config.wasp_ase_candidate_region_file
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    vcf_valid = []
    
    for chrom in module_config.chroms:
        for gene in gff.chromosomes[chrom].sorted_children:
            for transcript in gene.sorted_children:
                if transcript.direction == "+":
                    start = transcript.start - 2000
                    end = transcript.start + 500
                else:
                    end = transcript.end + 2000
                    start = transcript.end - 500
                start = start - 1000
                end = end + 1000
                vcf_valid.append(["chr" + chrom, start, end])
    
    vcf_valid = pd.DataFrame(vcf_valid, columns=["chr", "start", "end"])
    vcf_valid.to_csv(out_file, sep="\t", index=False)
