from mb_analysis.config import module_config
from nanoepitools.annotations.annotations import GFFAnnotationsReader

"""
Little helper script creating a bed file with promoter regions for all transcripts
"""


if __name__=="__main__":
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()

    for chrom in module_config.chroms:
        for gene in gff.chromosomes[chrom].children.values():
            for transcript in gene.children.values():
                if transcript.direction == "+":
                    print(f"{chrom}\t{transcript.start-2000}\t{transcript.start+500}")
                else:
                    print(f"{chrom}\t{transcript.end - 500}\t{transcript.end + 2000}")