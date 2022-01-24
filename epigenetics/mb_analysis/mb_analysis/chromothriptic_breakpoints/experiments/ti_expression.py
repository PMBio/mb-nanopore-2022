import pandas as pd
from mb_analysis.config import module_config
from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature, GeneNameToEnsemblID
from mb_analysis.summary.summarize import get_diffexp

"""Check differential expression for genes near templated insertions"""

if __name__ == "__main__":
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    df = pd.read_csv(module_config.templated_insertions_bed_file_nov2021, sep="\t")
    df["chr"] = df["chr"].map(lambda x: x.replace("chr", ""))
    
    rna_diff = get_diffexp()
    
    for index, ti in df.iterrows():
        for gene in gff.chromosomes[ti["chr"]].get_in_range(ti["start"] - 1e4, ti["end"] + 1e4):
            geneid = gene.sanitized_id()
            if geneid in rna_diff.index:
                print(ti["chr"], ti["start"], ti["end"], geneid, gene.name, rna_diff.loc[geneid])
