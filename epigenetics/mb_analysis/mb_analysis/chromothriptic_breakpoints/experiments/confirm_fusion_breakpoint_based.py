import pysam
import pandas as pd
import matplotlib.pyplot as plt
from nanoepitools.plotting.general_plotting import PlotArchiver
from mb_analysis.config import module_config
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.chromothriptic_breakpoints.breakpoints import ChromothripticBreakpoints
import numpy as np


pa = PlotArchiver("summary", config={"plot_archive_dir": "/homes/snajder/data1/plots_medulloblastoma"})

gff = GFFAnnotationsReader()
gff.read(module_config.gff_file, only_protein_coding=False)
gff.build_index()



def validate_by_genomic_bp(row, max_dist):
    chrom_a, pos_a = row["breakpoint1"].split(":")
    chrom_b, pos_b = row["breakpoint2"].split(":")
    pos_a = int(pos_a)
    pos_b = int(pos_b)
    path = bps.get_shortest_path(chrom_a, pos_a, chrom_b, pos_b, max_dist=max_dist)
    return path


def validate_by_read_support(bam, row, window):
    chrom_a, pos_a = row["breakpoint1"].split(":")
    chrom_b, pos_b = row["breakpoint2"].split(":")
    pos_a = int(pos_a)
    pos_b = int(pos_b)
    reads_a = {q.query_name for q in bam.fetch(chrom_a, max(0, pos_a - window // 2), pos_a + window // 2)}
    reads_b = {q.query_name for q in bam.fetch(chrom_b, max(0, pos_b - window // 2), pos_b + window // 2)}
    return len(reads_a.intersection(reads_b)) > 0

def get_reads_for_gene(chrom, pos, window):
    if chrom_a == chrom_b and abs(pos_a - pos_b) < 5e4:
        # These are just nearby genes, obviously we'll have reads covering both
        return False
    
    best_gene = None
    mindist = window
    for gene in gff.chromosomes[chrom].get_in_range(pos - window // 2, pos + window // 2, max_recursion=0):
        dist = min(abs(gene.start - pos), abs(gene.end - pos))
        if dist < mindist:
            mindist = dist
            best_gene = gene
    if best_gene is not None:
        reads = {q.query_name for q in f.fetch(chrom, best_gene.start, best_gene.end)}
    else:
        reads = {q.query_name for q in f.fetch(chrom, pos - window // 2, pos + window // 2)}
    return reads


def validate_by_read_support_in_gene(row, window):
    chrom_a, pos_a = row["breakpoint1"].split(":")
    chrom_b, pos_b = row["breakpoint2"].split(":")
    pos_a = int(pos_a)
    pos_b = int(pos_b)
    if chrom_a == chrom_b and abs(pos_a - pos_b) < 5e4:
        # These are just nearby genes, obviously we'll have reads covering both
        return False
    reads_a = get_reads_for_gene(chrom_a, pos_a, window)
    reads_b = get_reads_for_gene(chrom_b, pos_b, window)
    return len(reads_a.intersection(reads_b)) > 0

bam = "/homes/snajder/data/medulloblastoma/from_reference_incl_supp/mapping/Primary.sorted.bam"


bps = ChromothripticBreakpoints(gff)
bps.load()


fusion_df = pd.read_csv(module_config.fusion_genes_primary_file, sep="\t")

fusion_df["X.gene1"] = fusion_df["X.gene1"].map(lambda x: ",".join([g.split("(")[0] for g in x.split(",")]))
fusion_df["gene2"] = fusion_df["gene2"].map(lambda x: ",".join([g.split("(")[0] for g in x.split(",")]))


ordered_confidence = ["low", "medium", "high"]
with pysam.AlignmentFile(bam, "r") as f, pa.open_multipage_pdf("number_of_fusions_allmethods"):
    for min_confidence in range(len(ordered_confidence)):
        allowed_confidence = ordered_confidence[min_confidence:]
        min_confidence = allowed_confidence[0]
        print(min_confidence, allowed_confidence)
        
        fusion_df_conf = fusion_df.loc[fusion_df["confidence"].map(lambda x: x in allowed_confidence)]
        print(fusion_df_conf.shape[0], "calls")
        grouped = fusion_df_conf.groupby(["X.gene1", "gene2"])
        print(len(grouped.groups), "pairs")
        n_bp_confirmed = 0
        n_read_support = 0
        n_both_supports = 0
        n_not_found = 0
        for group in grouped.groups:
            sub_df = grouped.get_group(group)
            found = False
            for _, row in sub_df.iterrows():
                read_support = validate_by_read_support(f, row, 2e3)
                path = validate_by_genomic_bp(row, 5e5)
                if path.is_connected and path.uses_breakpoint and read_support:
                    #print("Both", row["X.gene1"], row["gene2"])
                    n_both_supports += 1
                    found = True
                    break
                elif path.is_connected and path.uses_breakpoint:
                    #print("Only path", row["X.gene1"], row["gene2"])
                    n_bp_confirmed += 1
                    found = True
                    break
                elif read_support:
                    #print("Only reads", row["X.gene1"], row["gene2"])
                    n_read_support += 1
                    found = True
                    break
            if not found:
                n_not_found += 1
        
        def absolute_value(val):
            total = n_not_found + n_read_support + n_bp_confirmed + n_both_supports
            """Because piechart in matplotlib cant show non-percentages"""
            a = int(np.round(val / 100.0 * total, 0))
            return a
        
        pa.figure()
        plt.title(f"{min_confidence} or higher")
        patches, _,_  = plt.pie(
            [n_both_supports, n_bp_confirmed, n_read_support, n_not_found],
            colors=["#A14D38", "#0492B5", "#A18DF2", "gray"],
            shadow=False,
            autopct=absolute_value,
        )
        plt.legend(patches, labels=[
                "High confidence read support",
                "Explainable using high confidence genomics breakpoints",
                "Supported by individual reads",
                "Unsupported",
            ], loc="best")
        pa.savefig()
