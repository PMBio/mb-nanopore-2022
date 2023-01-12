import tqdm
import gzip
from matplotlib.patches import Patch
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mb_analysis.config import module_config
from mb_analysis.reference_cpgs import ReferenceCpGs, format_cpg
from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.figures.dmr_overview import load_asm_hits, annotate_asm_hits

def classify_dmr(cgi_df, chrom, positions):
    chrom_df = cgi_df.get_group(chrom)
    num_cgi = 0
    num_shore = 0
    num_shelf = 0
    num_opensea = 0
    for pos in positions:
        pos = int(pos.split(":")[1])
        if any((chrom_df["start"] <= pos) & (chrom_df["end"] > pos)):
            num_cgi += 1
        elif any((chrom_df["start"]-2000 <= pos) & (chrom_df["end"]+2000 > pos)):
            num_shore += 1
        elif any((chrom_df["start"]-4000 <= pos) & (chrom_df["end"]+4000 > pos)):
            num_shelf += 1
        else:
            num_opensea += 1
    return num_cgi, num_shore, num_shelf, num_opensea
    
types = ["CGI", "Shore", "Shelf", "Open sea"]
type_colors = ["#721817", "#FA9F42", "#2B4162", "#0B6E4F"]
fig_kwargs={"figsize":(8,6)}

def plot_counts(counts, title):
    pa.figure(**fig_kwargs)
    plt.title(title)
    sums = counts.sum(axis=0)
    for i, s in enumerate(sums):
        plt.bar(i, s)
    plt.xticks([0,1,2,3], types)
    plt.ylabel("Number of CpGs")
    pa.savefig()

def classify_row(counts_row):
    counts_row = counts_row / counts_row.sum()
    for i in range(len(counts_row)):
        if counts_row[i] > 0.9:
            return i
    return -1
    
def plot_counts_per_dmr(hits, counts, title):
    pa.figure(**fig_kwargs)
    plt.title(title)
    sums_cgrelated = np.zeros(4)
    sums_cgunrelated = np.zeros(4)
    for hit, counts_row in zip(hits, counts):
        i = classify_row(counts_row)
        if hit["is_cg_variant_related"]:
            sums_cgrelated[i] += 1
        else:
            sums_cgunrelated[i] += 1
    
    for i in range(len(sums_cgrelated)):
        plt.bar(i, sums_cgrelated[i], color="r")
        plt.bar(i, sums_cgunrelated[i], bottom=sums_cgrelated[i], color="grey")
    
    plt.xticks([0,1,2,3], types)
    legend_elements = [Patch(facecolor='r', edgecolor='r', label='CpG-variant related'),
                       Patch(facecolor='gray', edgecolor='gray', label='Not CpG-variant related')]
    plt.gca().legend(handles=legend_elements)
    plt.ylabel("Number of segments")
    pa.savefig()


def plot_effect_distribution(hits, counts, title):
    pa.figure(**fig_kwargs)
    plt.title(title)
    for i in range(4):
        y = []
        for hit, count_row in zip(hits, counts):
            y = y + [abs(hit["diff"])] * count_row[i]
        plt.boxplot(y, positions=[i])
    plt.xticks([0, 1, 2, 3], types)
    pa.savefig()


def vcf_annotate_cg(vcf_df):
    from pyfaidx import Fasta
    ref = Fasta(module_config.reference_fasta_file, "r")
    
    def is_variant_cg_related(row):
        ref_seq = ref[row["chr"]][row["pos"] - 2: row["pos"] + len(row["REF"])].seq
        if "CG" in ref_seq:
            return True
        else:
            return False
    
    vcf_df["is_cg_related"] = vcf_df.apply(is_variant_cg_related, axis=1)

if __name__ == "__main__":
    pa = PlotArchiver("summary", config=module_config)
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    ref_cpgs = ReferenceCpGs()
    
    min_diff = 0.5
    
    samplecomp_hits = []
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        samplecomp_hits += [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"], "diff": line["diff"]} for line in
            pm.read_file(b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=min_diff)]
    
    samplecomp_hits = merge_duplicate_diffmet_hits(samplecomp_hits)
    
    for hit in tqdm.tqdm(samplecomp_hits):
        hit["CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=True, formatted=True)
        hit["complex_CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=False, formatted=True)
    
    """
    Read and init variants
    """
    vcf_blood = PhasedVCF(module_config.vcf_blood_file)
    vcf_blood.read()
    vcf_df = vcf_blood.vcf_df.reset_index()
    
    bs_masked_variants = vcf_df.apply(
        lambda x: (x["ALT"] == "T" and x["REF"] == "C") or (x["ALT"] == "C" and x["REF"] == "T"), axis=1)
    
    useable_variants = vcf_df.loc[~bs_masked_variants]

    vcf_annotate_cg(vcf_df)
    
    vcf_df_grouped = vcf_df.groupby("chr")
    
    """ Load ASM"""
    sample_order = ("Primary", "Relapse")
    asm_hits = load_asm_hits(sample_order, min_diff)
    asm_hits = annotate_asm_hits(asm_hits, vcf_df_grouped, ref_cpgs)
    
    cgi_df = pd.read_csv(module_config.goldenpath_cgi_annotation_file, sep="\t", usecols=[1,2,3], names=["chrom", "start", "end"])

    asm_counts = {sample: np.array([classify_dmr(cgi_df.groupby("chrom"), hit["chrom"], hit["CpGs"]) for hit in asm_hits[sample]]) for sample in asm_hits}
    samplecomp_counts = np.array([classify_dmr(cgi_df.groupby("chrom"), hit["chrom"], hit["CpGs"]) for hit in samplecomp_hits])


    def annotate_hits_cg_variant_related(vcf_df, hits):
        vcf_only_cg_related = vcf_df.loc[vcf_df["is_cg_related"]]
        vcf_only_cg_related_grouped = vcf_only_cg_related.groupby("chr")
        for hit in hits:
            vcf_chrom = vcf_only_cg_related_grouped.get_group(hit["chrom"])
            hit["is_cg_variant_related"] = any((vcf_chrom["pos"] > hit["start"]) & (vcf_chrom["pos"] < hit["end"]))


    annotate_hits_cg_variant_related(vcf_df, samplecomp_hits)
    annotate_hits_cg_variant_related(vcf_df, asm_hits["Primary"])
    annotate_hits_cg_variant_related(vcf_df, asm_hits["Relapse"])

    print(sum(c["is_cg_variant_related"] for c in samplecomp_hits) / len(samplecomp_hits))
    print(sum(c["is_cg_variant_related"] for c in asm_hits["Primary"]) / len(asm_hits["Primary"]))
    print(sum(c["is_cg_variant_related"] for c in asm_hits["Relapse"]) / len(asm_hits["Relapse"]))
    
    with pa.open_multipage_pdf("dmr_classification"):
        plot_counts(samplecomp_counts, "DMR Primary vs Relapse")
        for sample in asm_counts:
            plot_counts(asm_counts[sample], f"ASM {sample}")

        plot_counts_per_dmr(samplecomp_hits, samplecomp_counts, "DMR Primary vs Relapse")
        for sample in asm_counts:
            plot_counts_per_dmr(asm_hits[sample], asm_counts[sample], f"ASM {sample}")
            
        plot_effect_distribution(samplecomp_hits, samplecomp_counts, "DMR Primary vs Relapse")
        for sample in asm_counts:
            plot_effect_distribution(asm_hits[sample], asm_counts[sample], f"ASM {sample}")

        
