import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scipy.stats
from pingouin import partial_corr

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.math import fdr_from_pvals
from mb_analysis.ase_asm_analysis.collective_asm import CollectiveAlleleSpecificMethylation
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from mb_analysis.ase_asm_analysis.ase import ASE
from mb_analysis.ase_asm_analysis.ase_vs_asm_plotting import plot_ase_vs_asm, plot_ase_vs_asm_nosplit
from mb_analysis.ase_asm_analysis.cnv import load_gene_cnv
from mb_analysis.config import module_config

"""
Creates ASEs ASM figure. At the time of writing this is Figure 4b.
"""

def compute_ase_asm(ase_gene_ratio_df, asm_result, gene_cnv_dict, gff, min_abs_diff, sample):
    promoters_hit = asm_result[sample]
    ase_sample = ase_gene_ratio_df.loc[ase_gene_ratio_df.index.map(lambda x: x in promoters_hit)].copy()
    ase_sample["fdr"] = fdr_from_pvals(ase_sample["pval"])
    plot_data = dict(ase=[], asm=[], cnv=[], color=[], gene_id=[], label=[], shape=[])
    
    color_mode = "expr_consistency"
    for gene, row in ase_sample.iterrows():
        gene_name = gff.get_gene(gene).name
        # If there are multiple effects in the promoter, pick the region with the largest effect
        diffmet = promoters_hit[gene]
        if abs(diffmet) < min_abs_diff:
            continue
        
        row = ase_sample.loc[gene]
        
        is_sig = row["fdr"] < 0.1
        ase_ratio = row["hp1_ratio"]
        plot_data["ase"].append(ase_ratio)
        plot_data["shape"].append("*" if is_sig else "o")
        
        plot_data["dmr"].append(diffmet)
        plot_data["label"].append(gene_name)
        plot_data["gene_id"].append(gene)
        
        if gene in gene_cnv_dict[sample]:
            gene_hp1_cnv_ratio = gene_cnv_dict[sample][gene][0] / sum(gene_cnv_dict[sample][gene])
            plot_data["color"].append(
                "r" if gene_hp1_cnv_ratio > 0.65 else "b" if gene_hp1_cnv_ratio < 0.35 else "gray"
            )
            plot_data["cnv"].append(gene_hp1_cnv_ratio)
        else:
            plot_data["color"].append("gray")
            plot_data["cnv"].append(np.nan)
    plot_data = {k: np.array(v) for k, v in plot_data.items()}
    
    return plot_data


def main():
    pa = PlotArchiver("mb_asm_ase", config=module_config)
    gene_to_variant_mapping_mode = "nearest"
    samples = ["Primary"]
    min_abs_diff = 0.25
    
    """ load data """
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    with CollectiveAlleleSpecificMethylation(gff, samples=samples) as asm:
        asm_result = asm.read_all_samples(annotation="promoters", min_diff=min_abs_diff)
    
    vcf_blood = PhasedVCF(module_config.vcf_blood_file)
    vcf_blood.read()
    
    gene_cnv_dict = {
        sample: load_gene_cnv(module_config.cnv_file_template.format(sample=sample), vcf_blood, gff)
        for sample in samples
    }
    
    ase = ASE(module_config.wasp_ase_file)
    ase.load(load_exon_annotation=False, load_gene_annotation=True, pval_thres=None)
    ase.assign_counts_to_hp(vcf_blood, "blood")
    ase.compute_hp1_ratio()
    ase_gene_ratio_df = ase.ase_df.groupby("geneid").agg(lambda x: next(iter(x)))
    
    """ Compute ASE/ASE/ASCNV for each sample separately """
    for sample in samples:
        plot_data = compute_ase_asm(ase_gene_ratio_df, asm_result, gene_cnv_dict, gff, min_abs_diff, sample)
        
        df = pd.DataFrame(
            {
                "gene_id": plot_data["gene_id"],
                "gene_name": plot_data["labels"],
                "ase_hp1": plot_data["ase"],
                "asm_hp1_minus_hp2": -plot_data["dmr"],
                "cn_hp1": plot_data["cnv"],
            }
        )
        df.to_csv(Path(module_config.asm_vs_ase_dir).joinpath(f"{sample}_significant.tsv"), sep="\t", index=False)
        
        with pa.open_multipage_pdf(f"ase_vs_asm_promoter_{sample}"):
            # Don't highlight significance
            plot_ase_vs_asm_nosplit(
                pa,
                plot_data["ase"],
                plot_data["dmr"],
                plot_data["labels"],
                colors=plot_data["color"],
                colorbar=False,
                vmin=0,
                vmax=1,
                hp2_minus_hp1=True,
                min_abs_diff=min_abs_diff,
            )
            plt.suptitle(f"Promoter ASM vs ASE {sample}", y=1)
            pa.savefig()
            
            # Highlight significance
            plot_ase_vs_asm_nosplit(
                pa,
                plot_ase,
                plot_asm,
                [l if s == "*" else "" for l, s in zip(plot_data["labels"], plot_data["shapes"])],
                shape=plot_data["shapes"],
                colors=plot_data["color"],
                colorbar=False,
                vmin=0,
                vmax=1,
                hp2_minus_hp1=True,
                min_abs_diff=min_abs_diff,
            )
            plt.suptitle(f"Promoter ASM vs ASE {sample}", y=1)
            pa.savefig()
        
        df = pd.DataFrame(
            {
                "gene_id": plot_data["gene_ids"],
                "gene_name": plot_data["labels"],
                "ase_hp1": plot_data["ase"],
                "asm_hp1_minus_hp2": -plot_data["dmr"],
                "cn_hp1": plot_data["cnv"],
            }
        )
        df.to_csv(Path(module_config.asm_vs_ase_dir).joinpath(f"{sample}_significant.tsv"), sep="\t", index=False)
    
    idx_sigase = plot_data["shapes"] == "*"
    print("Pearson R all: ", scipy.stats.pearsonr(plot_data["ase"], -plot_data["dmr"]))
    print("Pearson R ASE sig: ", scipy.stats.pearsonr(plot_data["ase"][idx_sigase], -plot_asm[idx_sigase]))
    
    print("Partial corr all: ")
    print(partial_corr(data=df, x="asm_hp1_minus_hp2", y="ase_hp1", covar="cn_hp1"))
    
    print("Partial corr ASE sig: ")
    print(partial_corr(data=df.loc[idx_sigase], x="asm_hp1_minus_hp2", y="ase_hp1", y_covar="cn_hp1"))


if __name__ == "__main__":
    main()
