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
from mb_analysis.summary.summarize import get_diffmet_promoter_primary_relapse, get_diffexp, load_gene_ase_primary
from mb_analysis.summary.outlier_analysis_result import OutlierAnalysisResults
from mb_analysis.config import module_config

"""
Creates ASEs ASM figure. At the time of writing this is Figure 4b.
"""


def compute_ase_asm(ase_gene_ratio_df, asm_result, gene_cnv_dict, gff, min_abs_diff, sample, only_protein_coding=False):
    promoters_hit = asm_result[sample]
    ase_sample = ase_gene_ratio_df.loc[ase_gene_ratio_df.index.map(lambda x: x in promoters_hit)].copy()
    #ase_sample["fdr"] = fdr_from_pvals(ase_sample["pval"])
    plot_data = dict(ase=[], asm=[], cnv=[], color=[], gene_ids=[], labels=[], shapes=[], is_pc=[], pval=[])
    
    for gene, row in ase_sample.iterrows():
        gene_feature = gff.get_gene(gene)
        protein_coding = gene_feature.info["biotype"] == "protein_coding"
        if only_protein_coding and not protein_coding:
            continue
        gene_name = gene_feature.name
        # If there are multiple effects in the promoter, pick the region with the largest effect
        diffmet = promoters_hit[gene]
        if abs(diffmet) < min_abs_diff:
            continue
        
        is_sig = row["pval"] < 0.05
        plot_data["pval"].append(row["pval"])
        ase_ratio = row["hp1_ratio"]
        plot_data["ase"].append(ase_ratio)
        plot_data["shapes"].append("*" if is_sig else ("o" if protein_coding else "x"))
        
        plot_data["asm"].append(diffmet)
        plot_data["labels"].append(gene_name)
        plot_data["gene_ids"].append(gene)
        plot_data["is_pc"].append(protein_coding)
        
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


if __name__ == "__main__":
    pa = PlotArchiver("mb_asm_ase", config=module_config)
    gene_to_variant_mapping_mode = "nearest"
    samples = ["Primary"]
    min_abs_diff = 0.25
    
    """ load data """
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False, keep_info=["biotype"])
    gff.build_index()
    
    with CollectiveAlleleSpecificMethylation(gff, samples=samples) as asm:
         asm_result= asm.read_all_samples(annotation="promoters", min_diff=min_abs_diff)
    
    vcf_blood = PhasedVCF(module_config.vcf_blood_file)
    vcf_blood.read()
    
    gene_cnv_dict = {
        sample: load_gene_cnv(module_config.cnv_file_template.format(sample=sample), vcf_blood, gff, replace_chr=False)
        for sample in samples
    }
    
    ase_gene_ratio_df = load_gene_ase_primary(vcf_blood, pval_thres=None)[1]
    
    """ Compute ASE/ASE/ASCNV for each sample separately """
    for sample in samples:
        plot_data = compute_ase_asm(ase_gene_ratio_df, asm_result, gene_cnv_dict, gff, min_abs_diff, sample)
        
        df = pd.DataFrame(
            {
                "gene_ids": plot_data["gene_ids"],
                "gene_name": plot_data["labels"],
                "ase_hp1": plot_data["ase"],
                "asm_hp1_minus_hp2": -plot_data["asm"],
                "cn_hp1": plot_data["cnv"],
            }
        )
        df.to_csv(Path(module_config.asm_vs_ase_dir).joinpath(f"{sample}_significant.tsv"), sep="\t", index=False)
        
        with pa.open_multipage_pdf(f"ase_vs_asm_promoter_{sample}"):
            # Don't highlight significance
            plot_ase_vs_asm_nosplit(
                pa,
                plot_data["ase"],
                plot_data["asm"],
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
            
            # Only significant and protein coding
            labels_printed = [l if s == "*" and pc else "" for l, s, pc in zip(plot_data["labels"], plot_data["shapes"], plot_data["is_pc"])]
            # Highlight significance
            plot_ase_vs_asm_nosplit(
                pa,
                plot_data["ase"],
                plot_data["asm"],
                labels_printed,
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
                "ase_hp1": plot_data["ase"], "ase_hp1_pval": plot_data["pval"],
                "asm_hp1_minus_hp2": -plot_data["asm"],
                "cn_hp1": plot_data["cnv"],
            }
        )
        df.to_csv(Path(module_config.asm_vs_ase_dir).joinpath(f"{sample}_significant.tsv"), sep="\t", index=False)
    

    # Compute statistics cited in paper
    print("Fisher exact test for enrichment of ASM/ASE overlap:")
    all_expr = OutlierAnalysisResults()
    all_expr.load(filter=False)

    idx_sigase = plot_data["shapes"] == "*"

    print("thresholding asm to min diff 0.5 and requiring ase to be significant")
    idx = idx_sigase & (np.abs(plot_data["asm"]) >= 0.5)
    print("Pearson R: ", scipy.stats.pearsonr(plot_data["ase"][idx], -plot_data["asm"][idx]))
    print("Spearman R: ", scipy.stats.spearmanr(plot_data["ase"][idx], -plot_data["asm"][idx]))

    print("Partial corr: ")
    print(partial_corr(data=df.loc[idx], x="asm_hp1_minus_hp2", y="ase_hp1", covar="cn_hp1", method="pearson"))
    print(partial_corr(data=df.loc[idx], x="asm_hp1_minus_hp2", y="ase_hp1", covar="cn_hp1", method="spearman"))

    print("Partial corr: ")
    print(partial_corr(data=df.loc[idx], x="asm_hp1_minus_hp2", y="ase_hp1", covar="cn_hp1", method="pearson"))
    print(partial_corr(data=df.loc[idx], x="asm_hp1_minus_hp2", y="ase_hp1", covar="cn_hp1", method="spearman"))

    num_genes = (all_expr.get_mb_expression("Primary") > 0).sum()  # use all expressed genes as background for test
    n_ase_dmr = idx.sum()
    n_ase_notdmr = (ase_gene_ratio_df["pval"] < 0.05).sum() - n_ase_dmr
    n_notase_dmr = sum(abs(diff) > 0.5 for diff in asm_result["Primary"].values()) - n_ase_dmr
    n_notase_notdmr = num_genes - n_notase_dmr - n_ase_notdmr + n_ase_dmr
    print(scipy.stats.fisher_exact([[n_ase_dmr, n_ase_notdmr], [n_notase_dmr, n_notase_notdmr]], alternative="two-sided"))

    df_asc = df.loc[idx & ~df["cn_hp1"].isna()]
    print("Pearson R: ", scipy.stats.pearsonr(df_asc["asm_hp1_minus_hp2"], df_asc["cn_hp1"]))
    num_pos = ((df_asc["cn_hp1"] > 0.65) & (df_asc["asm_hp1_minus_hp2"] > 0.5)).sum()
    num_neg = ((df_asc["cn_hp1"] < 0.35) & (df_asc["asm_hp1_minus_hp2"] < 0.5)).sum()



    print("without thresholding for figure legend")
    idx = idx_sigase
    print("Pearson R: ", scipy.stats.pearsonr(plot_data["ase"][idx], -plot_data["asm"][idx]))
    print("Spearman R: ", scipy.stats.spearmanr(plot_data["ase"][idx], -plot_data["asm"][idx]))
    num_genes = (all_expr.get_mb_expression("Primary") > 0).sum()  # use all expressed genes as background for test
    n_ase_dmr = idx.sum()
    n_ase_notdmr = (ase_gene_ratio_df["pval"] < 0.05).sum() - n_ase_dmr
    n_notase_dmr = ase_gene_ratio_df.shape[0] - n_ase_dmr
    n_notase_notdmr = num_genes - n_notase_dmr - n_ase_notdmr + n_ase_dmr
    print(
        scipy.stats.fisher_exact([[n_ase_dmr, n_ase_notdmr], [n_notase_dmr, n_notase_notdmr]], alternative="two-sided")
    )
    
    df_asc = df.loc[idx & ~df["cn_hp1"].isna()]
    print("Pearson R: ", scipy.stats.pearsonr(df_asc["asm_hp1_minus_hp2"], df_asc["cn_hp1"]))
    num_pos = ((df_asc["cn_hp1"] > 0.65) & (df_asc["asm_hp1_minus_hp2"] > 0.5)).sum()
    num_neg = ((df_asc["cn_hp1"] < 0.35) & (df_asc["asm_hp1_minus_hp2"] < 0.5)).sum()
