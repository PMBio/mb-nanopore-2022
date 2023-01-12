import numpy as np
import scipy.stats

from pingouin import partial_corr
import pandas as pd
import matplotlib.pyplot as plt
from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.summary.summarize import get_diffmet_promoter_primary_relapse, get_diffexp
from mb_analysis.ase_asm_analysis.cnv import load_tumor_sample_cnv

"""
Plotting differential expression versus differential methylation for supplementary figure
"""

if __name__ == "__main__":
    pa = PlotArchiver("supplementary_figures", config=module_config)
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    diffmet = get_diffmet_promoter_primary_relapse(gff, 2000, 500)
    diffexp = get_diffexp()
    
    genes = list(set(diffmet.keys()).intersection(diffexp.keys()))
    diffmet_intersect = np.array([diffmet[k] for k in genes])
    diffexp_intersect = np.array([diffexp[k] for k in genes])
    
    r, p = scipy.stats.pearsonr(diffexp_intersect, diffmet_intersect)
    regression_slope, regression_offset = np.polyfit(diffmet_intersect, diffexp_intersect, 1)
    
    pa.figure()
    plt.xlim(-1,1)
    
    plot_2d_density(diffmet_intersect, diffexp_intersect, cmap="jet")
    x = np.arange(-1, 1, 0.01)
    plt.plot([plt.xlim()[0], plt.xlim()[1]], [0,0], c="k")
    plt.title("Differential methylation versus differential expression (Relapse - Primary)")
    plt.xlabel("Methylation rate difference")
    plt.ylabel("Log fold change expression")
    pa.savefig("diffmet_diffexp_correlation")
    
    # Highly differentially expressed
    highly_differentially_expressed = diffexp.loc[abs(diffexp) > 2]
    highly_expressed_with_diffmet = list(set(genes).intersection(set(highly_differentially_expressed.index)))
    diffmet_hewdm = np.array([diffmet[k] for k in highly_expressed_with_diffmet])
    diffexp_hewdm = np.array([diffexp[k] for k in highly_expressed_with_diffmet])
    print("Spearman-R for highly differentially expressed versus promoter methylation")
    print(scipy.stats.spearmanr(diffmet_hewdm, diffexp_hewdm))

    pa.figure()
    plt.xlim(-1, 1)
    plot_2d_density(diffmet_intersect, diffexp_intersect, cmap="jet")
    plt.scatter(diffmet_hewdm, diffexp_hewdm, marker="x", color="white")
    plt.title("Differential methylation versus differential expression (Relapse - Primary)")
    plt.xlabel("Methylation rate difference")
    plt.ylabel("Log fold change expression")
    pa.savefig("diffmet_diffexp_correlation_high_diffexp")

    print("Fisher exact test for enrichment of Diffmet/Diffex overlap (highly expressed):")
    n_diffexp_dmr = len(diffexp_hewdm)
    n_diffexp_notdmr = len(highly_differentially_expressed) - n_diffexp_dmr
    n_notdiffexp_dmr = len(diffmet) - n_diffexp_dmr
    num_genes = sum(1  for chrom in module_config.chroms if chrom in gff.chromosomes for gene in gff.chromosomes[chrom].children )
    num_genes = len(diffexp)
    n_notdiffexp_notdmr = num_genes - n_notdiffexp_dmr - n_diffexp_notdmr + n_diffexp_dmr
    print(scipy.stats.fisher_exact([[n_diffexp_dmr, n_diffexp_notdmr], [n_notdiffexp_dmr, n_notdiffexp_notdmr]], alternative="two-sided"))
    
    gene_copy_number = load_tumor_sample_cnv(module_config.tumor_coverage_file, gff)
    cn_diff = gene_copy_number["Relapse"] - gene_copy_number["Primary"]
    cn_diff = cn_diff.loc[~np.isinf(cn_diff) & ~np.isnan(cn_diff)]
    
    genes = list(set(diffexp.keys().intersection(set(cn_diff.keys()))))
    diffexp_intersect_cn = np.array([diffexp[k] for k in genes])
    cn_intersect_diffexp = np.array([cn_diff.loc[k] for k in genes])
    
    highly_diffcn = cn_diff.loc[np.abs(cn_diff) >= 1]
    highly_expressed_with_cn = list(set(highly_diffcn.index).intersection(set(highly_differentially_expressed.index)))
    diffcn_hewcn = np.array([cn_diff.loc[k] for k in highly_expressed_with_cn])
    diffexp_hewcn = np.array([diffexp[k] for k in highly_expressed_with_cn])
    print("Spearman-R for highly differentially expressed versus copy number variation")
    print(scipy.stats.spearmanr(diffcn_hewcn, diffexp_hewcn))
    
    print("Intersection between copy-number related diffexp and DMR related: ", len(set(highly_expressed_with_diffmet).intersection(set(highly_expressed_with_cn))))

    def plot_2d_density(x, y, nbins=50, cmap=plt.cm.BuGn_r, contour_colors="k"):
        k = scipy.stats.gaussian_kde((x, y))
        xi, yi = np.mgrid[x.min(): x.max(): nbins * 1j, y.min(): y.max(): nbins * 1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        print(min(zi), max(zi))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading="gouraud", cmap=cmap)
        plt.colorbar()
        plt.contour(xi, yi, zi.reshape(xi.shape), colors=contour_colors)



    pa.figure()
    plt.xlim(-1, 1)
    plot_2d_density(diffmet_intersect, diffexp_intersect, cmap="jet")
    diffcn_hewdm_also_cn = np.array([cn_diff.loc[k] if k in cn_diff.index else 0 for k in highly_expressed_with_diffmet])
    plt.scatter(diffmet_hewdm, diffexp_hewdm, marker="x", color="white")
    idx_cn_low = diffcn_hewdm_also_cn < -1
    idx_cn_high= diffcn_hewdm_also_cn > 1
    plt.scatter(diffmet_hewdm[idx_cn_low], diffexp_hewdm[idx_cn_low], marker="x", color="red")
    plt.scatter(diffmet_hewdm[idx_cn_high], diffexp_hewdm[idx_cn_high], marker="x", color="yellow")
    plt.title("Differential methylation versus differential expression (Relapse - Primary)")
    plt.xlabel("Methylation rate difference")
    plt.ylabel("Log fold change expression")
    pa.savefig("diffmet_diffexp_correlation_high_diffexp")


    df = pd.DataFrame({"cn":diffcn_hewdm_also_cn, "exp":diffexp_hewdm, "met":diffmet_hewdm })

    print(partial_corr(data=df, x="met", y="exp", covar="cn", method="spearman"))
    scipy.stats.spearmanr(diffmet_hewdm, diffexp_hewdm)
    
    
    genes = list(set(diffmet.keys()).intersection(set(cn_diff.keys())))
    