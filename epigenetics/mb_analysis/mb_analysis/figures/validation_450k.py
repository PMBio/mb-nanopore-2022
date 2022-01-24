from pathlib import Path

import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tqdm

from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.plotting.general_plotting import plot_2d_density
from meth5.meth5 import MetH5File


from mb_analysis.config import module_config


def get_sample_correlation(sample, llr_threshold=2.0):
    """
    Reads 450k array data and then matches probe methylation rates with methylation rates from nanopore
    :param sample: Primary or Relapse
    :param llr_threshold:
    :return: generator of tuples (450k methylation rate, nanopore methylation rate)
    """
    met_450k_grouped = met_450k.groupby("CpG_chrm")
    with MetH5File(module_config.meth5_template_file.format(sample=sample), "r") as f, tqdm.tqdm(
        total=met_450k.shape[0]
    ) as pbar:
        for chrom in f.get_chromosomes():
            if not chrom in met_450k_grouped.groups:
                continue
            met_450k_chr = met_450k_grouped.get_group(chrom)
            all_bs, all_ranges = f[chrom].get_all_values().get_llr_site_rate(llr_threshold=llr_threshold)
            for _, row in met_450k_chr.iterrows():
                idx = (all_ranges[:, 1] >= row["CpG_beg"]) & (row["CpG_end"] >= all_ranges[:, 0])
                bs = np.nanmean(all_bs[idx])
                pbar.update(1)
                yield row[sample_dict[sample]], bs


if __name__ == "__main__":
    path_450k_data = Path(module_config.methylation450karray_normed_path)
    path_450k_metadata = Path(module_config.methylation450karray_metadata_path)
    sample_dict = module_config.methylation450karray_sample_dict
    
    met_450k = pd.read_csv(path_450k_data, sep="\t", index_col=0)
    met_450k_manifest = pd.read_csv(path_450k_metadata, sep="\t").set_index("probeID", drop=True)
    
    cpg_coord_mapping = met_450k_manifest[["CpG_chrm", "CpG_beg", "CpG_end"]].copy()
    cpg_coord_mapping["CpG_chrm"] = cpg_coord_mapping["CpG_chrm"].map(
        lambda x: x.replace("chr", "") if isinstance(x, str) else x
    )
    
    met_450k = met_450k.merge(cpg_coord_mapping, left_index=True, right_index=True)
    
    met_450k = met_450k.loc[met_450k["CpG_chrm"].map(lambda x: ~np.isnan(x) if isinstance(x, float) else True)]
    
    sample_xy = {}
    for sample in sample_dict:
        sample_xy[sample] = np.array(list(get_sample_correlation(sample)))
    
    pa = PlotArchiver("supplementary_figures", module_config)
    
    with pa.open_multipage_pdf("np_450k_bs_comparison"):
        for sample in sample_xy:
            idx = ~(np.isnan(sample_xy[sample][:, 1]) | np.isnan(sample_xy[sample][:, 0]))
            pa.figure()
            plt.title(f"Nanopore vs 450k Array methylation rate {sample}")
            plot_2d_density(sample_xy[sample][idx, 1], sample_xy[sample][idx, 0], cmap="jet", contour_colors="w")
            plt.xlabel("Nanopore")
            plt.ylabel("450k Array")
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.axes().set_aspect("equal")
            pa.savefig()
            
            r = scipy.stats.pearsonr(sample_xy[sample][idx, 0], sample_xy[sample][idx, 1])[0]
            print(f"{sample} pearson R:", r)
            print(f"{sample} pearson R2:", r ** 2)
