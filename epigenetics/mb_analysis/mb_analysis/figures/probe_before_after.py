from meth5.meth5 import MetH5File
import pandas as pd
from nanoepitools.plotting.general_plotting import PlotArchiver
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import tqdm

from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.config import module_config

"""
Creates plots that show how diminished differential methylation effect in the nearest probes
of an 850k array
"""


def load_matching_array(hits):
    """
    For each DMR, we find if there is a probe in the 850k array, and if not, where the nearest probe before and after are
    """
    array_found = []
    array_missing = []
    nearest_to_missing = []
    for hit in tqdm.tqdm(hits):
        subset = manifest_850k.get_group(hit["chrom"])
        subset = subset.loc[(subset["start"] <= hit["end"]) & (hit["start"] <= subset["end"])]
        if subset.shape[0] == 0:
            array_missing.append(hit)
            subset = manifest_850k.get_group(hit["chrom"])
            subset = subset.loc[subset["end"] < hit["start"]]
            if subset.shape[0] == 0:
                print(hit)
                nearest_before = None
            else:
                dist = np.abs(subset["end"] - hit["start"])
                nearest_before = subset.iloc[np.argmin(dist)].to_dict()
            
            subset = manifest_850k.get_group(hit["chrom"])
            subset = subset.loc[subset["start"] > hit["end"]]
            if subset.shape[0] == 0:
                print(hit)
                nearest_after = None
            else:
                dist = np.abs(subset["start"] - hit["end"])
                nearest_after = subset.iloc[np.argmin(dist)].to_dict()
            nearest_to_missing.append([nearest_before, nearest_after])
        else:
            array_found.append(hit)
    return array_found, array_missing, nearest_to_missing


def compute_three_part_diffmet():
    """ Computes methylation rates from nanopore for (nearest array probe before, DMR, nearest array probe after"""
    diff_compare = []
    for hit, nearest_array in zip(tqdm.tqdm(notfound), notfound_nearest):
        has_np_calls_for_probe = 0
        all_nearest = {}
        for sample in m5_files:
            bs_nearest = []
            for nearest in nearest_array:
                if nearest is not None:
                    rates, _ = (
                        m5_files[sample][hit["chrom"]]
                        .get_values_in_range(nearest["start"] - 5, nearest["end"] + 5)
                        .get_llr_site_rate()
                    )
                    if len(rates) > 0:
                        has_np_calls_for_probe += 1
                        bs = rates.mean()
                        bs_nearest.append(bs)
            all_nearest[sample] = bs_nearest
        if has_np_calls_for_probe == 4:
            # That means we have for primary and relapse measurements for the upstream and downstream probe
            diff_before = all_nearest["Relapse"][0] - all_nearest["Primary"][0]
            diff_after = all_nearest["Relapse"][1] - all_nearest["Primary"][1]
            diff_np = hit["diff"]
        diff_compare.append([diff_before, diff_np, diff_after])
    
    diff_compare = np.array(diff_compare)
    return diff_compare


def plot_met_diff_of_surrounding_probes(diff_compare):
    """Creates a plot showing the methylation rate differences at the location of the nearest 850k array probe to
    demonstrate whether "missed" DMRs by probes would have been detected in the nearest probes"""
    
    dir_correct_before = (diff_compare[:, 1] > 0) == (diff_compare[:, 0] > 0)
    dir_correct_after = (diff_compare[:, 1] > 0) == (diff_compare[:, 2] > 0)
    found_before = (np.abs(diff_compare[:, 0]) > 0.5) & dir_correct_before
    found_after = (np.abs(diff_compare[:, 2]) > 0.5) & dir_correct_after
    
    with pa.open_multipage_pdf("Surrounding 850k calls"):
        pa.figure()
        idx = (~np.isnan(diff_compare[:, 0])) & (~np.isnan(diff_compare[:, 2]))
        plt.violinplot([diff_compare[:, i][idx] for i in range(3)], positions=[0, 1, 2])
        plt.xticks([0, 1, 2], ["Probe before", "Segment", "Probe after"])
        
        for i in range(len(diff_compare)):
            plt.plot([0, 1, 2], diff_compare[i], c="grey", linewidth=0.1, alpha=0.1)
        plt.ylim(-1, 1)
        plt.ylabel("Differential methylation Relapse")
        pa.savefig()


if __name__ == "__main__":
    min_diff = 0.5
    
    pa = PlotArchiver("supplementary_figures", config=module_config)
    
    """
    Read sample diffmet
    """
    samplecomp_hits = []
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        samplecomp_hits += [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"], "diff": line["diff"]}
            for line in pm.read_file(b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=min_diff)
        ]
    
    samplecomp_hits = merge_duplicate_diffmet_hits(samplecomp_hits)
    met_850k_manifest_path = module_config.methylation850karray_metadata_path
    manifest_850k = pd.read_csv(met_850k_manifest_path, sep="\t").set_index("probeID", drop=True)
    manifest_850k = manifest_850k[["CpG_chrm", "CpG_beg", "CpG_end"]].copy()
    manifest_850k = manifest_850k.rename({"CpG_chrm": "chrom", "CpG_beg": "start", "CpG_end": "end"}, axis=1)
    manifest_850k = manifest_850k.groupby("chrom")
    
    m5_files = {sample: module_config.meth5_template_file.format(sample=sample) for sample in ["Primary", "Relapse"]}
    m5_files = {sample: MetH5File(m5_files[sample], mode="r") for sample in m5_files}
    
    found, notfound, notfound_nearest = load_matching_array(samplecomp_hits)
    
    diff_compare = compute_three_part_diffmet()
    
    plot_met_diff_of_surrounding_probes(diff_compare)
