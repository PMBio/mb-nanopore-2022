import pandas as pd
from nanoepitools.plotting.general_plotting import PlotArchiver

import tqdm
from matplotlib_venn import venn3


from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.config import module_config

"""
Creates plots that show how diminished differential methylation effect in the nearest probes
of an 850k array
"""

if __name__ == "__main__":
    min_diff = 0.5
    
    """
    Read sample diffmet
    """
    
    samplecomp_hits = []
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        samplecomp_hits += [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"], "diff": line["diff"]}
            for line in pm.read_file(b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=min_diff)
        ]
    
    samplecomp_hits = merge_duplicate_diffmet_hits(samplecomp_hits)
    met_850k_manifest_path = module_config.methylation850karray_metadata_path
    manifest_850k = pd.read_csv(met_850k_manifest_path, sep="\t").set_index("probeID", drop=True)
    manifest_850k = manifest_850k[["CpG_chrm", "CpG_beg", "CpG_end"]].copy()
    manifest_850k = manifest_850k.rename({"CpG_chrm": "chrom", "CpG_beg": "start", "CpG_end": "end"}, axis=1)
    manifest_850k["chrom"] = manifest_850k["chrom"].map(lambda x: x.replace("chr", "") if isinstance(x, str) else x)
    manifest_850k = manifest_850k.groupby("chrom")
    
    path_450k_metadata = module_config.methylation450karray_metadata_path
    manifest_450k = pd.read_csv(path_450k_metadata, sep="\t").set_index("probeID", drop=True)
    manifest_450k["CpG_chrm"] = manifest_450k["CpG_chrm"].map(
        lambda x: x.replace("chr", "") if isinstance(x, str) else x
    )
    manifest_450k = manifest_450k.rename({"CpG_chrm": "chrom", "CpG_beg": "start", "CpG_end": "end"}, axis=1)
    manifest_450k = manifest_450k.groupby("chrom")
    
    def load_matching_array(npset, array_manifest):
        array_found = []
        array_missing = []
        for i, hit in tqdm.tqdm(enumerate(npset)):
            subset = array_manifest.get_group(hit["chrom"])
            subset = subset.loc[(subset["start"] <= hit["end"]) & (hit["start"] <= subset["end"])]
            if subset.shape[0] == 0:
                array_missing.append(i)
            else:
                array_found.append(i)
        return set(array_found), set(array_missing)
    
    found_450k, notfound_450k = load_matching_array(samplecomp_hits, manifest_450k)
    found_850k, notfound_850k = load_matching_array(samplecomp_hits, manifest_850k)
    
    pa = PlotArchiver("supplementary_figures", config=module_config)
    
    pa.figure()
    venn3(
        (found_450k, found_850k, found_850k.union(notfound_850k)),
        set_labels=("Illumina 450k", "MethylationEPIC (850k)", "ONT"),
        alpha=0.2,
    )
    pa.savefig("venn_array_discovered_samplediffmet")
