import tqdm
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.annotations.enhancers import Enhancers
from mb_analysis.summary.plot_gene import Plotter

from mb_analysis.genotyper.genotyper import find_base_in_alignment, read_map_ref_alt_other_alignments, BatchedBamFile

"""
Related to genotype_mb_alignments.py, for a more visual approach we plot the SVs that correlate with differential
methylation.
"""


def plot_sv(pa, pl, sv, batched_bam):
    
    sample_colors = {
        "Primary (ALT)": "r",
        "Primary (REF)": "r",
        "Primary": "r",
        "Germline (ALT)": "g",
        "Germline (REF)": "g",
        "Germline": "g",
        "Relapse (ALT)": "b",
        "Relapse (REF)": "b",
        "Relapse": "b",
    }
    sample_hatch = {
        "Primary (ALT)": "//",
        "Primary (REF)": "\\\\",
        "Primary": "",
        "Germline (ALT)": "//",
        "Germline (REF)": "\\\\",
        "Germline": "",
        "Relapse (ALT)": "//",
        "Relapse (REF)": "\\\\",
        "Relapse": "",
    }
    
    sample_order = [
        "Germline",
        "Germline (ALT)",
        "Germline (REF)",
        "Primary",
        "Primary (ALT)",
        "Primary (REF)",
        "Relapse",
        "Relapse (ALT)",
        "Relapse (REF)",
    ][::-1]
    
    sample_order = [
        "Primary",
        "Primary (ALT)",
        "Primary (REF)",
    ][::-1]

    ######################

    sample_order_combined = [
        "Primary (OTHER)", "Primary (H1) (OTHER)", "Primary (H2) (OTHER)",
        "Primary (ALT)", "Primary (H1) (ALT)","Primary (H2) (ALT)",
        "Primary (REF)", "Primary (H1) (REF)","Primary (H2) (REF)",
    ][::-1]

    sample_colors_combined = {key: "gray" if "OTHER" in key else "r" if "ALT" in key else "b" for key in sample_order_combined}
    sample_hatch_combined = {key: "//" if "H1" in key else "\\\\" if "H2" in key else "" for key in sample_order_combined}
    ######################
    figure_kwargs = {"figsize": (20, 7)}
    key = f"{sv['chrom']}_{sv['pos']}"
    with pa.open_multipage_pdf(f"diffmet_sv_with_annotations_{key}"):
        read_mapping = read_map_ref_alt_other_alignments(batched_bam, sv["chrom"], sv["pos"], sv["ref"], sv["alt"])
        
        def altref_sample_fun(sample, met_matrix):
            return np.array(
                [
                    f"{sample} ({read_mapping[rn].upper()})" if rn in read_mapping else sample
                    for rn in met_matrix.read_names
                ]
            )
        
        start = min([sv["pos"], sv["diffmet_start"], sv["diffmet_end"]]) - 5000
        end = max([sv["pos"], sv["diffmet_start"], sv["diffmet_end"]]) + 5000
        annotations = [
            {"region": [sv["diffmet_start"], sv["diffmet_end"]], "text": "Diffmet", "color": "r"},
            {"region": [sv["pos"] - 10, sv["pos"] + 10], "text": "SNP", "color": "b"},
        ]
        pl.plot_region_custom_samples(
            sv["chrom"],
            start,
            end,
            altref_sample_fun,
            sample_colors,
            sample_hatch,
            sample_order,
            title=f"{sv['chrom']}:{sv['pos']} Ref:{sv['ref']} Alt:{sv['alt']}",
            annotations=annotations,
            ws=0,
            figure_kwargs=figure_kwargs,
            coordinate_space=False,
            must_overlap=sv["pos"],
            requires_read_names_for_sample_name_fun=True,
        )

        def altref_hp_sample_fun(sample, met_matrix):
            haplotypes = [f" ({rg})" if rg in ["H1", "H2"] else "" for rg in met_matrix.read_samples]
            gts = [f" ({read_mapping.get(rn, 'OTHER').upper()})" for rn in met_matrix.read_names]
            return np.array([f"{sample}{hp}{gt}" for hp, gt in zip(haplotypes, gts)])

        pl.plot_region_custom_samples(
            sv["chrom"],
            start,
            end,
            altref_hp_sample_fun,
            sample_colors_combined,
            sample_hatch_combined,
            sample_order_combined,
            title=f"{sv['chrom']}:{sv['pos']} Ref:{sv['ref']} Alt:{sv['alt']}",
            annotations=annotations,
            ws=0,
            figure_kwargs=figure_kwargs,
            coordinate_space=False,
            must_overlap=sv["pos"],
            requires_read_names_for_sample_name_fun=True,
        )
        
        pl.plot_region(
            sv["chrom"],
            start,
            end,
            title=f"{sv['chrom']}:{sv['pos']} Ref:{sv['ref']} Alt:{sv['alt']}",
            annotations=annotations,
            ws=0,
            figure_kwargs=figure_kwargs,
            must_overlap=sv["pos"],
            with_relapse=False,
            coordinate_space=False,
            with_germline=False,
        )


if __name__ == "__main__":
    # Output from genotype_mb_alignments.py
    outfile = Path(module_config.sv_vs_met_dir).joinpath("sv_vs_met_diffmet_regions.tsv")
    svs = pd.read_csv(outfile, sep="\t", dtype={"chrom": str})
    
    sample = "Primary"
    bam_dir = Path(module_config.bam_template_dir.format(sample=sample))
    bam_files = [bf for bf in bam_dir.iterdir() if bf.name.endswith("sorted.filtered.bam")]
    
    print("Loading GFF", flush=True)
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
    enhancers.load(replace_chr=False)
    enhancers.annotate_nearest_gene(gff, maxdist=5000)
    enhancers.filter_nearest_gene_none()

pa = PlotArchiver("sv_vs_met", config=module_config)
pl = Plotter(gff, pa, enhancers)

def compute_diffmet(llrs, samples, sample_a, sample_b, extra_index):
    idx_a = (samples == sample_a) & extra_index
    idx_b = (samples == sample_b) & extra_index
    llrs_a = np.array(llrs[idx_a, :].todense()).flatten()
    llrs_b = np.array(llrs[idx_b, :].todense()).flatten()
    llrs_a = llrs_a[llrs_a != 0]
    llrs_b = llrs_b[llrs_b != 0]
    if len(llrs_a) < 3 or len(llrs_b) < 3:
        return np.nan, np.nan, np.nan
    
    bs_a = (llrs_a > 2).sum() / (np.abs(llrs_a) > 2).sum()
    bs_b = (llrs_b > 2).sum() / (np.abs(llrs_b) > 2).sum()
    
    stat, pval = scipy.stats.mannwhitneyu(llrs_a, llrs_b, alternative="two-sided")
    return np.nanmean(bs_a) - np.nanmean(bs_b), stat, pval


with BatchedBamFile(bam_files) as batched_bam:
    for i, sv in tqdm.tqdm(list(svs.iterrows())):
        ref_alt_mapping = read_map_ref_alt_other_alignments(batched_bam, sv["chrom"], sv["pos"], sv["ref"], sv["alt"])
        met_matrix = pl.loader.get_merged_matrix(
            sv["chrom"],
            sv["diffmet_start"] - 2,
            sv["diffmet_end"] + 2,
            with_relapse=False,
            with_germline=False,
            requires_read_names_for_sample_name_fun=True,
        )
        ref_alt_samples = np.array([ref_alt_mapping.get(r, "Other") for r in met_matrix.read_names])
        
        reads_allowed = [
            read
            for read, hp, refalt in zip(met_matrix.read_names, met_matrix.read_samples, ref_alt_samples)
            if hp != "Primary" and refalt != "Other"
        ]
        extra_index = np.array([r in reads_allowed for r in met_matrix.read_names])
        
        diffmet_haplotype, hp_stat, hp_pval = compute_diffmet(
            met_matrix.met_matrix, met_matrix.read_samples, "Primary (HP1)", "Primary (HP2)", extra_index
        )
        diffmet_sv, sv_stat, sv_pval = compute_diffmet(
            met_matrix.met_matrix, ref_alt_samples, "ref", "alt", extra_index
        )
        
        if np.isnan(diffmet_sv):
            continue
        
        if np.isnan(diffmet_haplotype):
            print("Warning: cant compute for haplotypes")
            continue
        
        if np.abs(diffmet_sv) <= np.abs(diffmet_haplotype) + 0.1:
            continue
        
        if hp_pval < sv_pval:
            continue
        print("############")
        print(sv["chrom"], sv["pos"])
        print(diffmet_sv, sv_stat, sv_pval)
        print(diffmet_haplotype, hp_stat, hp_pval)
        plot_sv(pa, pl, sv, batched_bam)
