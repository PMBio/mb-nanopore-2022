import tqdm
from pathlib import Path

import numpy as np
import pandas as pd

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


def plot_sv(pa, sv, batched_bam):
    
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
            must_overlap=sv["pos"],
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
        )


if __name__ == "__main__":
    # Output from genotype_mb_alignments.py
    outfile = Path(module_config.sv_vs_met_dir).joinpath("sv_vs_met_diffmet_regions.tsv")
    svs = pd.read_csv(outfile, sep="\t", dtype={"chrom": str})
    
    sample = "Primary"
    bam_dir = Path(module_config.bam_template_dir.format(sample=sample))
    bam_files = [bf for bf in bam_dir.iterdir() if bf.name.endswith(".bam")]
    
    print("Loading GFF", flush=True)
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
    enhancers.load()
    enhancers.annotate_nearest_gene(gff, maxdist=5000)
    enhancers.filter_nearest_gene_none()
    
    pa = PlotArchiver("sv_vs_met", config=module_config)
    pl = Plotter(gff, pa, enhancers)
    with BatchedBamFile(bam_files) as batched_bam:
        for _, sv in tqdm.tqdm(list(svs.iterrows())):
            plot_sv(pa, sv, batched_bam)
