import tqdm

from nanoepitools.plotting.general_plotting import PlotArchiver
import matplotlib.pyplot as plt
import numpy as np

from mb_analysis.config import module_config
from mb_analysis.reference_cpgs import ReferenceCpGs, format_cpg
from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from nanoepitools.annotations.annotations import GFFAnnotationsReader


"""
Creates plots that show how allele specific methylation could not have been detected if it were
simply done with short read bisulfite sequencing, due to distance from SNVs
"""


def load_asm_hits(samples, min_diff):
    asm_files = {
        sample: module_config.pycometh_haplotype_sample_template_file.format(sample=sample)
        for sample in module_config.samples
    }
    asm_hits = {}
    for sample in samples:
        asm_file = asm_files[sample]
        pm_asm = PycomethOutput(met_comp_file=asm_file)
        hits = [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"], "diff": line["diff"]}
            for line in pm_asm.read_file(
                b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=min_diff
            )
        ]
        hits = merge_duplicate_diffmet_hits(hits)
        asm_hits[sample] = hits
    return asm_hits


def annotate_asm_hits(asm_hits, vcf_df_grouped, ref_cpgs):
    new_asm_hits = {}
    for sample, hits in asm_hits.items():
        # Check distance to nearest SNV
        new_hits = []
        for hit in tqdm.tqdm(hits):
            hit["CpGs_pos"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=True, formatted=False)
            hit["complex_CpGs_pos"] = ref_cpgs.get_CGs(
                hit["chrom"], hit["start"], hit["end"], upper=False, formatted=False
            )
            hit["CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=True, formatted=True)
            hit["complex_CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=False, formatted=True)
            
            vcf_df_chrom = vcf_df_grouped.get_group(hit["chrom"])
            hit["CpGs_in_snv_range"] = []
            for cg in hit["CpGs_pos"]:
                dists = np.abs(vcf_df_chrom["pos"] - cg)
                if any(dists < 150):
                    hit["CpGs_in_snv_range"].append(format_cpg(hit["chrom"], cg))
            
            hit["complex_CpGs_in_snv_range"] = []
            for cg in hit["complex_CpGs_pos"]:
                dists = np.abs(vcf_df_chrom["pos"] - cg)
                if any(dists < 150):
                    hit["complex_CpGs_in_snv_range"].append(format_cpg(hit["chrom"], cg))
            
            new_hits.append(hit)
        new_asm_hits[sample] = new_hits
    return new_asm_hits



if __name__ == "__main__":
    pa = PlotArchiver("summary", config=module_config)
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    ref_cpgs = ReferenceCpGs()
    
    min_diff = 0.5
    
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
        lambda x: (x["ALT"] == "T" and x["REF"] == "C") or (x["ALT"] == "C" and x["REF"] == "T"), axis=1
    )
    
    useable_variants = vcf_df.loc[~bs_masked_variants]
    vcf_df_grouped = vcf_df.groupby("chr")
    
    """ Load ASM"""
    sample_order = ("Primary", "Relapse")
    asm_hits = load_asm_hits(sample_order, min_diff)
    asm_hits = annotate_asm_hits(asm_hits, vcf_df_grouped, ref_cpgs)
    
    with pa.open_multipage_pdf("snv_dist_to_asm"):
        num_bps_sample_diffmet = len({cpg for hit in samplecomp_hits for cpg in hit["CpGs"]})
        num_bps_sample_diffmet_highcomplex = len({cpg for hit in samplecomp_hits for cpg in hit["complex_CpGs"]})
        num_bps_sample_diffmet_lowcomplex = num_bps_sample_diffmet - num_bps_sample_diffmet_highcomplex
        
        # Could be phased in BS-seq
        bar_bps_asm_can_find_bs = np.array(
            [len({cpg for hit in asm_hits[s] for cpg in set(hit["complex_CpGs_in_snv_range"])}) for s in sample_order]
        )
        # Could be phased in BS-seq but is repeat-masked
        bar_bps_asm_low_complex = np.array(
            [
                len(
                    {
                        cpg
                        for hit in asm_hits[s]
                        for cpg in set(hit["CpGs_in_snv_range"]).difference(set(hit["complex_CpGs_in_snv_range"]))
                    }
                )
                for s in sample_order
            ]
        )
        # Could not be phased in BS-seq
        bar_bps_asm_unphaseable = np.array(
            [
                len(
                    {
                        cpg
                        for hit in asm_hits[s]
                        for cpg in set(hit["complex_CpGs"]).difference(set(hit["complex_CpGs_in_snv_range"]))
                    }
                )
                for s in sample_order
            ]
        )
        # Could not be phased in BS-seq and is repeat-masked
        bar_bps_asm_low_complex_unphaseable = np.array(
            [
                len(
                    {
                        cpg
                        for hit in asm_hits[s]
                        for cpg in set(hit["CpGs"])
                        .difference(set(hit["complex_CpGs"]))
                        .difference(set(hit["CpGs_in_snv_range"]))
                    }
                )
                for s in sample_order
            ]
        )
        
        """ Plotting """
        plt.figure(figsize=(8, 3))
        y = [num_bps_sample_diffmet_highcomplex, *bar_bps_asm_can_find_bs]
        plt.bar([1, 2.5, 3.5], y, color="#603E95", label="Detectable in BS-seq")
        bottom = np.array(y)
        y = [num_bps_sample_diffmet_lowcomplex, *bar_bps_asm_low_complex]
        plt.bar([1, 2.5, 3.5], y, bottom=bottom, color="#009DA1", label="Low complexity")
        bottom += np.array(y)
        y = [0, *bar_bps_asm_low_complex_unphaseable]
        plt.bar([1, 2.5, 3.5], y, bottom=bottom, color="#D7255D", label="Low complexity & Unphaseable")
        bottom += np.array(y)
        y = [0, *bar_bps_asm_unphaseable]
        plt.bar([1, 2.5, 3.5], y, bottom=bottom, color="#FAC22B", label="Unphaseable")
        plt.xticks([1, 2.5, 3.5], ["Primary vs Relapse", "ASM Primary", "ASM Relapse"])
        plt.ylim(0, plt.ylim()[1] * 1.1)
        plt.ylabel("CpGs in regions")
        plt.title(f"Number of CpGs in differentially methylated regions (difference >= {min_diff})")
        plt.legend()
        pa.savefig()
    
    with pa.open_multipage_pdf("diffmet_segment_length"):
        lengths = np.log10([h["end"] - h["start"] for h in samplecomp_hits])
        pa.figure()
        plt.hist(lengths, bins=200)
        plt.xlabel("Log 10 number of basepairs")
        plt.ylabel("Frequency")
        plt.xlim(0, 5)
        plt.title(f"Primary vs Relapse diffmet segment length")
        pa.savefig()
        
        lengths = np.log10([h["end"] - h["start"] for h in asm_hits["Primary"]])
        pa.figure()
        plt.hist(lengths, bins=200)
        plt.xlabel("Log 10 number of basepairs")
        plt.ylabel("Frequency")
        plt.xlim(0, 5)
        plt.title(f"Primary ASM segment length")
        pa.savefig()
