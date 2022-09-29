import re

import itertools
from typing import Dict, Optional
from types import FunctionType
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm

from meth5.meth5 import MetH5File, ChromosomeContainer, MethlyationValuesContainer, compute_betascore

from mb_analysis.config import module_config
from mb_analysis.ase_asm_analysis.collective_asm import CollectiveAlleleSpecificMethylation
from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature, GeneNameToEnsemblID
from mb_analysis.ase_asm_analysis.cnv import load_gene_cnv, load_tumor_sample_cnv
from nanoepitools.pycometh_result import PycomethOutput
from nanoepitools.annotations.enhancers import Enhancers
from mb_analysis.summary.plot_gene import Plotter
from mb_analysis.summary.fusion_genes import FusionGenes, FusionPartner
from mb_analysis.summary.double_minute_parts import DoubleMinuteParts
from mb_analysis.summary.outlier_analysis_result import OutlierAnalysisResults
from mb_analysis.chromothriptic_breakpoints.breakpoints import ChromothripticBreakpoints
from mb_analysis.ase_asm_analysis.ase import ASE
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from nanoepitools.plotting.general_plotting import PlotArchiver

"""This module attempts to create a grand summary of genes affected by various affects across all samples"""


def isnotnull(x):
    if isinstance(x, list):
        return True
    return not pd.isnull(x)


def pretty_print_iterable(x):
    if isinstance(x, float):
        if pd.isnull(x):
            return ""
    return ", ".join([format_float(xi) if isinstance(xi, float) else str(xi) for xi in x])


def format_float(x):
    return f"{x:.4f}"


def pretty_print_default(x):
    if isnotnull(x):
        if isinstance(x, float):
            return format_float(x)
        else:
            return str(x)
    else:
        return ""


class SummaryTable:
    def __init__(self, gff):
        self.summary = pd.DataFrame(
            index={
                gene.sanitized_id() for chrom in gff.chromosomes for gene in gff.chromosomes[chrom].children.values()
            }
        )
        self.summary["gene_name"] = pd.Series(
            {
                gene.sanitized_id(): gene.name
                for chrom in gff.chromosomes
                for gene in gff.chromosomes[chrom].children.values()
            }
        )
        self.column_importance_fun = {}
        self.column_pretty_print_fun = {"gene_name": lambda x: x}
        self.combined_importance_functions = {}
        self.loc = self.summary.loc
    
    def add_column(
        self,
        colname: str,
        rows: Dict,
        dtype: type,
        importance: Optional[float] = None,
        importance_fun: Optional[FunctionType] = None,
        pretty_print_fun: Optional[FunctionType] = None,
    ):
        self.summary[colname] = pd.Series(rows, dtype=dtype)
        if importance is not None:
            importance_fun = lambda x: (1 if isnotnull(x) else 0) * importance
        if importance_fun is not None:
            self.column_importance_fun[colname] = importance_fun
        if pretty_print_fun is None:
            pretty_print_fun = pretty_print_default
        self.column_pretty_print_fun[colname] = pretty_print_fun
    
    def compute_importance(self):
        importance = pd.Series(0, index=self.summary.index, dtype=float)
        for colname, importance_fun in self.column_importance_fun.items():
            importance += self.summary[colname].apply(importance_fun)
        
        for importance_fun in self.combined_importance_functions.values():
            importance += self.summary.apply(importance_fun, axis=1)
        
        self.column_pretty_print_fun["importance"] = lambda x: f"{x:.2f}"
        self.summary["importance"] = importance
    
    def get_importance_filtered(self, min_importance: float = 1.0):
        self.compute_importance()
        return self.summary.loc[self.summary["importance"] >= min_importance].sort_values("importance", ascending=False)
    
    def to_csv(self, filename: str, min_importance: float = 1.0, sep: str = "\t"):
        if min_importance > 0:
            self.compute_importance()
            df = self.summary.loc[self.summary["importance"] >= min_importance].sort_values(
                "importance", ascending=False
            )
        else:
            df = self.summary
        
        with open(filename, "w") as f:
            columns = list(df.columns)
            f.write("EnsemblID")
            f.write(sep)
            f.write(sep.join(columns))
            
            for geneid, row in df.iterrows():
                f.write("\n")
                f.write(geneid)
                for column in columns:
                    f.write(sep)
                    pretty_print_fun = self.column_pretty_print_fun[column]
                    f.write(str(pretty_print_fun(row[column])))
    
    def add_combined_importance_function(self, key, fun):
        self.combined_importance_functions[key] = fun
    
    def __getitem__(self, item):
        return self.summary[item]
    
    def __repr__(self):
        return repr(self.get_importance_filtered())


def pretty_print_int64(x):
    return "" if pd.isnull(x) else x


def get_diffmet_promoter_primary_relapse(gff, before_tss, after_tss, aggregate=True):
    rows = {}
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        promoters_hit = pm.load_promoters_hit(
            gff, before_tss, after_tss, b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=0.5
        )
        for geneid, new_entries in promoters_hit.items():
            if geneid in rows:
                rows[geneid] += new_entries
            else:
                rows[geneid] = new_entries
    if aggregate:
        rows = {k: np.nanmean([hit["diffmet"] for hit in v]) for k, v in tqdm.tqdm(list(rows.items()))}
    return rows


def add_diffmet_promoter_primary_relapse(summary, gff):
    print("=== Differential promoter methylation Primary vs Relapse ===")
    rows = get_diffmet_promoter_primary_relapse(gff, before_tss=2000, after_tss=500)
    summary.add_column(colname="Promoters diffmet (Relapse-Primary)", rows=rows, dtype=float, importance=2)
    rows = get_diffmet_promoter_primary_relapse(gff, before_tss=5000, after_tss=5000)
    summary.add_column(colname="5k from TSS diffmet (Relapse-Primary)", rows=rows, dtype=float, importance=1)


def add_diffmet_enhancer_primary_relapse(summary, gff):
    print("=== Differential enhancer methylation Primary vs Relapse ===")
    enhancers = Enhancers(module_config.enhancer_cerebellum_file)
    enhancers.load(replace_chr=False)
    enhancers.annotate_nearest_gene(gff, maxdist=3e4)
    enhancers.filter_nearest_gene_none()
    
    rows = {}
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        enhancers_hit = pm.load_enhancers_hit(
            enhancers, b_minus_a=True, drop_insignificant=False, pval_threshold=0.1, min_diff=0.5
        )
        rows.update({k: np.nanmean([hit["diffmet"] for hit in v]) for k, v in enhancers_hit.items()})
    summary.add_column(colname="Enhancer 30k from TSS diffmet (Relapse-Primary)", rows=rows, dtype=float, importance=2)


class HaplotypeImportanceFunction:
    def __init__(self, sample):
        self.weight = 0 if sample == "Germline" else 1 if sample == "Relapse" else 2
    
    def __call__(self, x):
        has_a_hit = 1 if isnotnull(x) else 0
        return has_a_hit * self.weight


def add_diffmet_haplotype(summary, gff):
    print("=== Haplotype differential promoter and enhancer methylation in each sample ===")
    
    with CollectiveAlleleSpecificMethylation(gff) as asm:
        asm_result = asm.read_all_samples(annotation="promoters")
        
        for sample in asm_result.keys():
            summary.add_column(
                colname=f"Promoters Diffmet Haplotype {sample} (HP 2-1)",
                rows=asm_result[sample],
                dtype=float,
                importance_fun=HaplotypeImportanceFunction(sample),
            )
        
        asm_result = asm.read_all_samples(annotation="enhancers")
        for sample in asm_result.keys():
            summary.add_column(
                colname=f"Enhancers Diffmet Haplotype {sample} (HP 2-1)",
                rows=asm_result[sample],
                dtype=float,
                importance_fun=HaplotypeImportanceFunction(sample),
            )


def get_diffexp():
    expr = OutlierAnalysisResults()
    expr.load()
    diff_exp = expr.get_mb_expression("Relapse") - expr.get_mb_expression("Primary")
    return diff_exp


def add_diffexp(summary, gff):
    print("=== Differential expression ===")
    diff_exp = get_diffexp()
    summary.add_column(
        colname="Log Differential expression (Relapse - Primary)",
        rows=diff_exp.to_dict(),
        dtype=float,
        importance_fun=lambda x: np.abs(x) > 2,
    )


def add_gene_fusion(summary, gff):
    print("=== Part of gene fusion ===")
    for colname, fusion_file, importance in (
        ("Fusion targets Primary", module_config.fusion_genes_primary_file, 1),
        ("Fusion targets Relapse", module_config.fusion_genes_relapse_file, 0.25),
    ):
        fusions = FusionGenes(gff).load(fusion_file).filter(mode="high_confidence")
        rows = {}
        for gene_a, gene_b in fusions:
            confidence = gene_b.confidence
            if gene_a not in rows:
                rows[gene_a] = set()
            rows[gene_a].update({gene_b})
        summary.add_column(
            colname=colname, rows=rows, dtype="object", importance=1, pretty_print_fun=pretty_print_iterable
        )


def add_breakpoints(summary, gff):
    print("=== Gene body in chromothriptic breakpoint ===")
    breakpoints = ChromothripticBreakpoints(gff)
    breakpoints.load(replace_chr=False)
    breakpoint_col_setting = [
        {"key": "Gene body in breakpoint", "offset": 0, "importance": 1},
        {"key": "Breakpoint within 10k of gene body", "offset": 1e4, "importance": 0.5},
    ]
    rows = {setting["key"]: {} for setting in breakpoint_col_setting}
    
    for bp in breakpoints:
        pairs = [(bp["query.chr"], bp["query.start"]), (bp["query.chr2"], bp["query.end"])]
        for (chr, pos), (partner_chr, partner_pos) in itertools.permutations(pairs):
            label = f"{chr}:{pos}-{partner_chr}:{partner_pos}"
            
            def update_bp_dict(key, offset):
                for gene in gff.chromosomes[chr].get_in_range(pos - offset, pos + offset, max_recursion=0):
                    geneid = gene.sanitized_id()
                    if geneid not in rows[key]:
                        rows[key][geneid] = []
                    rows[key][gene.sanitized_id()].append(label)
            
            for setting in breakpoint_col_setting:
                update_bp_dict(setting["key"], setting["offset"])
    
    for setting in breakpoint_col_setting:
        summary.add_column(
            colname=setting["key"],
            rows=rows[setting["key"]],
            dtype="object",
            importance=setting["importance"],
            pretty_print_fun=pretty_print_iterable,
        )


def add_breakpoint_fusion_support(summary, gff):
    print("=== Breakpoint supports fusion ===")
    breakpoints = ChromothripticBreakpoints(gff)
    breakpoints.load(replace_chr=False)
    fusion_series = summary["Fusion targets Primary"]
    idx = ~pd.isnull(fusion_series)
    rows = {}
    for gene_a in fusion_series.index[idx]:
        gene_a_feature = gff.get_gene(gene_a)
        if gene_a_feature is None:
            print("Can't find ", gene_a, " in gff")
            continue
        confirmed_partners = set()
        for gene_a_partner in fusion_series.loc[gene_a]:
            gene_b = gene_a_partner.gene
            gene_b_feature = gff.get_gene(gene_b)
            if gene_b_feature is None:
                print("Can't find ", gene_b, " in gff")
                continue
            shortest_path = breakpoints.get_shortest_path_between_two_genes(
                gene_a_feature, gene_b_feature, max_dist=5e5
            )
            if shortest_path.uses_breakpoint and shortest_path.is_connected:
                if len(shortest_path.path) > 4:
                    genes_involved = {
                        g.id
                        for node in shortest_path.path[2:-2]
                        for g in gff.chromosomes[node.node["chrom"]].get_in_range(
                            node.node["start"], node.node["start"] + 1, max_recursion=0
                        )
                    }
                    genes_involved = genes_involved.difference({gene_a_feature.id, gene_b_feature.id})
                    if len(genes_involved) > 0:
                        print("Fusion ", gene_a, gene_b, " also involves ", genes_involved)
                
                confirmed_partners.update(
                    {f"{gene_b} ({shortest_path.distance} bp away, {(len(shortest_path.path) - 2) // 2} jumps)"}
                )
        if len(confirmed_partners) > 0:
            rows[gene_a] = confirmed_partners
    
    summary.add_column(
        colname="Breakpoint support for Fusion",
        rows=rows,
        dtype="object",
        importance=1,
        pretty_print_fun=pretty_print_iterable,
    )


def add_double_minute(summary, gff):
    print("=== Is on double minute ===")
    rows = {}
    for coords in DoubleMinuteParts().load(replace_chr=False):
        chrom: GFFFeature = gff.chromosomes[coords["chr"]]
        for gene in chrom.get_in_range(coords["start"], coords["end"], max_recursion=0):
            rows[gene.sanitized_id()] = 1
    
    summary.add_column(
        colname="Is on Double Minute",
        rows=rows,
        dtype=pd.Int64Dtype(),
        importance=2,
    )


def add_expression_outliers(summary, gff):
    print("=== Expression Outliers ===")
    rows = {gene: direction for gene, direction in OutlierAnalysisResults().load().get_outlier_direction()}
    summary.add_column(
        colname="Expression Outlier direction",
        rows=rows,
        dtype=pd.Int64Dtype(),
        importance=0.5,
        pretty_print_fun=lambda x: {1: "+", -1: "-"}.get(x, ""),
    )


def add_tumor_copy_number(summary, gff):
    copy_numbers = load_tumor_sample_cnv(module_config.tumor_coverage_file, gff)
    rows_primary = copy_numbers["Primary"].to_dict()
    rows_relapse = copy_numbers["Relapse"].to_dict()
    summary.add_column(colname="Normalized Copy Number Primary", rows=rows_primary, dtype=float, importance=0)
    summary.add_column(colname="Normalized Copy Number Relapse", rows=rows_relapse, dtype=float, importance=0)


def add_cnv(summary, gff, vcf_blood):
    print("=== CNV ===")
    cnv = load_gene_cnv(module_config.cnv_primary_file, vcf_blood, gff, replace_chr=False)
    rows = {gene: counts[0] / sum(counts) for gene, counts in cnv.items() if sum(counts) > 0}
    summary.add_column(colname="Allelic CN ratio Primary (HP1)", rows=rows, dtype=float)
    
    cnv = load_gene_cnv(module_config.cnv_relapse_file, vcf_blood, gff, replace_chr=False)
    rows = {gene: counts[0] / sum(counts) for gene, counts in cnv.items() if sum(counts) > 0}
    summary.add_column(colname="Allelic CN ratio Relapse (HP1)", rows=rows, dtype=float)


def load_gene_ase_primary(vcf_blood, pval_thres=0.05):
    from nanoepitools.math import fdr_from_pvals
    
    ase = ASE(module_config.wasp_ase_file)
    ase.load(load_exon_annotation=True, pval_thres=None, add_chr=True)
    ase.assign_counts_to_hp(vcf_blood, "blood")
    gene_ase = ase.get_per_exon_hp_ratio().groupby("geneid")
    
    all_ase = {gene: list(gene_ase.get_group(gene)["hp1_ratio"]) for gene in gene_ase.groups.keys()}
    
    ase_maxdev = ase.get_per_exon_hp_ratio_max_effect()
    ase_maxdev["fdr"] = fdr_from_pvals(ase_maxdev["pval"])
    if pval_thres is not None:
        ase_maxdev = ase_maxdev.loc[ase_maxdev["fdr"] < pval_thres]
    
    return all_ase, ase_maxdev


def add_ase(summary, gff, vcf_blood):
    print("=== ASE ===")
    
    all_ase, ase_maxdev = load_gene_ase_primary(vcf_blood)
    summary.add_column(
        colname="ASE HP1 ratio",
        rows=all_ase,
        dtype=object,
        importance=1,
        pretty_print_fun=pretty_print_iterable,
    )
    
    summary.add_column(
        colname="ASE HP1 ratio maxdev",
        rows=ase_maxdev["hp1_ratio"].to_dict(),
        dtype=object,
        importance=0,
    )


def ase_asm_importance_function(row):
    """Provides additional importance to rows which have ASE and ASM effects, where the
    ASE is not sufficiently explained by the CNV

    Two conditions are checked:
     A) does ASE and ASM point in the "expected" direction? (less methylation -> more expression)
     B) Does CNV sufficiently explain ASE?
    returns:
      * 0 if there is no ASM
      * 0.5 if B but not A
      * 1 if not B and not A
      * 1 if A and B
      * 2 if A but not B

    """
    ase = row["ASE HP1 ratio maxdev"]
    cnv = row["Allelic CN ratio Primary (HP1)"]
    if np.isnan(ase):
        return 0
    
    importance_candidates = []
    for asm in row["Promoters Diffmet Haplotype Primary (HP 2-1)"], row["Enhancers Diffmet Haplotype Primary (HP 2-1)"]:
        if pd.isnull(asm) or np.isnan(asm):
            importance_candidates.append(0)
        
        if np.sign(asm) == np.sign(ase - 0.5):
            # dmr positive -> HP1 is demethylated -> HP1 expression is expected to be higher
            multiplier = 3
        else:
            # dmr is different direction from ase
            multiplier = 1
        
        if abs(ase - cnv) < 0.25:
            # cnv explains ase sufficiently
            importance_candidates.append(1 * multiplier)
        else:
            # cnv does not explain ase
            importance_candidates.append(2 * multiplier)
    
    return max(importance_candidates)


def add_cosmic_annotation(summary, gff):
    gene_name_translator = GeneNameToEnsemblID(gff)
    subtype_re = re.compile("^cosmic_mutations_([^_]*)_subtype.tsv$")
    rows = {}
    for cosmic_tsv in Path(module_config.literature_annotation_dir).iterdir():
        mitch = subtype_re.match(cosmic_tsv.name)
        if mitch is None:
            continue
        gene_names_in_annotation = list(set([l.split("\t")[0].split("_")[0].strip() for l in open(cosmic_tsv, "r")]))
        print(gene_names_in_annotation)
        gene_ids_in_annotation = gene_name_translator.translate_gene_names_to_ids(gene_names_in_annotation)
        subtype = mitch.group(1)
        for gene in gene_ids_in_annotation:
            if gene in rows:
                rows[gene] = f"{rows[gene]},{subtype}"
            else:
                rows[gene] = subtype
    
    def importance_fun(val):
        if pd.isnull(val):
            return 0
        elif "ssh" in val.split(","):
            return 2
        else:
            return len(val.split(",")) * 0.25
    
    summary.add_column("cosmic_has_mutation", rows, pd.Int64Dtype, importance_fun=importance_fun)


def add_literature_annotation(summary, gff):
    annotation_files = {
        "Main Drivers": {"filename": "pap1_main_driver_genes.txt", "importance": 0.5},
        "SHH Amplified": {"filename": "pap1_pap2_shh_copy_number_amp_genes.txt", "importance": 0.1},
        "SHH Deleted": {"filename": "pap1_pap2_shh_copy_number_del_genes.txt", "importance": 0.1},
        "Drug target": {"filename": "pap2_drug_target_genes.txt", "importance": 0.5},
        "MB Amplified": {"filename": "pap2_mb_copy_number_amp_genes.txt", "importance": 0.1},
        "MB Deleted": {"filename": "pap2_mb_copy_number_del_genes.txt", "importance": 0.1},
        # "SNVs COSMIC": {"filename": "pap3_coding_snv_cosmic_genes.txt", "importance": 0.1},
        # "Significant SNV": {"filename": "pap3_most_significant_mutation_genes.txt", "importance": 0.25},
    }
    gene_name_translator = GeneNameToEnsemblID(gff)
    for colname, col_settings in annotation_files.items():
        
        gene_list_file = Path(module_config.literature_annotation_dir).joinpath(col_settings["filename"])
        gene_names_in_annotation = list(set([l.strip().replace('"', "") for l in open(gene_list_file, "r")]))
        gene_ids_in_annotation = gene_name_translator.translate_gene_names_to_ids(gene_names_in_annotation)
        for name, id in zip(gene_names_in_annotation, gene_ids_in_annotation):
            if id is None:
                print("Unknown gene: ", name)
        
        rows = {id: 1 for id in gene_ids_in_annotation if id is not None}
        summary.add_column(colname, rows, pd.Int64Dtype, importance=col_settings["importance"])


def load_vcf_blood():
    vcf_blood = PhasedVCF(module_config.vcf_blood_file)
    vcf_blood.read()
    return vcf_blood


def build_summary(gff, columns="all"):
    summary = SummaryTable(gff)
    
    vcf_blood = None
    
    
    
    if columns == "all" or "diffmet_promoter" in columns:
        add_diffmet_promoter_primary_relapse(summary, gff)
    
    if columns == "all" or "diffmet_enhancer" in columns:
        add_diffmet_enhancer_primary_relapse(summary, gff)
    
    if columns == "all" or "diff_exp" in columns:
        add_diffexp(summary, gff)
    
    if columns == "all" or "diffmet_haplotype" in columns:
        add_diffmet_haplotype(summary, gff)
    
    if columns == "all" or "ase" in columns:
        if vcf_blood is None:
            vcf_blood = load_vcf_blood()
        add_ase(summary, gff, vcf_blood)
    
    if columns == "all" or "cnv" in columns:
        if vcf_blood is None:
            vcf_blood = load_vcf_blood()
        add_cnv(summary, gff, vcf_blood)
    
    if columns == "all" or "copy_number" in columns:
        add_tumor_copy_number(summary, gff)
        
    if columns == "all" or "expression_outlier" in columns:
        add_expression_outliers(summary, gff)
    
    if columns == "all" or "gene_fusion" in columns:
        add_gene_fusion(summary, gff)
    
    if columns == "all" or "breakpoints" in columns:
        add_breakpoints(summary, gff)
    
    if columns == "all" or "breakpoint_fusion_support" in columns:
        add_breakpoint_fusion_support(summary, gff)
    
    if columns == "all" or "double_minute" in columns:
        add_double_minute(summary, gff)
    
    if columns == "all" or "literature" in columns:
        add_literature_annotation(summary, gff)
    
    if columns == "all" or "cosmic" in columns:
        add_cosmic_annotation(summary, gff)
    
    if columns == "all" or all([x in columns for x in ("cnv", "ase", "diffmet_haplotype")]):
        summary.add_combined_importance_function("asm_ase_cnv", ase_asm_importance_function)
    
    return summary


if __name__ == "__main__":
    print(
        """
    ==========================================================
    Loading annotations...
    ==========================================================
    """
    )
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    summary_met = build_summary(
        gff,
        columns=[
            "diffmet_promoter",
            "diffmet_enhancer",
            "diffmet_haplotype",
            "diff_exp",
            "ase",
            "literature",
            "cosmic",
        ],
    )
    summary_met.to_csv(module_config.grand_summary_file_methylation, sep="\t", min_importance=-100)
    
    summary = build_summary(gff)
    summary.to_csv(module_config.grand_summary_file, sep="\t", min_importance=-100)
    
    print("Quickly compute enrichment of main drivers in dmr:")
    n_dmr_driver = (
        ~summary_met.summary["5k from TSS diffmet (Relapse-Primary)"].isna()
        & ~summary_met.summary["Main Drivers"].isnull()
    ).sum()
    n_drivers = (~summary_met.summary["Main Drivers"].isnull()).sum()
    n_dmr = (~summary_met.summary["5k from TSS diffmet (Relapse-Primary)"].isnull()).sum()
    num_genes = sum(1 for c in gff.chromosomes.values() for g in c.children)
    
    import scipy
    
    n_dmr_notdrivr = n_dmr - n_dmr_driver
    n_notdmr_driver = n_drivers - n_dmr_driver
    n_notdmr_notdriver = num_genes - n_notdmr_driver - n_dmr_notdrivr + n_dmr_driver
    print(
        scipy.stats.fisher_exact(
            [[n_dmr_driver, n_dmr_notdrivr], [n_notdmr_driver, n_notdmr_notdriver]], alternative="two-sided"
        )
    )
    
    indicated_in_medulloblastoma = (
        ~summary.summary["Main Drivers"].isna()
        | ~summary.summary["SHH Amplified"].isna()
        | ~summary.summary["SHH Deleted"].isna()
        | ~summary.summary["Drug target"].isna()
        | ~summary.summary["MB Amplified"].isna()
        | ~summary.summary["MB Deleted"].isna()
        | ~summary.summary["cosmic_has_mutation"].isna()
    )
    print(
        "Primary ASE indicated in MB: ",
        (~summary.summary.loc[indicated_in_medulloblastoma]["ASE HP1 ratio maxdev"].isna()).sum(),
    )
    
    idx = ~summary.summary["ASE HP1 ratio maxdev"].isna()
    
    def cnv_explains_ase(row):
        is_pos = (row["Allelic CN ratio Primary (HP1)"] > 0.65) & (row["ASE HP1 ratio maxdev"] > 0.5)
        is_neg = (row["Allelic CN ratio Primary (HP1)"] < 0.35) & (row["ASE HP1 ratio maxdev"] < 0.5)
        return is_pos or is_neg
    
    ase_cnv = summary.summary.loc[idx][["ASE HP1 ratio maxdev", "Allelic CN ratio Primary (HP1)"]].apply(cnv_explains_ase, axis=1)
    
    num_ase = idx.sum()
    num_cnv = summary.summary["Allelic CN ratio Primary (HP1)"].map(lambda x: abs(x - 0.5) > 0.15).sum()
    num_ase_cnv = ase_cnv.sum()
    num_ase_not_cnv = num_ase - num_ase_cnv
    num_cnv_not_ase = num_cnv - num_ase_cnv
    num_genes = summary.summary.shape[0]
    num_not_ase_not_cnv = num_genes - num_ase_not_cnv - num_cnv_not_ase - num_ase_cnv
    
    scipy.stats.fisher_exact(
        [[num_ase_cnv, num_ase_not_cnv], [num_cnv_not_ase, num_not_ase_not_cnv]], alternative="two-sided"
    )
