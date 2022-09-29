import matplotlib.pyplot as plt

import tqdm
from nanoepitools.plotting.general_plotting import PlotArchiver
from mb_analysis.summary.plot_gene import Plotter
from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.annotations.enhancers import Enhancers
from mb_analysis.config import module_config
from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits

"""
Class for plotting methylation profiles of interesting genes and regions
"""


class PromoterDiffmetPlotter:
    """
    This class loads promoter mapped DMRs from pycometh outputs into memory so we can quickly plot
    the sections that contain them
    """
    
    def __init__(self):
        files = {
            "sc_hmm": module_config.pycometh_primary_relapse_file_hmm,
            "sc_cgi": module_config.pycometh_primary_relapse_file_cgi,
            "asm_primary": module_config.pycometh_haplotype_sample_template_file.format(sample="Primary"),
        }
        self.promoters_hit = {
            key: PycomethOutput(met_comp_file=val).load_promoters_hit(
                gff,
                2000,
                500,
                min_diff=0.5,
                b_minus_a=True,
                drop_insignificant=False,
                pval_threshold=0.05,
                progress=True,
            )
            for key, val in files.items()
        }
    
    def get_gene_hits(self, gene_id, sample_comp):
        hit_set = (
            [self.promoters_hit["sc_hmm"], self.promoters_hit["sc_cgi"]]
            if sample_comp
            else [self.promoters_hit["asm_primary"]]
        )
        for hits in hit_set:
            if gene_id in hits:
                for hit in hits[gene_id]:
                    yield hit
    
    def plot(
        self,
        genes,
        sample_comp=True,
        additional_annotations=[],
        figure_kwargs=None,
        bpbefore=3000,
        bpafter=3000,
        **kwargs
    ):
        if figure_kwargs is None:
            figure_kwargs = {"figsize": (8, 4)}
        
        for geneid in genes:
            hits = list(self.get_gene_hits(geneid, sample_comp))
            if len(hits) == 0:
                continue
            start = min(hit["start"] for hit in hits) - bpbefore
            end = max(hit["end"] for hit in hits) + bpafter
            chrom = hits[0]["chrom"]
            annotations = [
                {"region": [hit["start"], hit["end"]], "text": "Differential methylation", "color": "r"} for hit in hits
            ] + additional_annotations
            print(annotations)
            pl.plot_region(chrom, start, end, figure_kwargs=figure_kwargs, ws=0, annotations=annotations, **kwargs)


class RegionDiffmetPlotter:
    """
    This class loads all DMRs without mapping them to a promoter (this loads much faster) and can be used for plotting
    regions of interest that aren't necessarily promoter linked
    """
    
    def __init__(self):
        files = {
            "sc_hmm": module_config.pycometh_primary_relapse_file_hmm,
            "sc_cgi": module_config.pycometh_primary_relapse_file_cgi,
            "asm_primary": module_config.pycometh_haplotype_sample_template_file.format(sample="Primary"),
        }
        
        self.hits = {
            key: list(
                PycomethOutput(met_comp_file=val).read_file(
                    min_diff=0.35,
                    b_minus_a=True,
                    drop_insignificant=False,
                    pval_threshold=0.05,
                    progress=True,
                )
            )
            for key, val in files.items()
        }
        for k, v in self.hits.items():
            for hit in v:
                hit["chrom"] = hit["chromosome"]
        
        self.hits = {
            "sc": merge_duplicate_diffmet_hits(self.hits["sc_hmm"] + self.hits["sc_cgi"]),
            "asm_primary": merge_duplicate_diffmet_hits(self.hits["asm_primary"]),
        }
    
    def plot(self, chrom, start, end, sample_comp=True, additional_annotations=[], figure_kwargs=None, **kwargs):
        hits = self.hits["sc"] if sample_comp else self.hits["asm_primary"]
        annotations = additional_annotations
        for hit in hits:
            if hit["chrom"] == chrom and start < hit["end"] and hit["start"] < end:
                annotations.append(
                    {"region": [hit["start"], hit["end"]], "text": "Differential methylation", "color": "r"}
                )
        print(annotations)
        pl.plot_region(chrom, start, end, figure_kwargs=figure_kwargs, ws=0, annotations=annotations, **kwargs)


if __name__ == "__main__":
    pa = PlotArchiver("summary", config=module_config)
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
    enhancers.load(replace_chr=False)
    enhancers.annotate_nearest_gene(gff, maxdist=5000)
    enhancers.filter_nearest_gene_none()
    
    # module_config.samples = ["Primary"]
    pl = Plotter(gff, pa, enhancers=enhancers)
    
    rplot = RegionDiffmetPlotter()
    pdplot = PromoterDiffmetPlotter()
    # default figure size
    figure_kwargs = {"figsize": (8, 4)}
    
    with pa.open_multipage_pdf("genes_promoter_diffmet_chromatin_organization"):
        genes = [
            "ENSG00000183207",
            "ENSG00000275379",
            "ENSG00000275126",
            "ENSG00000276410",
            "ENSG00000197153",
            "ENSG00000105447",
        ]
        pdplot.plot(genes)
    
    with pa.open_multipage_pdf("genes_promoter_diffmet_tp53_binding_related"):
        genes = ["ENSG00000174775", "ENSG00000141736", "ENSG00000153707", "ENSG00000196367"]
        pdplot.plot(genes)
    
    with pa.open_multipage_pdf("genes_diffmet_on_dm_downmet_region"):
        pl.plot_gene(
            "ENSG00000166436",
            figure_kwargs=figure_kwargs,
            annotations=[],
            ws=0,
            title="TRIM66",
            min_marker_width_relative=0.02,
            with_germline=False,
        )
        
        pl.plot_gene(
            "ENSG00000130413",
            figure_kwargs=figure_kwargs,
            annotations=[],
            ws=0,
            title="STK33",
            min_marker_width_relative=0.02,
            with_germline=False,
        )
    
    with pa.open_multipage_pdf("tp53"):
        pl.plot_region("chr17", 7668944-50000, 7692092+50000, figure_kwargs=figure_kwargs, annotations=[], ws=0, title="TP53",
                       aggregate_samples=False, with_germline=True, with_no_hp=True, coordinate_space=True,
                       marker_height=0.9, fontsize=8, show_has_value_indicator=False)
        
        
    with pa.open_multipage_pdf("diffmet_example_genes"):
        rplot.plot("chr14",41726997-2000, 41728947+2000,  figure_kwargs=figure_kwargs, title="Chr11",
                       aggregate_samples=False, with_germline=True, with_no_hp=True, coordinate_space=False,
                       marker_height=0.9, fontsize=8, show_has_value_indicator=False, additional_annotations=[], sample_comp=True)
        
        pl.plot_region("chr10", 75398609-2000, 75398999+2000, figure_kwargs=figure_kwargs, annotations=[], ws=0,
            title="MAP3K14", aggregate_samples=False, with_germline=False, with_no_hp=True, coordinate_space=False,
            marker_height=0.9, fontsize=8, show_has_value_indicator=False)
    

        genes = ["ENSG00000166436",  # TRIM66
                 "ENSG00000006062" # SPATA32
        ]
        pdplot.plot(genes, sample_comp=False, aggregate_samples=True, with_germline=False, with_no_hp=False,
                    coordinate_space=False, marker_height=0.9, fontsize=8, show_has_value_indicator=True, )
        genes = [
            "ENSG00000124785", # NRN1
            "ENSG00000185920", # PTCH1
        ]
        pdplot.plot(genes, sample_comp=True, aggregate_samples=True, with_germline=False, with_no_hp=False,
                    coordinate_space=False, marker_height=0.9, fontsize=8, show_has_value_indicator=False,
        )
        genes = [
            "ENSG00000118946" # PCDH17
        ]
        pdplot.plot(genes, sample_comp=False,
            aggregate_samples=True,
            with_germline=False,
            with_no_hp=False,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=False,
        )

    
    with pa.open_multipage_pdf("genes_of_interest"):
        figure_kwargs = {"figsize": (8, 4)}
        annotations = [{'region': [6002148, 6003069], 'text': 'Differential methylation', 'color': 'r'}]
        pl.plot_region(
            "chr6",
            5999000,
            6009050,
            figure_kwargs=figure_kwargs, annotations=annotations,
            ws=0,
            title="NRN1",
            aggregate_samples=True,
            with_germline=False,
            with_no_hp=False,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=True,
        )
        
        # PCDH17
        annotations = [{"region": [57633174, 57635093], "text": "Differential methylation", "color": "r"}]
        pl.plot_region(
            "chr13",
            57630000,
            57634998,
            figure_kwargs=figure_kwargs,
            annotations=annotations,
            ws=0,
            title="PCDH17",
            aggregate_samples=True,
            coordinate_space=False,
            with_germline=False,
            with_no_hp=False,
            marker_height=0.9,
            fontsize=8, show_has_value_indicator=True,
        )
        
        # PTCH1
        annotations = [{"region": [95503519, 95504075], "text": "Differential methylation", "color": "r"}]
        pl.plot_region(
            "chr9",
            95498525,
            95508336,
            figure_kwargs=figure_kwargs,
            annotations=annotations,
            ws=0,
            title="PTCH1",
            aggregate_samples=True,
            coordinate_space=False,
            with_germline=False,
            with_no_hp=False,
            marker_height=0.9,
            fontsize=8, show_has_value_indicator=True,
        )

    
    with pa.open_multipage_pdf("mb_known_genes"):
        pdplot.plot(["ENSG00000185920"], with_germline=True)
        pdplot.plot(["ENSG00000205212"], with_germline=True, sample_comp=False)
    
    with pa.open_multipage_pdf("BASP1"):
        pl.plot_region(
            "chr17",
            17229000 - 5000,
            17229000 + 5000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[],
        )

    with pa.open_multipage_pdf("chr11_17_ctg2_met"):
        pl.plot_region("chr11", 8430000-25000, 8640000+25000, ws=0, fontsize=8, aggregate_samples=True,
            with_germline=True, with_relapse=True, with_no_hp=False, hold=True, figure_kwargs=dict(figsize=(16,6)))
        plt.vlines([8437500, 8636500], 0, plt.ylim()[1])
        pa.savefig()
        
    with pa.open_multipage_pdf("chr11_17_ctg2_adjacent_met"):
        pl.plot_region("chr17", 22923858, 22950292, ws=0, fontsize=8, aggregate_samples=False,
            with_germline=True, with_relapse=True, with_no_hp=True, figure_kwargs=dict(figsize=(16,6)))
        pl.plot_region("chr17", 26608681, 26631450, ws=0, fontsize=8, aggregate_samples=False, with_germline=True,
                       with_relapse=True, with_no_hp=True, figure_kwargs=dict(figsize=(16, 6)))


    with pa.open_multipage_pdf("MYPOP"):
        pl.plot_region(
            "chr19",
            45900356 - 5000,
            45900356 + 5000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [45900356 - 10, 45900356 + 10], "text": "Breakpoint", "color": "red"}],
        )
        pl.plot_region(
            "chr19",
            58586563 - 1000,
            58586648 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [58586563, 58586648], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "chr19",
            8057491 - 1000,
            8057746 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [8057491, 8057746], "text": "Insertion", "color": "red"}],
        )
    
    with pa.open_multipage_pdf("TLL1_and_insertions"):
        pl.plot_region(
            "4",
            166012532 - 5000,
            166012532 + 5000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [166012532 - 10, 166012532 + 10], "text": "Breakpoint", "color": "red"}],
        )
        pl.plot_region(
            "4",
            123601849 - 1000,
            123602312 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [123601849, 123602312], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "4",
            123225185 - 1000,
            123225497 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [123225185, 123225497], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "4",
            123116140 - 1000,
            123127028 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=True,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [123116140, 123127028], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "5",
            96174945 - 1000,
            96175172 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [96174945, 96175172], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "7",
            18822957 - 1000,
            18823294 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [18822957, 18823294], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "7",
            24780465 - 1000,
            24780660 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [24780465, 24780660], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "7",
            29857180 - 1000,
            29857578 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [29857180, 29857578], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "7",
            30812714 - 1000,
            30813149 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [30812714, 30813149], "text": "Insertion", "color": "red"}],
        )
        pl.plot_region(
            "7",
            7803367 - 1000,
            7803801 + 1000,
            ws=0,
            fontsize=8,
            aggregate_samples=False,
            with_germline=True,
            with_relapse=True,
            with_no_hp=True,
            annotations=[{"region": [7803367, 7803801], "text": "Insertion", "color": "red"}],
        )
    
    with pa.open_multipage_pdf("PPAR-G"):
        additional_annotations = [
            {"region": [12286828, 12295611], "text": "Promoter", "color": "brown"},
            {"region": [12351000, 12351427], "text": "Promoter", "color": "brown"},
            {"region": [12191943, 12195814], "text": "Enhancer", "color": "green"},
            {"region": [12307262, 12309889], "text": "Enhancer", "color": "green"},
            {"region": [12344801, 12346199], "text": "Enhancer", "color": "green"},
        ]
        pdplot.plot(
            ["ENSG00000132170"],
            figure_kwargs=figure_kwargs,
            aggregate_samples=True,
            coordinate_space=False,
            with_germline=True,
            with_no_hp=True,
            marker_height=0.9,
            fontsize=8,
            sample_comp=True,
            show_has_value_indicator=False,
            additional_annotations=additional_annotations,
        )
    
    with pa.open_multipage_pdf("TBX1"):
        figure_kwargs = {"figsize": (10, 4)}
        additional_annotations = [
            {"region": [19757000, 19757201], "text": "Ensembl promoter", "color": "g"},
            {"region": [19757800, 19758201], "text": "Ensembl promoter", "color": "g"},
            {"region": [19759600, 19764200], "text": "Ensembl promoter", "color": "g"},
        ]
        pl.draw_transcripts_separately = True
        rplot.plot(
            "chr22",
            19755570,
            19765181,
            figure_kwargs=figure_kwargs,
            title="TBX1",
            aggregate_samples=True,
            with_germline=False,
            with_no_hp=False,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=True,
            sample_comp=True,
            additional_annotations=additional_annotations,
        )
        pl.draw_transcripts_separately = False
        
        pl.draw_transcripts_separately = True
        rpl.plot(
            ["ENSG00000184058"],
            with_germline=False,
            aggregate_samples=True,
            with_no_hp=False,
            coordinate_space=False,
            figure_kwargs={"figsize": (10, 5)},
            additional_annotations=additional_annotations,
        )
        pl.draw_transcripts_separately = False
    
    with pa.open_multipage_pdf("PCDH17_full"):
        rplot.plot(
            "13",
            57630104 - 2000,
            57729311 + 2000,
            figure_kwargs={"figsize": (20, 20)},
            title="PCDH17",
            aggregate_samples=False,
            with_germline=False,
            with_no_hp=False,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=False,
            sample_comp=False,
            additional_annotations=[{"region": (57706143 - 100, 57706143 + 100), "text": "SNP", "color": "r"}],
        )
    
    with pa.open_multipage_pdf("asm_ase_genes"):
        idx = plot_data["shapes"] == "*"
        for gene in plot_data["gene_ids"][idx]:
            gene = gff.get_gene(gene)
            rplot.plot(
                gene.parent.name,
                gene.start - 5000,
                gene.end + 5000,
                figure_kwargs={"figsize": (20, 20)},
                title=gene.name,
                aggregate_samples=False,
                with_germline=False,
                with_no_hp=True,
                coordinate_space=False,
                marker_height=0.9,
                fontsize=8,
                show_has_value_indicator=False,
                sample_comp=False,
                additional_annotations=[],
            )
    
    with pa.open_multipage_pdf("HLA-DQA1"):
        rplot.plot(
            "6",
            32624450 - 2000,
            32652057 + 2000,
            figure_kwargs={"figsize": (20, 20)},
            title="HLA-DQA1",
            aggregate_samples=False,
            with_germline=True,
            with_no_hp=True,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=False,
            sample_comp=False,
        )
    
    with pa.open_multipage_pdf("AC138761.1_deletion"):
        pl.plot_region(
            "17",
            22278000,
            22300000,
            figure_kwargs=dict(figsize=(20, 20)),
            ws=0,
            fontsize=8,
            title="deletion around AC138761.1",
            aggregate_samples=False,
            with_germline=True,
            with_no_hp=True,
            coordinate_space=False,
        )
