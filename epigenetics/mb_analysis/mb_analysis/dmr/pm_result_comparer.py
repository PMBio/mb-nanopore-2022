from typing import Dict, Optional
from matplotlib_venn import venn3, venn2
import matplotlib
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.config import module_config
from nanoepitools.pycometh_result import PycomethOutput
from nanoepitools.annotations.enhancers import Enhancers
from mb_analysis.summary.plot_gene import Plotter

matplotlib.use("Agg")


def merge_duplicate_diffmet_hits(oldhits):
    # Merge so we don't count them double
    newhits = []
    a = 0
    b = 0
    for i, hiti in enumerate(oldhits):
        a += 1
        duplicate = -1
        for j, hitj in enumerate(newhits):
            if hiti["start"] < hitj["end"] and hitj["start"] < hiti["end"] and hiti["chrom"] == hitj["chrom"]:
                duplicate = j
                break
        if duplicate >= 0:
            newhits[duplicate] = {
                "chrom": hiti["chrom"],
                "start": min(hiti["start"], newhits[duplicate]["start"]),
                "end": max(hiti["end"], newhits[duplicate]["end"]),
            }
        else:
            b += 1
            newhits.append(hiti)
    return newhits


def find_overlap(lista, listb):
    overlap = []
    diff = []
    for a in lista:
        found = False
        for b in listb:
            if a["start"] - 1000 < b["end"] and b["start"] - 1000 < a["end"]:
                overlap += [a, b]
                found = True
                break
        if not found:
            diff.append(a)
    return overlap, diff


class PMComparer:
    def __init__(self, samples: Dict[str, str], project: str, gene_maxdist=5000, sample_settings=None):
        print("Reading GFF")
        self.gff = GFFAnnotationsReader()
        self.gff.read(module_config.gff_file, only_protein_coding=False)
        self.gff.build_index()
        print("Reading Enhancers")
        self.enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
        self.enhancers.load()
        self.enhancers.annotate_nearest_gene(self.gff, maxdist=gene_maxdist)
        self.enhancers.filter_nearest_gene_none()
        print("Preparing plotter")
        self.pa = PlotArchiver(project, config={"plot_archive_dir": "/homes/snajder/data1/plots_medulloblastoma"})
        self.pl = Plotter(self.gff, self.pa, enhancers=self.enhancers)
        self.samples = samples
        self.sample_settings = sample_settings
        self.pm: Optional[Dict[str, PycomethOutput]] = None
        self.promoters_hit = None
        self.hits = None
    
    def load(self):
        self.pm = {}
        for sample in self.samples:
            self.pm[sample] = PycomethOutput(self.samples[sample])
    
    def load_hits(self, **kwargs):
        if self.pm is None:
            self.load()
        
        self.hits = {}
        for sample in self.pm:
            hits = [
                {
                    "chrom": line["chromosome"],
                    "start": line["start"],
                    "end": line["end"],
                    "diff": line["diff"],
                    "fdr": line["adj_pvalue"],
                }
                for line in self.pm[sample].read_file(**kwargs)
            ]
            self.hits[sample] = merge_duplicate_diffmet_hits(hits)
    
    def load_promoters_hit(self, *args, **kwargs):
        if self.sample_settings is None:
            """backwards compatibility"""
            return self.load_promoters_hit_override_settings(*args, **kwargs)
        
        if self.pm is None:
            self.load()
        
        self.promoters_hit = {}
        for sample, settings in self.sample_settings.items():
            self.promoters_hit[sample] = self.pm[sample].load_promoters_hit(
                self.gff,
                settings["promoter_before_tss"],
                settings["promoter_after_tss"],
                b_minus_a=settings["b_minus_a"],
                drop_insignificant=settings["drop_insignificant"],
                pval_threshold=settings["pval_threshold"],
                progress=settings["progress"],
                min_diff=settings["min_diff"],
            )
    
    def load_promoters_hit_override_settings(
        self,
        before=2000,
        after=500,
        b_minus_a=True,
        drop_insignificant=True,
        pval_threshold=0.05,
        progress=True,
        min_diff=0.25,
    ):
        if self.pm is None:
            self.load()
        
        self.promoters_hit = {}
        for sample in self.samples:
            
            self.promoters_hit[sample] = self.pm[sample].load_promoters_hit(
                self.gff,
                before,
                after,
                b_minus_a=b_minus_a,
                drop_insignificant=drop_insignificant,
                pval_threshold=pval_threshold,
                progress=progress,
                min_diff=min_diff,
            )
    
    def load_enhancers_hit(
        self, b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, progress=True, min_diff=0.25
    ):
        if self.pm is None:
            self.load()
        
        self.enhancers_hit = {}
        for sample in self.samples:
            self.enhancers_hit[sample] = self.pm[sample].load_enhancers_hit(
                self.enhancers,
                b_minus_a=b_minus_a,
                drop_insignificant=drop_insignificant,
                pval_threshold=pval_threshold,
                progress=progress,
                min_diff=min_diff,
            )
    
    def plot_promoter_venn(self, fig_key=None):
        if self.promoters_hit is None:
            self.load_promoters_hit()
        
        self.pa.figure()
        samples = list(self.samples.keys())
        if len(samples) == 3:
            venn3([set(self.promoters_hit[sample].keys()) for sample in samples], samples)
        elif len(samples) == 2:
            venn2([set(self.promoters_hit[sample].keys()) for sample in samples], samples)
        if fig_key is None:
            self.pa.savefig()
        else:
            self.pa.savefig(fig_key)
    
    def plot_enhancer_venn(self, fig_key=None):
        if self.promoters_hit is None:
            self.load_promoters_hit()
        
        self.pa.figure()
        samples = list(self.samples.keys())
        if len(samples) == 3:
            venn3([set(self.enhancers_hit[sample].keys()) for sample in samples], samples)
        elif len(samples) == 2:
            venn2([set(self.enhancers_hit[sample].keys()) for sample in samples], samples)
        if fig_key is None:
            self.pa.savefig()
        else:
            self.pa.savefig(fig_key)
