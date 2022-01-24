from mb_analysis.config import module_config
from mb_analysis.sample_met_matrix import MetMatrixLoader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from nanoepitools.plotting.plot_methylation_profile import plot_met_profile
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.annotations.enhancers import Enhancers
from meth5.meth5 import MetH5File, ChromosomeContainer


class Plotter:
    """
    Helper class for plotting methylation matrices including genomic annotation and separated by sample
    """
    
    def __init__(
        self, gff: GFFAnnotationsReader, pa: PlotArchiver, enhancers: Enhancers, draw_transcripts_separately=False
    ):
        self.loader = MetMatrixLoader()
        self.pa = pa
        self.gff = gff
        self.enhancers = enhancers
        self.draw_transcripts_separately = draw_transcripts_separately
        self.sample_colors = {
            "Primary (HP1)": "#DF566E",
            "Primary (HP2)": "#DF566E",
            "Primary": "#DF566E",
            "Germline (HP1)": "g",
            "Germline (HP2)": "g",
            "Germline": "g",
            "Relapse (HP1)": "#9747CF",
            "Relapse (HP2)": "#9747CF",
            "Relapse": "#9747CF",
        }
        self.sample_hatch = {
            "Primary (HP1)": "//",
            "Primary (HP2)": "\\\\",
            "Primary": "",
            "Germline (HP1)": "//",
            "Germline (HP2)": "\\\\",
            "Germline": "",
            "Relapse (HP1)": "//",
            "Relapse (HP2)": "\\\\",
            "Relapse": "",
        }
    
    def draw_region_annotation_track(
        self, text, regions, y_offset, fontsize=8, linewidth=2, barcolor="b", coordinate_translator=None
    ):
        if len(regions) == 0:
            return
        min_region_x = min([r[0] for r in regions])
        max_region_x = max([r[1] for r in regions])
        if coordinate_translator is not None:
            min_region_x = coordinate_translator(min_region_x)
            max_region_x = coordinate_translator(max_region_x)
        
        xlim = plt.xlim()
        plt.xlim(xlim)
        
        axis_range = xlim[1] - xlim[0]
        min_text_x = xlim[0] + axis_range * 0.1
        max_text_x = xlim[1] - axis_range * 0.3
        text_x = min_region_x - axis_range * 0.0075
        
        if text_x < min_text_x:
            ha = "left"
            text_x = max_region_x + axis_range * 0.0075
            if text_x > max_text_x:
                text_x = min_text_x
        else:
            ha = "right"
        for region in regions:
            if coordinate_translator is not None:
                region = [coordinate_translator(r) for r in region]
            rect = patches.Rectangle(
                (region[0], y_offset - linewidth / 2),
                region[1] - region[0],
                linewidth,
                linewidth=0,
                facecolor=barcolor,
            )
            plt.gca().add_patch(rect)
        
        plt.text(
            text_x,
            y_offset,
            text,
            ha=ha,
            va="center",
            fontsize=fontsize,
            bbox=dict(boxstyle="square,pad=0", fc="w", ec="w", alpha=0.85),
        )
    
    def draw_region_annotation(
        self, text, region, y_offset, fontsize=8, linewidth=2, barcolor="b", coordinate_translator=None
    ):
        xlim = plt.xlim()
        plt.xlim(xlim)
        
        if coordinate_translator is not None:
            region = [coordinate_translator(r) for r in region]
        
        region = [max(xlim[0], region[0]), min(xlim[1], region[1])]
        axis_range = xlim[1] - xlim[0]
        min_text_x = xlim[0] + axis_range * 0.1
        max_text_x = xlim[1] - axis_range * 0.3
        text_x = region[0] - axis_range * 0.0075
        
        if text_x < min_text_x:
            ha = "left"
            text_x = region[1] + axis_range * 0.0075
            if text_x > max_text_x:
                text_x = min_text_x
        else:
            ha = "right"
        
        plt.text(
            text_x,
            y_offset,
            text,
            ha=ha,
            va="center",
            fontsize=fontsize,
            bbox=dict(boxstyle="square,pad=0", fc="w", ec="w", alpha=0.85),
        )
        
        rect = patches.Rectangle(
            (region[0], y_offset - linewidth / 2), region[1] - region[0], linewidth, linewidth=0, facecolor=barcolor
        )
        plt.gca().add_patch(rect)
    
    def plot_annotations(
        self,
        chrom,
        start,
        end,
        additional_annotations=[],
        coordinate_translator=None,
        promoter_before_tss=2000,
        promoter_after_tss=500,
        fontsize=10,
    ):
        ylim = plt.ylim()
        annotation_offset_step = (ylim[1] - ylim[0]) * 0.08
        annotation_offset = -annotation_offset_step
        
        xlim = plt.xlim()
        plt.gca().margins(x=0.05)
        
        annotation_offset -= annotation_offset_step
        linewidth = annotation_offset_step * 0.5
        
        for annotation_properties in additional_annotations:
            annotation_offset -= annotation_offset_step
            if annotation_properties["region"][0] > end:
                continue
            if annotation_properties["region"][1] < start:
                continue
            self.draw_region_annotation(
                annotation_properties["text"],
                annotation_properties["region"],
                annotation_offset,
                linewidth=linewidth,
                fontsize=fontsize,
                barcolor=annotation_properties["color"],
                coordinate_translator=coordinate_translator,
            )
        
        # Draw genes and exons
        chr_genes = self.gff.chromosomes[chrom]
        
        genes_in_view = [f for f in chr_genes.get_in_range(start, end, max_recursion=0)]
        for gene in genes_in_view:
            region = [gene.start, gene.end]
            annotation_offset -= annotation_offset_step
            direction_symbol = "<<" if gene.direction == "-" else ">>"
            text_label = f"{direction_symbol} {gene.name} {direction_symbol}"
            if not self.draw_transcripts_separately:
                self.draw_region_annotation(
                    text_label,
                    region,
                    annotation_offset,
                    barcolor="b",
                    linewidth=linewidth * 0.25,
                    fontsize=fontsize,
                    coordinate_translator=coordinate_translator,
                )
            for transcript in gene.sorted_children:
                for exon in transcript.sorted_children:
                    if exon.end < start or exon.start > end:
                        continue
                    
                    region = [exon.start, exon.end]
                    
                    self.draw_region_annotation(
                        "",
                        region,
                        annotation_offset,
                        barcolor="b",
                        linewidth=linewidth,
                        fontsize=fontsize,
                        coordinate_translator=coordinate_translator,
                    )
                
                if self.draw_transcripts_separately:
                    self.draw_region_annotation(
                        text_label,
                        (transcript.start, transcript.end),
                        annotation_offset,
                        barcolor="b",
                        linewidth=linewidth * 0.25,
                        fontsize=fontsize,
                        coordinate_translator=coordinate_translator,
                    )
                    annotation_offset -= annotation_offset_step
        
        # Draw promoters
        transcripts_in_view = [
            f for f in chr_genes.get_in_range(start - promoter_before_tss, end + promoter_before_tss, max_recursion=1)
        ]
        annotation_offset -= annotation_offset_step
        annotation_regions = []
        for transcript in transcripts_in_view:
            if transcript.direction == "+":
                region = [transcript.start - promoter_before_tss, transcript.start + promoter_after_tss]
            else:
                region = [transcript.end - promoter_after_tss, transcript.end + promoter_before_tss]
            
            # Checking again that promoter is still visible
            if region[1] > start and region[0] < end:
                annotation_regions.append(region)
        
        self.draw_region_annotation_track(
            "Promoters",
            annotation_regions,
            annotation_offset,
            barcolor="brown",
            linewidth=linewidth,
            fontsize=fontsize,
            coordinate_translator=coordinate_translator,
        )
        
        annotation_offset -= annotation_offset_step
        annotation_regions = []
        # Draw enhancers
        for _, enhancer in self.enhancers.enhancers_df.iterrows():
            if enhancer["chr"] == chrom and enhancer["end"] > start and enhancer["start"] < end:
                region = [enhancer["start"], enhancer["end"]]
                annotation_regions.append(region)
        self.draw_region_annotation_track(
            "Enhancers",
            annotation_regions,
            annotation_offset,
            barcolor="darkgreen",
            coordinate_translator=coordinate_translator,
            fontsize=fontsize,
            linewidth=linewidth,
        )
        
        plt.ylim(annotation_offset - annotation_offset_step, ylim[1])
    
    def plot_region_custom_samples(
        self,
        chrom,
        start,
        end,
        custom_sample_fun,
        sample_colors,
        sample_hatch,
        sample_order,
        with_relapse=True,
        with_germline=True,
        ws=10000,
        title=None,
        annotations=[],
        figure_kwargs={},
        must_overlap=None,
        aggregate_samples=False,
        marker_height=0.8,
        min_marker_width_relative=None,
        fontsize=10,
        coordinate_space=True,
        show_has_value_indicator=True,
        segment=None,
        vlines_in_coord_space=None,
        hold=False,
    ):
        if min_marker_width_relative is None:
            if coordinate_space:
                min_marker_width_relative = 0.002
            else:
                min_marker_width_relative = 0
        print(chrom, start, end)
        if ws <= 0:
            ws = end - start
        
        for window_start in range(start, end, ws):
            print("Window: ", window_start, window_start + ws)
            self.pa.figure(**figure_kwargs)
            if title is not None:
                plt.title(title)
            
            merged_matrix = self.loader.get_merged_matrix(
                chrom,
                window_start,
                window_start + ws,
                with_relapse=with_relapse,
                with_germline=with_germline,
                must_overlap_position=must_overlap,
            )
            
            if merged_matrix is None:
                continue
            plot_met_profile(
                np.array(merged_matrix.met_matrix.todense()),
                samples=merged_matrix.read_samples,
                sample_order=sample_order,
                sample_colors=sample_colors,
                sample_hatch=sample_hatch,
                site_genomic_pos=merged_matrix.genomic_coord if coordinate_space else None,
                site_genomic_pos_end=merged_matrix.genomic_coord_end if coordinate_space else None,
                aggregate_samples=aggregate_samples,
                min_marker_width_relative=min_marker_width_relative,
                marker_height=marker_height,
                show_has_value_indicator=show_has_value_indicator,
                segment=segment,
            )
            
            if not coordinate_space:
                
                def translate_to_coordinate_index(x):
                    start = np.where(merged_matrix.genomic_coord <= x)[0]
                    start = start[-1] if len(start) > 0 else 0
                    end = np.where(merged_matrix.genomic_coord_end >= x)[0]
                    end = end[0] if len(end) > 0 else len(merged_matrix.genomic_coord_end) - 1
                    genomic_width = max(1, merged_matrix.genomic_coord[end] - merged_matrix.genomic_coord_end[start])
                    index_width = end - start
                    offset = index_width * ((x - merged_matrix.genomic_coord_end[start]) / genomic_width)
                    return start + offset
                
                if vlines_in_coord_space is not None:
                    vlines = [translate_to_coordinate_index(v) for v in vlines_in_coord_space]
                    plt.vlines(vlines, 0, plt.ylim()[1], linewidth=0.1, color="green")
                
                """ Now plot x-axis in coordinate space """
                x_start, x_end = plt.xlim()
                tick_dist = x_end // 5
                tick_candidates_genomic = np.array(
                    list(range(merged_matrix.genomic_coord[0], merged_matrix.genomic_coord_end[-1]))
                )
                tick_candidates = np.array([translate_to_coordinate_index(x) for x in tick_candidates_genomic])
                ticks = []
                tick_labels = []
                # import pdb
                # pdb.set_trace()
                for tick in np.arange(0, x_end, tick_dist):
                    tick_index = np.argmin(np.abs(tick_candidates - tick))
                    
                    ticks.append(tick)
                    
                    tick_labels.append(f"{tick_candidates_genomic[tick_index]}")
                plt.xticks(ticks, tick_labels)
            else:
                translate_to_coordinate_index = None
            
            self.plot_annotations(
                chrom,
                min(merged_matrix.genomic_coord),
                max(merged_matrix.genomic_coord_end),
                annotations,
                coordinate_translator=translate_to_coordinate_index,
                fontsize=fontsize,
            )
            plt.tick_params(left=False, right=False, labelleft=False)
            plt.xlabel(f"Location on chr{chrom}")
            try:
                plt.ticklabel_format(axis="x", style="plain", useOffset=False)
            except:
                # Will raise an error if ticklabels dont contain proper floats
                pass
            plt.gca().spines["top"].set_visible(False)
            plt.gca().spines["left"].set_visible(False)
            plt.gca().spines["right"].set_visible(False)
            if not hold:
                self.pa.savefig()
    
    def plot_region(self, chrom, start, end, with_germline=True, with_relapse=True, with_no_hp=True, **kwargs):
        sample_colors = self.sample_colors
        sample_hatch = self.sample_hatch
        
        sample_order = []
        if with_germline:
            if with_no_hp:
                sample_order.append("Germline")
            sample_order += ["Germline (HP1)", "Germline (HP2)"]
        if with_no_hp:
            sample_order.append("Primary")
        sample_order += ["Primary (HP1)", "Primary (HP2)"]
        if with_relapse:
            if with_no_hp:
                sample_order.append("Relapse")
            sample_order += ["Relapse (HP1)", "Relapse (HP2)"]
        
        sample_order = sample_order[::-1]
        
        def sample_fun(sample, met_matrix):
            return np.array([f"{sample} (HP{hp})" if hp != -1 else sample for hp in met_matrix.read_samples])
        
        self.plot_region_custom_samples(
            chrom,
            start,
            end,
            sample_fun,
            sample_colors,
            sample_hatch,
            sample_order,
            with_germline=with_germline,
            with_relapse=with_relapse,
            **kwargs,
        )
    
    def plot_promoter(self, gene_id: str, before=5000, after=5000, annotations=[], once_per_annotation=False, **kwargs):
        gene = self.gff.get_gene(gene_id)
        print("Gene: ", gene)
        if gene is None:
            return
        chrom = gene.parent.name
        annotations_not_yet_plotted = {i for i in range(len(annotations))}
        
        transcripts = gene.sorted_children
        regions_covered = []
        for transcript in transcripts:
            direction = transcript.direction
            if direction == "+":
                start = transcript.start - before
                end = transcript.start + after
            elif direction == "-":
                start = transcript.end - after
                end = transcript.end + before
            
            already_plotted = False
            for r in regions_covered:
                if abs(r[0] - start) < before / 2 or abs(r[1] - end) < after / 2:
                    already_plotted = True
            regions_covered.append((start, end))
            if already_plotted:
                continue
            
            if once_per_annotation:
                # plot only if it has an annotation thats not yet there
                has_annotation_in_range = False
                for i, annotation in enumerate(annotations):
                    if i not in annotations_not_yet_plotted:
                        continue
                    if annotation["region"][1] > start and annotation["region"][0] < end:
                        annotations_not_yet_plotted = annotations_not_yet_plotted.difference({i})
                        has_annotation_in_range = True
                if not has_annotation_in_range:
                    continue
            
            self.plot_region(chrom, start, end, title=f"{gene.name}", annotations=annotations, **kwargs)
    
    def plot_gene(self, gene_id: str, before=5000, after=5000, **kwargs):
        gene = self.gff.get_gene(gene_id)
        chrom, start, end = gene.parent.name, gene.start, gene.end
        self.plot_region(chrom, start - before, end + after, **kwargs)
    
    def plot_gene_metrate(self, gene_id: str, before=5000, after=5000, **kwargs):
        gene = self.gff.get_gene(gene_id)
        chrom, start, end = gene.parent.name, gene.start, gene.end
        for sample in module_config.samples:
            chrom_container: ChromosomeContainer = self.h5[sample][chrom]
            chrom_container.get_values_in_range(before, after)
