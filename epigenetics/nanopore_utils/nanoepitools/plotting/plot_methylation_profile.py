import queue
from typing import List, Tuple, Union, Optional, Dict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch, Circle
from matplotlib.transforms import ScaledTranslation

from meth5.sparse_matrix import SparseMethylationMatrixContainer
from nanoepitools.math import p_to_llr


def llr_to_color_value(val):
    return 1 - np.exp(-np.abs(val) * 0.5)


def compute_aggregate_matrix(matrix, threshold=2.5):
    matrix = (matrix > threshold).sum(axis=0) / (np.abs(matrix) > threshold).sum(axis=0)
    matrix = p_to_llr(matrix)
    matrix[np.isnan(matrix)] = 0
    matrix = matrix[np.newaxis, :]
    return matrix


def default_color_map(llr):
    if llr < 0:
        # A slightly matte blue
        return (24 / 255, 3 / 255, 248 / 255, llr_to_color_value(-llr))
    else:
        # A slightly matte orange
        return (252 / 255, 127 / 255, 44 / 255, llr_to_color_value(llr))


def plot_vertical_highlights(
    highlights: List[Tuple], highlight_color="b", site_genomic_pos_start=None, highlights_in_genomic_space=False
):
    for i, highlight_range in enumerate(highlights):
        if site_genomic_pos_start is not None and not highlights_in_genomic_space:
            highlight_range = [site_genomic_pos_start[p] for p in highlight_range]
        
        if isinstance(highlight_color, list):
            color = highlight_color[i]
        else:
            color = highlight_color
        
        plt.axvspan(
            highlight_range[0],
            highlight_range[1],
            color=color,
            alpha=0.25,
        )


def plot_segment_lines(segment: np.array, site_genomic_pos_start=None):
    ymax = plt.gca().get_ylim()[1]
    for i in range(1, len(segment)):
        if segment[i] == segment[i - 1]:
            continue
        if site_genomic_pos_start is not None:
            x = (site_genomic_pos_start[i - 1] + site_genomic_pos_start[i]) / 2
        else:
            x = i - 1 + 0.5
        plt.plot((x, x), (0, ymax - 1), c="k")


class DelayedExecution:
    def __init__(self):
        self.delayed_functions = []
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        for fn, args, kwargs in self.delayed_functions:
            fn(*args, **kwargs)
    
    def execute_on_exit(self, fn, *args, **kwargs):
        self.delayed_functions.append((fn, args, kwargs))


class MethylationPlot:
    def __init__(
        self,
        marker_height: float = 0.75,
        min_marker_width_relative: float = 0.002,
        sample_colors: Optional[Dict] = None,
        sample_hatch: Optional[Dict] = None,
        aggregate_samples: bool = False,
        color_map=default_color_map,
        show_has_value_indicator=True
    ):
        self.marker_height = marker_height
        self.min_marker_width_relative = min_marker_width_relative
        self.sample_colors = sample_colors
        self.sample_hatch = sample_hatch
        self.aggregate_samples = aggregate_samples
        self.color_map = color_map
        self.show_has_value_indicator = show_has_value_indicator
    
    def plot_legend(self, sample_marker_range, sample_order):
        if self.sample_colors is not None:
            xlim = plt.xlim()
            marker_width = (xlim[1] - xlim[0]) * 0.015
            for s, marker in sample_marker_range.items():
                kwargs = {}
                if self.sample_hatch is not None:
                    kwargs.update({"hatch": self.sample_hatch[s], "edgecolor": "w"})
                
                patch = Rectangle(
                    (xlim[1] - marker_width * 2, marker[0]),
                    marker_width,
                    marker[1] - marker[0],
                    facecolor=self.sample_colors[s],
                    **kwargs,
                )
                plt.gca().add_patch(patch)
            if self.sample_hatch is not None:
                legend_elements = [
                    Patch(facecolor=self.sample_colors[s], edgecolor="w", hatch=self.sample_hatch[s], label=s)
                    for s in sample_order
                ]
            else:
                legend_elements = [Patch(facecolor=self.sample_colors[s], edgecolor="w", label=s) for s in sample_order]
            
            plt.legend(handles=legend_elements[::-1], bbox_to_anchor=(1.05, 1), loc="upper left")
    
    def plot_met_patches(self, x, x_end, y, color):
        patches = [Rectangle((x[i], y[i]), x_end[i] - x[i], self.marker_height) for i in range(len(x))]
        patch_collection = PatchCollection(patches)
        patch_collection.set_color(color)
        patch_collection.set_edgecolor(None)
        plt.gca().add_collection(patch_collection)
    
    def plot_has_value_indicator(self, x_start, x_end, y):
        if not self.show_has_value_indicator:
            return
        x = (x_end + x_start) / 2
        y = y + self.marker_height / 2
        radius = 0.5
        patches = [
            Circle(
                (0, 0),
                radius=radius,
                transform=ScaledTranslation(x[i], y[i], plt.gca().transData),
                edgecolor=None,
                facecolor="k",
            )
            for i in range(len(x))
        ]
        for patch in patches:
            plt.gca().add_patch(patch)
    
    def plot_met_profile(
        self,
        matrix: np.ndarray,
        samples: Optional[np.ndarray] = None,
        site_genomic_pos_start: Optional[np.ndarray] = None,
        site_genomic_pos_end: Optional[np.ndarray] = None,
        sample_order: Optional[List[str]] = None,
    ):
        if samples is None:
            samples = np.array(["_" for _ in range(matrix.shape[0])])
            sample_order = ["_"]
        
        if sample_order is None:
            sample_order = sorted(list(set(samples)))
        
        if site_genomic_pos_end is not None:
            x_range = np.max(site_genomic_pos_end) - np.min(site_genomic_pos_start)
            min_marker_width = self.min_marker_width_relative * x_range
        
        y_off = 0
        start = 0
        x_max = 0
        end = matrix.shape[1]
        sample_marker_range = {}
        
        with DelayedExecution() as delayed:
            for s in sample_order:
                x = np.arange(start, end)
                part_matrix = matrix[:, x][(samples == s)]
                
                if part_matrix.shape[0] <= 1:
                    continue
                
                if self.aggregate_samples:
                    part_matrix = compute_aggregate_matrix(part_matrix)
                
                active_reads = np.array((part_matrix != 0).sum(axis=1)).flatten() > 0
                
                part_matrix = part_matrix[active_reads]
                hasval = np.array(part_matrix != 0).flatten()
                y = np.arange(part_matrix.shape[0]) + y_off
                
                x, y = np.meshgrid(x, y)
                x = x.flatten()[hasval]
                y = y.flatten()[hasval]
                matrix_data = np.array(part_matrix).flatten()[hasval]
                
                color = [self.color_map(v) for v in matrix_data]
                
                if site_genomic_pos_start is not None:
                    if site_genomic_pos_end is not None:
                        x_end = site_genomic_pos_end[x]
                    x = site_genomic_pos_start[x]  # Translate to actual pos on chrom
                    if site_genomic_pos_end is not None:
                        # Makes it so very short blocks are still visible
                        marker_adjust = (min_marker_width - x_end + x) / 2
                        marker_adjust[marker_adjust < 0] = 0
                        x = x - marker_adjust
                        x_end = x_end + marker_adjust
                else:
                    x = x - 0.4
                    x_end = x + 0.8
                
                self.plot_met_patches(x, x_end, y, color)
                plt.gca().autoscale_view()
                
                # The has_value indicators don't look good in coordinate space, so for now
                # only plot them if we are not in coordinate space
                if site_genomic_pos_start is None:
                    # Executing this delayed, because it needs to know the x and y limits,
                    # and at the moment we are still adding stuff to the plot that will extend
                    # the limits
                    delayed.execute_on_exit(self.plot_has_value_indicator, x, x_end, y)
                
                if self.sample_colors is not None:
                    sample_marker_range[s] = (y_off, part_matrix.shape[0] + y_off)
                    x_max = max(x.max(), x_max)
                
                y_off += part_matrix.shape[0]
        
        self.plot_legend(sample_marker_range, sample_order)
    
    def plot_met_profile_from_matrix(self, matrix: SparseMethylationMatrixContainer, **kwargs):
        plot_met_profile(
            self,
            matrix.met_matrix.todense(),
            samples=matrix.read_samples,
            site_genomic_pos=matrix.genomic_coord,
            site_genomic_pos_end=matrix.genomic_coord_end,
            **kwargs,
        )


def plot_met_profile(
    matrix: np.ndarray,
    samples: np.ndarray = None,
    sample_order: List[str] = None,
    sample_colors: Dict = None,
    site_genomic_pos=None,
    site_genomic_pos_end=None,
    marker_height=0.75,
    segment: np.array = None,
    segments_in_coord_space=False,
    highlights: List[Tuple] = None,
    highlight_color: Union[str, List[str]] = None,
    highlights_in_genomic_space: bool = False,
    min_marker_width_relative: float = 0.002,
    sample_hatch: Dict = None,
    aggregate_samples=False,
    show_has_value_indicator=True,
):
    """Deprecated but here for backwards compatibility"""
    p = MethylationPlot(
        marker_height=marker_height,
        min_marker_width_relative=min_marker_width_relative,
        sample_colors=sample_colors,
        sample_hatch=sample_hatch,
        aggregate_samples=aggregate_samples,
        color_map=default_color_map,
        show_has_value_indicator=show_has_value_indicator,
    )
    p.plot_met_profile(
        matrix=matrix,
        samples=samples,
        site_genomic_pos_start=site_genomic_pos,
        site_genomic_pos_end=site_genomic_pos_end,
        sample_order=sample_order,
    )
    if highlights is not None:
        plot_vertical_highlights(
            highlights=highlights,
            highlight_color=highlight_color,
            site_genomic_pos_start=site_genomic_pos,
            highlights_in_genomic_space=highlights_in_genomic_space,
        )
    if segment is not None:
        plot_segment_lines(
            segment=segment, site_genomic_pos_start=site_genomic_pos, segments_in_coord_space=segments_in_coord_space
        )


def plot_met_profile_from_matrix(matrix: SparseMethylationMatrixContainer, **kwargs):
    """Deprecated but here for backwards compatibility"""
    plot_met_profile(
        matrix.met_matrix.todense(),
        samples=matrix.read_samples,
        site_genomic_pos=matrix.genomic_coord,
        site_genomic_pos_end=matrix.genomic_coord_end,
        **kwargs,
    )
