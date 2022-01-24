from pathlib import Path
from typing import List

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

from nanoepitools.plotting.syncbot import SyncBot


def set_if_not_in_dict(d, k, v):
    if k not in d.keys():
        d[k] = v


def set_figure_defaults(fig):
    fig.set_tight_layout(True)
    fig.patch.set_facecolor("w")
    fig.autolayout = False


def plot_1d_density(x, y, *argc, **kwargs):
    density = scipy.stats.gaussian_kde(y)
    plt.plot(x, density(x), *argc, **kwargs)


def plot_2d_density(x, y, nbins=50, cmap=plt.cm.BuGn_r, contour_colors="k"):
    k = scipy.stats.gaussian_kde((x, y))
    xi, yi = np.mgrid[x.min() : x.max() : nbins * 1j, y.min() : y.max() : nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading="gouraud", cmap=cmap)
    plt.contour(xi, yi, zi.reshape(xi.shape), colors=contour_colors)


def plot_multiple_histograms(
    data: List[np.ndarray],
    bins=200,
    bound=np.inf,
    ubound=np.inf,
    lbound=-np.inf,
    alpha=0.5,
    colors=None,
    labels=None,
    normalize_histograms=False,
    title="",
    xlabel="",
    ylabel="Frequency",
):
    if not np.isinf(bound):
        ubound = np.abs(bound)
        lbound = -np.abs(bound)
    
    if not np.isinf(ubound) and not np.isinf(lbound):
        bins = np.arange(-bound, bound + 1, (2 * (bound + 1)) / bins)
    
    for i, val in enumerate(data):
        val = np.clip(val, lbound, ubound)
        weights = np.ones(len(val))
        if normalize_histograms:
            weights = weights / len(val)
        plt.hist(
            val,
            weights=weights,
            bins=bins,
            alpha=alpha,
            color=(colors[i] if colors is not None else None),
            label=(labels[i] if labels is not None else None),
        )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()


class PDFPagesWrapper:
    def __init__(self, pa, *argv):
        self.pa = pa
        self.pdf = PdfPages(*argv)
    
    def __enter__(self):
        self.pdf.__enter__()
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pa.pdf = None
        self.pdf.__exit__(exc_type, exc_val, exc_tb)
        self.pa.syncplots()


class PlotArchiver:
    def __init__(self, project, config=None):
        matplotlib.use("Agg")
        self.project = project
        if config is None:
            config = dict()
        else:
            config = dict(config) # make a copy
        set_if_not_in_dict(config, "plot_archive_dir", Path.home().joinpath("nanoepitools_plots"))
        set_if_not_in_dict(config, "syncplot_config_file", None)  # let syncbot choose the default
        set_if_not_in_dict(config, "filetype", "pdf")
        self.config = config
        
        self.project_path = Path(self.config["plot_archive_dir"]).joinpath(self.project)
        self.pdf: PDFPagesWrapper = None
        self.syncbot = SyncBot(config["syncplot_config_file"])
    
    def ensure_project_path_exists(self):
        self.project_path.mkdir(parents=True, exist_ok=True)
    
    def get_plot_path(self, key, filetype=None):
        if filetype is None:
            filetype = self.config["filetype"]
        filename = "{key}.{ft}".format(key=key, ft=filetype)
        return self.project_path.joinpath(filename)
    
    def savefig(self, key="figure", fig=None, close=True):
        self.ensure_project_path_exists()
        legends = [c for c in plt.gca().get_children() if isinstance(c, matplotlib.legend.Legend)]
        if fig is None:
            fig = plt.gcf()
        if self.pdf is not None:
            self.pdf.pdf.savefig(fig, bbox_inches="tight", bbox_extra_artists=legends)
        else:
            path = self.get_plot_path(key)
            fig.savefig(path, bbox_inches="tight", bbox_extra_artists=legends)
            self.syncplots()
        if close:
            plt.close(fig)
    
    def saveandshow(self, key="figure", fig=None):
        if fig is None:
            fig = plt.gcf()
        self.savefig(key, fig=fig, close=False)
        fig.show()
        # This is a workaround needed to ensure the figure is actually
        # displayed in jupyter notebook. If we don't wait for a bit and close
        # it right away, it may not be displayed
        plt.pause(0.0001)
        plt.close(fig)
    
    def open_multipage_pdf(self, key):
        self.ensure_project_path_exists()
        path = self.get_plot_path(key, "pdf")
        self.pdf = PDFPagesWrapper(self, path)
        return self.pdf
    
    def subplots(self, *args, **kwargs):
        if "dpi" not in kwargs:
            kwargs["dpi"] = 200
        fig, axes = plt.subplots(*args, **kwargs)
        set_figure_defaults(fig)
        return fig, axes
    
    def syncplots(self):
        self.syncbot.sync()
    
    def figure(self, **kwargs):
        """
        simply calls matplotlib.pyplot.figure() but then sets some default
        parameters that I find handy
        :param argc: parameters for plt.figure
        :returns: return value if plt.figure
        """
        if "dpi" not in kwargs:
            kwargs["dpi"] = 200
        fig = plt.figure(**kwargs)
        set_figure_defaults(fig)
        return fig
