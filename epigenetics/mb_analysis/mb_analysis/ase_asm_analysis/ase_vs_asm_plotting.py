import matplotlib.pyplot as plt
import numpy as np

from nanoepitools.plotting.general_plotting import PlotArchiver


def plot_ase_vs_asm(
    pa: PlotArchiver,
    ase,
    asm,
    labels,
    colors=None,
    colorbar=False,
    hp2_minus_hp1=False,
    shape=None,
    min_abs_diff=0.5,
    **kwargs,
):
    """ Plot ASE vs ASM with a broken X-axis"""
    f, (ax1, ax2) = pa.subplots(1, 2, sharey=True)
    # for type in {"promoter", "genebody"}:
    x = np.array(asm)
    if hp2_minus_hp1:
        x = -x
    y = np.array(ase)
    for xval, yval, label in zip(x, y, labels):
        ax = ax1 if xval < 0 else ax2
        ax.text(xval + 0.01, yval, label, fontsize=10, va="center")
    
    if colors is None:
        colors = ["b" for _ in x]
    colors = np.array(colors)
    
    if shape is None:
        shape = ["o" for _ in x]
    shape = np.array(shape)
    
    for marker in set(shape[x < 0]):
        idx = (x < 0) & (shape == marker)
        ax1.scatter(x[idx], y[idx], c=colors[idx], marker=marker, s=14, **kwargs)
    for marker in set(shape[x > 0]):
        idx = (x > 0) & (shape == marker)
        sc = ax2.scatter(x[idx], y[idx], c=colors[idx], marker=marker, s=14, **kwargs)
    
    if colorbar:
        f.colorbar(sc)
    
    ax1.set_xlim(-1, -min_abs_diff)
    ax2.set_xlim(min_abs_diff, 1)
    ax1.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    plt.ylim(0, 1)
    plt.xlabel(f"Differential methylation (r_HP1 - r_HP2)")
    ax1.set_ylabel("Allele specific expression ratio")


def plot_ase_vs_asm_nosplit(
    pa: PlotArchiver,
    ase,
    asm,
    labels,
    colors=None,
    colorbar=False,
    hp2_minus_hp1=False,
    shape=None,
    min_abs_diff=None,
    **kwargs,
):
    """ Plot ASE vs ASM with continuous X-axis"""
    pa.figure()
    x = np.array(asm)
    if hp2_minus_hp1:
        x = -x
    y = np.array(ase)
    for xval, yval, label in zip(x, y, labels):
        plt.text(xval + 0.01, yval, label, fontsize=10, va="center", fontname="Arial")
    
    if colors is None:
        colors = ["b" for _ in x]
    colors = np.array(colors)
    
    if shape is None:
        shape = ["o" for _ in x]
    shape = np.array(shape)
    
    for marker in set(shape):
        idx = shape == marker
        plt.scatter(x[idx], y[idx], c=colors[idx], marker="o", linewidths=0.5 if marker == "*" else 0, edgecolors="k", s=14 if marker == "*" else 6, **kwargs)
    plt.xlim(-1, 1)
    plt.ylim(0, 1)
    plt.gca().set_aspect(2)
    plt.xlabel(f"Differential methylation (r_HP1 - r_HP2)")
    plt.ylabel("Allele specific expression ratio")
