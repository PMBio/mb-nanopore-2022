import numpy as np
import matplotlib.pyplot as plt

from .general_plotting import set_figure_defaults


def compare_models(nome_model, other_model, include_unchanged=True, include_untrained=False, other_model_also_nome=False, fixed_bins=True, title=None):
    nome_model = nome_model.copy() # working on a copy

    nome_model['nome_kmer'] = nome_model.index
    if other_model_also_nome:
        nome_model['legacy_kmer'] = nome_model['nome_kmer']
    else:
        nome_model['legacy_kmer'] = nome_model['nome_kmer'].map(lambda x: x.replace('Q','T').replace('X','A').replace('Z','C'))

    merged = nome_model.merge(other_model, left_on='legacy_kmer', right_index=True, how='left')

    if not include_unchanged:
        merged = merged.loc[(merged['level_mean_x'] - merged['level_mean_y'])!=0]
    if not include_untrained:
        if 'was_trained_x' in merged.columns:
            merged = merged.loc[(merged['was_trained_x'] & merged['was_trained_y'])]
        elif 'was_trained' in merged.columns:
            merged = merged.loc[merged['was_trained']]

    merged.index = merged.apply(lambda x: '%s_%s'%(x.index, x['legacy_kmer']), axis=1)
    print(merged.shape)

    def plot_hist(idx,ax):
        if fixed_bins:
            bins = np.arange(-15,15,0.5)
            ax.hist(merged.loc[idx]['level_mean_x'] - merged.loc[idx]['level_mean_y'], bins=bins)
        else:
            ax.hist(merged.loc[idx]['level_mean_x'] - merged.loc[idx]['level_mean_y'])


    idx_5mc = np.array(merged.nome_kmer.map(lambda x: x.find('Z') >= 0))
    idx_6ma = np.array(merged.nome_kmer.map(lambda x: x.find('X') >= 0 or x.find('Q') >= 0))
    idx_n5mc = np.array(merged.nome_kmer.map(lambda x: x.find('Z') == -1))
    idx_n6ma = np.array(merged.nome_kmer.map(lambda x: x.find('X') == -1 and x.find('Q') == -1))

    idx_5mc_n6ma = idx_5mc & idx_n6ma
    idx_n5mc_6ma = idx_n5mc & idx_6ma
    idx_n5mc_n6ma = idx_n5mc & idx_n6ma
    idx_5mc_6ma = idx_5mc & idx_6ma

    fig, ax = plt.subplots(nrows=3,ncols=3,figsize=(15,10))
    fig.patch.set_facecolor("w")
    if not title is None:
        fig.suptitle(title, fontsize=16)
    plot_hist(idx_n5mc_n6ma, ax[0,0])
    plot_hist(idx_5mc_n6ma, ax[0,1])
    plot_hist(idx_n6ma, ax[0,2])
    plot_hist(idx_n5mc_6ma, ax[1,0])
    plot_hist(idx_5mc_6ma, ax[1,1])
    plot_hist(idx_6ma, ax[1,2])
    plot_hist(idx_n5mc, ax[2,0])
    plot_hist(idx_5mc, ax[2,1])
    fig.tight_layout(rect=[0, 0.07, 1, 0.93])

    cols = ['No 5mc', '5mc', '']
    rows = ['No 6ma', '6ma', '']

    for axi, col in zip(ax[0], cols):
        axi.set_title(col)

    for axi, row in zip(ax[:,0], rows):
        axi.set_ylabel(row, rotation=90, size='large')