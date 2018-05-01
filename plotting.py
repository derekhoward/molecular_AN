import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats.kde import gaussian_kde
import numpy as np
import itertools


def make_dist_plot(exp_matrix, brain_area, gene_list, ax, gene_text=False):
    """
    Plots the distributions of zscored levels of expression of genes of
    interest compared to the rest of genes

    Parameters
    ----------
    exp_matrix : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    brain_area : str
        select brain area to be plotted
    gene_list : series
        list of gene symbols of interest
    Returns
    -------

    """
    brain_area_exp = exp_matrix.loc[:, brain_area].reset_index()

    # get the data for genes of interest
    brain_area_exp['In gene list'] = brain_area_exp.gene_symbol.isin(gene_list)
    # select gene names/exp vals from rows with genes from gene_list
    exp_vals = list(brain_area_exp.loc[brain_area_exp['In gene list'], brain_area])
    gene_names = list(brain_area_exp.loc[brain_area_exp['In gene list'], 'gene_symbol'])

    background_vals = brain_area_exp[~brain_area_exp['In gene list']].loc[:, brain_area]

    kde = gaussian_kde(background_vals)

    # these are the values over which your kernel will be evaluated
    dist_space = np.linspace(min(background_vals), max(background_vals), 100)

    # plot the results
    with sns.plotting_context('notebook', font_scale = 1.25):
        distplot = ax.plot(dist_space, kde(dist_space), color='k', lw=1.5)
        #ax.set_xlabel('Expression (z-scored)')
        #distplot.set_xlabel('Expression (z-scored)')
        ax.set_ylabel('Density')
        #distplot.set_ylabel('Density')
        #ax.set_title('Expression of 6 gene hitlist within {}'.format(brain_area))
        #ax.set_title('{} in '.format(brain_area.capitalize()))

        # plot lines over the distribution for the disease genes
        palette = itertools.cycle(sns.color_palette())
        i = 0

        for name, position, colour in zip(gene_names, exp_vals, palette):
            ax.vlines(position, ymin=-0.1, ymax=1, lw=1.5, colors=colour)
            #plt.text(position, 0.75, name, rotation=90, fontsize=12, verticalalignment='center')
            if gene_text is True:
                #plt.text(6, 5.25-(i/10), name, fontsize=16, color=colour)
                plt.text(6, 6.75-(i/8), name, fontsize=16, color=colour)
            i+=1
    sns.despine()
    return distplot


def make_violins(exp, brain_areas, gene_list, brain_area_labels=None):
    subset = exp.loc[:, brain_areas]
    subset['in_gene_list'] = subset.index.isin(gene_list)
    tidy = subset.reset_index().melt(id_vars=['gene_symbol', 'in_gene_list'],
                                     var_name='brain_area',
                                     value_name='expression')
    tidy['in_gene_list'] = tidy['in_gene_list'].map({True: 'AN genes',
                                                     False: 'Background genes'})
    with sns.plotting_context('notebook', font_scale=1.25):
        fig, ax = plt.subplots(figsize=(12, 12))
        sns.violinplot(y='brain_area', x='expression', hue='in_gene_list',
                       palette={'AN genes': '#D9D4D3', 'Background genes': '#636364'},
                       split=True, hue_order=['AN genes', 'Background genes'],
                       inner='quartiles', data=tidy, ax=ax)
        ax.set_xlabel('Expression (z-scored)')
        ax.set_ylabel('')
        if brain_area_labels:
            ax.set_yticklabels(brain_area_labels)
        legend = ax.get_legend()
        legend.set_title('')
        legend._loc = 7

        sns.despine()
