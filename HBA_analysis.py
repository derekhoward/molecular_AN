import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.sandbox.stats.multicomp
from sklearn import metrics
import rpy2
import rpy2.robjects

def calc_pvalues(exp_series, gene_list):
    X = exp_series[exp_series.index.isin(gene_list)]
    y = exp_series[~exp_series.index.isin(gene_list)]
    assert(X.shape[0] + y.shape[0] == exp_series.shape[0])

    return mannwhitneyu(X, y, alternative='two-sided')[1]


def calc_AUC(exp_series, gene_list):
    # simplified function that calcs AUC for a single brain area series
    y_true = exp_series.index.isin(gene_list)
    y_score = exp_series.values

    return metrics.roc_auc_score(y_true, y_score)


def r_mannwhitneyu(sample1, sample2, exact=True, alternative="g"):
    sample1 = "c({})".format(str(list(sample1))[1:-1])
    sample2 = "c({})".format(str(list(sample2))[1:-1])
    if exact:
        exactness = ", exact=TRUE"
    else:
        exactness = ""
    rpy2.robjects.R()("""wres <- wilcox.test({}, {}, alternative="{}"{});
                       rm(sample1);
                       rm(sample2);""".format(sample1, sample2, alternative, exactness))
    wres = rpy2.robjects.r['wres']
    uval = wres[0][0]
    pval = wres[2][0]
    return uval, pval


def calc_R_pvalues(exp_series, gene_list):
    X = exp_series[exp_series.index.isin(gene_list)]
    y = exp_series[~exp_series.index.isin(gene_list)]
    assert(X.shape[0] + y.shape[0] == exp_series.shape[0])

    return r_mannwhitneyu(X, y)[1]


def generate_stats_table(exp_df, gene_list, verbose=True):
    """
    Creates a table of summary stats for each brain area

    Parameters
    ----------
    exp_df : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list : series
        list of gene symbols of interest
    Returns
    -------
    table
        results dataframe with brain areas as index

    """
    count = len(gene_list)
    n_genes_in_matrix = gene_list.isin(exp_df.index).sum()
    genes_not_found = gene_list[~gene_list.isin(exp_df.index)].values

    if verbose:
        print('You submitted a gene list with {} genes.\n\
    {} of those genes are present in the reference dataset.\n\
    Genes not found in our reference data: {}'.format(
                count, n_genes_in_matrix, genes_not_found))

    pvalues = exp_df.apply(lambda col: calc_pvalues(exp_series=col, gene_list=gene_list))
    fdr_corrected = statsmodels.sandbox.stats.multicomp.multipletests(pvalues, method="fdr_bh")[1]
    fdr_corrected = pd.Series(fdr_corrected, index=pvalues.index)
    auc = exp_df.apply(lambda col: calc_AUC(exp_series=col, gene_list=gene_list))

    table = pd.concat([auc, pvalues, fdr_corrected],
                      keys=['AUROC', 'p', 'pFDR', ], axis=1)
    table.set_index(exp_df.columns, inplace=True)
    return table.sort_values('AUROC', ascending=False)


def generate_Rstats_table(exp_df, gene_list, verbose=False):
    """
    Creates a table of summary stats for each brain area using R's wilcox test
    to account for approximation of python's mann whitney u

    Parameters
    ----------
    exp_df : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list : series
        list of gene symbols of interest
    Returns
    -------
    table
        dataframe with brain areas as index

    """
    count = len(gene_list)
    n_genes_in_matrix = gene_list.isin(exp_df.index).sum()
    genes_not_found = gene_list[~gene_list.isin(exp_df.index)].values

    if verbose:
        print('You submitted a gene list with {} genes.\n\
    {} of those genes are present in the reference dataset.\n\
    Genes not found in our reference data: {}'.format(
                count, n_genes_in_matrix, genes_not_found))

    pvalues = exp_df.apply(lambda col: calc_R_pvalues(exp_series=col, gene_list=gene_list))
    fdr_corrected = statsmodels.sandbox.stats.multicomp.multipletests(pvalues, method="fdr_bh")[1]
    fdr_corrected = pd.Series(fdr_corrected)
    fdr_corrected.index = pvalues.index
    auc = exp_df.apply(lambda col: calc_AUC(exp_series=col, gene_list=gene_list))

    table = pd.concat([auc, pvalues, fdr_corrected],
                      keys=['AUROC', 'p', 'pFDR'], axis=1)
    table.set_index(exp_df.columns, inplace=True)
    return table.sort_values('AUROC', ascending=False)
