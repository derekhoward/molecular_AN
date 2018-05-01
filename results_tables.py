import pandas as pd
import data_processing as data
import HBA_analysis as hba
import math
from pathlib import Path

adult_exp = data.get_dataset(dataset='adult', probes_strategy='reannotator')
fetal_exp = data.get_dataset(dataset='fetal', probes_strategy='reannotator')

negraes = data.get_genelist('negraes')
duncan = data.get_genelist('duncan')
lutterAN = data.get_genelist('lutterAN')
lutterBN = data.get_genelist('lutterBN')

results_dir = Path('./results')
results_dir.mkdir(exist_ok=True)

def add_sig_marks(df):
    """adds markers to brain names column: ** for pFDR < 0.05 and * p<0.05"""
    # add ** if FDR < 0.05
    df.loc[df['pFDR'] < 0.05, 'brain area'] = df.loc[df['pFDR']<0.05, 'brain area'].apply(lambda x: str(x)+'**')
    # add * if FDR>0.05 but p<0.05
    df.loc[(df['pFDR'] > 0.05) & (df['p'] < 0.05), 'brain area'] = df.loc[(df['pFDR'] > 0.05) & (df['p'] < 0.05), 'brain area'].apply(lambda x: str(x)+'*')

    return df


def round_sigfigs(num, sig_figs):
#http://code.activestate.com/recipes/578114-round-number-to-specified-number-of-significant-di/
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0 


def process_table(results, brain_areas, sig_figs, filename):
    output = results.loc[brain_areas, :]
    output.index.name = 'brain area'
    output = output.applymap(lambda x: round_sigfigs(x, sig_figs))
    output.reset_index(inplace=True)
    output = add_sig_marks(output)
    output_loc = results_dir / filename
    print(f'Writing results to {output_loc}')
    output.to_csv(output_loc, index=None)


# Results with Negraes gene list
# Adult allen brain data
brain_areas = ['lateral parabrachial nucleus', 'medial parabrachial nucleus',
               'paraventricular nucleus of the hypothalamus', 'arcuate nucleus of the hypothalamus',
               'pontine raphe nucleus', 'subcallosal cingulate gyrus', 'nucleus accumbens', 'ventral tegmental area',
               'central nucleus', 'bed  nucleus of stria terminalis']

results_adult_negraes = hba.generate_stats_table(exp_df=adult_exp, gene_list=negraes)
results_adult_negraes.applymap(lambda x: round_sigfigs(x, 3)).to_csv(results_dir / 'adult_negraes.csv')
results_adult_negraes['Rank'] = results_adult_negraes.AUROC.rank(ascending=False)
process_table(results_adult_negraes, brain_areas, 3, 'adult_negraes_ROI.csv')

# Fetal allen brain data
subgenual_cingulate_cortex = ['IZ in subgenual (subcallosal) cingulate cortex', 'VZ in subgenual cingulate neocortex',
                              'SZ in subgenual cingulate cortex', 'SP in subgenual (subcallosal) cingulate cortex',
                              'outer CP in subgenual (subcallosal) cingulate cortex',
                              'inner CP in subgenual (subcallosal) cingulate cortex']

central_amygdala = ['central nuclear group',
                    'lateral subdivision of central nucleus', # not present in fetal data
                    'medial subdivision of central nucleus']

# generate new columns that are aggregates of smaller areas sampled
fetal_exp['subgenual_cingulate_cortex'] = fetal_exp.loc[:, subgenual_cingulate_cortex].mean(axis=1)
fetal_exp['central_amygdala'] = fetal_exp.loc[:, central_amygdala].mean(axis=1)

# slight differences in naming conventions in comparison to adult...
fetal_brain_areas = ['lateral parabrachial nucleus', 'medial parabrachial nucleus',
                     'paraventricular nucleus of hypothalamus', 'arcuate nucleus of hypothalamus',
                     'raphe magnus nucleus', 'raphe obscurus nucleus', 'solitary nucleus',
                     'core of nucleus accumbens', 'ventral tegmental area',
                     'bed nucleus of stria terminalis', 'subgenual_cingulate_cortex', 'central_amygdala']

results_fetal_negraes = hba.generate_stats_table(exp_df=fetal_exp, gene_list=negraes)
results_fetal_negraes.applymap(lambda x: round_sigfigs(x, 3)).to_csv(results_dir / 'fetal_negraes.csv')
results_fetal_negraes['Rank'] = results_fetal_negraes.AUROC.rank(ascending=False)
process_table(results_fetal_negraes, fetal_brain_areas, 3, 'fetal_negraes_ROI.csv')


# Results with Duncan gene list
# - repeat same procedure but using rpy2 to generate table since there are only 6 genes of interest

# Adult brain data

results_adult_duncan = hba.generate_Rstats_table(exp_df=adult_exp, gene_list=duncan)
process_table(results_adult_duncan, brain_areas, 3, 'adult_duncan_ROI.csv')

# Fetal brain data

results_fetal_duncan = hba.generate_Rstats_table(exp_df=fetal_exp, gene_list=duncan)
process_table(results_fetal_duncan, fetal_brain_areas, 3,  'fetal_duncan_ROI.csv')

# Results with LutterAN and LutterBN gene lists

# Adult data

results_adult_lutterAN = hba.generate_stats_table(exp_df=adult_exp, gene_list=lutterAN)
results_adult_lutterBN = hba.generate_stats_table(exp_df=adult_exp, gene_list=lutterBN)

process_table(results_adult_lutterAN, brain_areas, 3, 'adult_lutterAN_ROI.csv')
process_table(results_adult_lutterBN, brain_areas, 3, 'adult_lutterBN_ROI.csv')

# Fetal expression data
results_fetal_lutterAN = hba.generate_stats_table(exp_df=fetal_exp, gene_list=lutterAN)
results_fetal_lutterBN = hba.generate_stats_table(exp_df=fetal_exp, gene_list=lutterBN)

process_table(results_fetal_lutterAN, fetal_brain_areas, 3, 'fetal_lutterAN_ROI.csv')
process_table(results_fetal_lutterBN, fetal_brain_areas, 3, 'fetal_lutterBN_ROI.csv')
