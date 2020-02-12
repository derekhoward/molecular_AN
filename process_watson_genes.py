from pathlib import Path
import pandas as pd


GENELISTS_OUTPUT = Path('./data/genelists')
sig_loci = ['NCKIPSD', 'CADM1', 'ASB3', 'ERLEC1',
            'MGMT', 'FOXP1', 'PTBP2', 'CDH10', 'NSUN3']

print(
    f'Writing 9 genes near the 8 significant loci to {GENELISTS_OUTPUT / "Watson et al.Table1.SIG_LOCI.txt"}')
pd.Series(sig_loci).to_csv(GENELISTS_OUTPUT /
                           'Watson et al.Table1.SIG_LOCI.txt', index=False, header=False)

protein_coding = pd.read_excel(
    './data/raw/41588_2019_439_MOESM3_ESM.xlsx', sheet_name=5, skiprows=1)


print(
    f'Writing protein coding genes to {GENELISTS_OUTPUT / "Watson et al.TableS6.protein_coding.ALL.txt"}')
protein_coding.GENENAME_GM.dropna().to_csv(
    GENELISTS_OUTPUT / 'Watson et al.TableS6.protein_coding.ALL.txt', index=False, header=False)


top_loc = protein_coding[protein_coding.INPUTID == 1].loc[:, 'GENENAME_GM']
print(
    f'Writing protein coding genes near the first SNP of interest to {GENELISTS_OUTPUT / "Watson et al.TableS6.protein_coding.loc1.txt"}')
top_loc.dropna().to_csv(
    GENELISTS_OUTPUT / 'Watson et al.TableS6.protein_coding.loc1.txt', index=False, header=False)


remaining_genes = protein_coding[protein_coding.INPUTID !=
                                 1].loc[:, 'GENENAME_GM']
print(
    f'Writing protein coding genes that are not associated to top SNP {GENELISTS_OUTPUT / "Watson et al.TableS6.protein_coding.remaining_loci.txt"}')
remaining_genes.dropna().to_csv(
    GENELISTS_OUTPUT / "Watson et al.TableS6.protein_coding.remaining_loci.txt", index=False, header=False)


magma_genes = pd.read_excel(
    './data/raw/41588_2019_439_MOESM3_ESM.xlsx', sheet_name=10, skiprows=1)
magma_genes = magma_genes[magma_genes['P (green =  Bonferroni-significant)']
                          < 0.0000025]
magma_genes = magma_genes.Name
print(
    f'Writing genes derived from genewise MAGMA analysis {GENELISTS_OUTPUT / "Watson et al.TableS11.genewiseMAGMA.txt"}')
magma_genes.dropna().to_csv(
    GENELISTS_OUTPUT / "Watson et al.TableS11.genewiseMAGMA.txt", index=False, header=False)
