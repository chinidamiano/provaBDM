import GEOparse
import pandas as pd

gse4 = GEOparse.get_GEO(filepath="./Data/Human/GSE48964_family.soft.gz")

plats_4 = list(gse4.gpls.keys())[0]

samples4 = gse4.phenotype_data["source_name_ch1"]
samples4 = pd.DataFrame(samples4); samples4.head()
samples4.rename(columns={"source_name_ch1":"cbmi"}, inplace=True); samples4.head()
samples4["cbmi"] = samples4["cbmi"].apply(lambda x: 'obese' if x.split(' ')[1].lower()=='obese' else 'lean'); samples4.head()

samples4.to_pickle('./Preprocessed_Data/Human/batch4_pheno.p')
with open('./Preprocessed_Data/Human/batch4_pheno.txt', 'w') as handle:
    samples4.to_csv(handle, sep='\t')

samples4_exprs = gse4.pivot_samples('VALUE')[list(samples4.index)]
samples4_exprs.index.name = 'ID'

ann_table4 = pd.read_csv('./Data/Human/GPL6244.annot', sep='\t', skiprows=27, dtype=str, na_values='NaN',usecols=[0,2], index_col=0)#.dropna()

samples4_exprs = samples4_exprs.dropna()

ann_table = ann_table4.dropna()

ann_table = ann_table[~ann_table['Gene symbol'].str.contains("///")]

samples4_exprs.index = samples4_exprs.index.astype(str)

exprs_4 = ann_table.merge(samples4_exprs, left_index=True, right_index=True, how='inner')

exprs_4 = exprs_4.groupby('Gene symbol').mean()
del exprs_4.index.name

exprs_4 = exprs_4.T
exprs_4.to_pickle('./Preprocessed_Data/Human/batch4_geno.p')
with open('./Preprocessed_Data/Human/batch4_geno.txt', 'w') as handle:
    exprs_4.to_csv(handle, sep='\t')
