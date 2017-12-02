import GEOparse
import pandas as pd

gse2 = GEOparse.get_GEO(filepath="./Data/Human/GSE26637_family.soft.gz")

plats_2 = list(gse2.gpls.keys())[0]

samples2 = gse2.phenotype_data[["characteristics_ch1.0.gender", "characteristics_ch1.2.stimulation", "characteristics_ch1.3.resistance status"]]
samples2 = samples2.rename(columns={'characteristics_ch1.0.gender':'gender', 'characteristics_ch1.2.stimulation':'fasting_status',
                         'characteristics_ch1.3.resistance status':'insulin_status'})
samples2['cbmi'] = samples2['insulin_status'].map(lambda x: 'lean' if x == 'sensitive' else 'obese')

samples2.to_pickle('./Preprocessed_Data/Human/batch2_pheno.p')
with open('./Preprocessed_Data/Human/batch2_pheno.txt', 'w') as handle:
    samples2.to_csv(handle, sep='\t')

samples2_exprs = gse2.pivot_samples('VALUE')[list(samples2.index)]

samples2_ann = samples2_exprs.reset_index().merge(gse2.gpls['GPL570'].table[["ID", "Gene Symbol"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
samples2_ann.drop('ID', inplace=True, axis=1)
samples2_ann['Gene Symbol'] = samples2_ann['Gene Symbol'].astype(str)

samples2_ann = samples2_ann[~samples2_ann["Gene Symbol"].str.contains("///")].dropna()
samples2_ann['Gene Symbol'].astype(str, inplace=True)
samples2_ann = samples2_ann.dropna()

exs_2 = samples2_ann.groupby('Gene Symbol').mean()

del exs_2.index.name
exprs_2 = exs_2.T
exprs_2.to_pickle('./Preprocessed_Data/Human/batch2_geno.p')
with open('./Preprocessed_Data/Human/batch2_geno.txt', 'w') as handle:
    exprs_2.to_csv(handle, sep='\t')
