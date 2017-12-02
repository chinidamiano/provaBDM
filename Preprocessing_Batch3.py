import GEOparse
import pandas as pd

gse3 = GEOparse.get_GEO(filepath="./Data/Human/GSE27949_family.soft.gz")

plats_3 = list(gse3.gpls.keys())[0]

samples3 = gse3.phenotype_data["characteristics_ch1.1.bmi"]
samples3 = pd.DataFrame(samples3); samples3.head()
samples3.rename(columns={"characteristics_ch1.1.bmi":"bmi"}, inplace=True)
samples3["cbmi"] = samples3["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )
samples3 = samples3[['cbmi', 'bmi']]

samples3.to_pickle('./Preprocessed_Data/Human/batch3_pheno.p')
with open('./Preprocessed_Data/Human/batch3_pheno.txt', 'w') as handle:
    samples3.to_csv(handle, sep='\t')

samples3_ = gse3.pivot_samples('VALUE')[list(samples3.index)]

samples3_ann = samples3_.reset_index().merge(gse3.gpls['GPL570'].table[["ID", "Gene Symbol"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
del samples3_ann["ID"]
samples3_ann['Gene Symbol'] = samples3_ann['Gene Symbol'].astype(str)


samples3_ann = samples3_ann[~samples3_ann['Gene Symbol'].str.contains("///")]

samples3_ann = samples3_ann.dropna()

exprs_3 = samples3_ann.groupby('Gene Symbol').median()

del exprs_3.index.name
exprs_3 = exprs_3.T
exprs_3.to_pickle('./Preprocessed_Data/Human/batch3_geno.p')
with open('./Preprocessed_Data/Human/batch3_geno.txt', 'w') as handle:
    exprs_3.to_csv(handle, sep='\t')
