import GEOparse
import pandas as pd

gse7 = GEOparse.get_GEO(filepath="./Data/Human/GSE33526_family.soft.gz")

plats_7 = list(gse7.gpls.keys())[0]

samples7 = gse7.phenotype_data[["characteristics_ch2.1.bmi", "characteristics_ch2.0.tissue"]]

samples7 = pd.DataFrame(samples7)

samples7.rename(columns={"characteristics_ch2.1.bmi":"bmi", "characteristics_ch2.0.tissue":"tissue"}, inplace=True)

samples7["cbmi"] = samples7["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )

samples7.to_pickle('./Preprocessed_Data/Human/batch7_pheno.p')
with open('./Preprocessed_Data/Human/batch7_pheno.txt', 'w') as handle:
    samples7.to_csv(handle, sep='\t')

samples7_ = gse7.pivot_samples('VALUE')[list(samples7.index)]

samples7_ann = samples7_.reset_index().merge(gse7.gpls[plats_7].table[["ID", "GENE_SYMBOL"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')

del samples7_ann["ID"]
samples7_ann['GENE_SYMBOL'] = samples7_ann['GENE_SYMBOL'].astype(str)

samples7_ann = samples7_ann[~samples7_ann['GENE_SYMBOL'].str.contains("///")]

samples7_ann = samples7_ann.dropna()

exprs_7 = samples7_ann.groupby('GENE_SYMBOL').mean()
del exprs_7.index.name

exprs_7 = exprs_7.T
exprs_7.to_pickle('./Preprocessed_Data/Human/batch7_geno.p')
with open('./Preprocessed_Data/Human/batch7_geno.txt', 'w') as handle:
    exprs_7.to_csv(handle, sep='\t')
