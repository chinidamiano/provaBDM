import GEOparse
import pandas as pd

gse8 = GEOparse.get_GEO(filepath="./Data/Human/GSE78958_family.soft.gz")

plats_8 = list(gse8.gpls.keys())[0]

samples8 = gse8.phenotype_data[["characteristics_ch1.1.bmi", "characteristics_ch1.0.patient ethnicity", "characteristics_ch1.2.tumor grade", "characteristics_ch1.4.tumor stage", "characteristics_ch1.3.tumor subtype (via breastprs)"]]

samples8 = pd.DataFrame(samples8)

samples8.rename(columns={"characteristics_ch1.1.bmi":"bmi", "characteristics_ch1.0.patient ethnicity":"patient ethnicity", "characteristics_ch1.2.tumor grade":"tumor grade", "characteristics_ch1.4.tumor stage":"tumor stage", "characteristics_ch1.3.tumor subtype (via breastprs)":"tumor subtype (via breastprs)"}, inplace=True)

samples8.to_pickle('./Preprocessed_Data/Human/batch8_pheno.p')
with open('./Preprocessed_Data/Human/batch8_pheno.txt', 'w') as handle:
    samples8.to_csv(handle, sep='\t')

samples8_ = gse8.pivot_samples('VALUE')[list(samples8.index)]

samples8_ann = samples8_.reset_index().merge(gse8.gpls[plats_8].table[["ID", "Gene Symbol"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')

del samples8_ann["ID"]
samples8_ann['Gene Symbol'] = samples8_ann['Gene Symbol'].astype(str)

samples8_ann = samples8_ann[~samples8_ann['Gene Symbol'].str.contains("///")]

samples8_ann = samples8_ann.dropna()

exprs_8 = samples8_ann.groupby('Gene Symbol').mean()
del exprs_8.index.name

exprs_8 = exprs_8.T
exprs_8.to_pickle('./Preprocessed_Data/Human/batch8_geno.p')
with open('./Preprocessed_Data/Human/batch8_geno.txt', 'w') as handle:
    exprs_8.to_csv(handle, sep='\t')
