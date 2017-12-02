import GEOparse
import pandas as pd
import os

gse5 = GEOparse.get_GEO(filepath="./Data/Human/GSE62117_family.soft.gz")

plats_5 = list(gse5.gpls.keys())[0]

d = pd.DataFrame(gse5.phenotype_data["title"])
d["keep"] = gse5.phenotype_data["title"].map(lambda x: x.split("-")[1] =="CNTRL")
ix_s = d[d["keep"] == True]

samples5 = gse5.phenotype_data.loc[list(ix_s.index)]
samples5 = samples5[["characteristics_ch1.1.age (yrs)", "characteristics_ch1.5.bmi", "characteristics_ch1.2.ethnicity", "characteristics_ch1.3.gender", "characteristics_ch1.0.subject id", "characteristics_ch1.7.time point", "characteristics_ch1.8.tissue", "characteristics_ch1.6.treated with", "characteristics_ch1.4.weight (kg)"]]
samples5.rename(columns={"characteristics_ch1.1.age (yrs)":"age (yrs)", "characteristics_ch1.5.bmi":"bmi", "characteristics_ch1.2.ethnicity":"ethnicity", "characteristics_ch1.3.gender":"gender", "characteristics_ch1.0.subject id":"subject id", "characteristics_ch1.7.time point":"time point", "characteristics_ch1.8.tissue":"tissue", "characteristics_ch1.6.treated with":"treated with", "characteristics_ch1.4.weight (kg)":"weight (kg)"}, inplace=True)
samples5["cbmi"] = samples5["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )

list(samples5["cbmi"]).count("overweight")
list(samples5["cbmi"]).count("lean")
list(samples5["cbmi"]).count("obese")

samples5 = samples5[~(samples5["cbmi"] == "overweight")]

samples5.to_pickle('./Preprocessed_Data/Human/batch5_pheno.p')
with open('./Preprocessed_Data/Human/batch5_pheno.txt', 'w') as handle:
    samples5.to_csv(handle, sep='\t')

samples5_ = gse5.pivot_samples('VALUE')[list(samples5.index)]


samples5_ann = samples5_.reset_index().merge(gse5.gpls[plats_5].table[["GENE", "GENE_SYMBOL"]],
                                left_on='ID_REF', right_index=True).set_index('ID_REF').dropna()

samples5_ann = samples5_ann[~samples5_ann['GENE_SYMBOL'].str.contains("///")]
exprs_5 = samples5_ann.groupby('GENE_SYMBOL').mean()
exprs_5.drop("GENE", axis=1, inplace=True)
del exprs_5.index.name

exprs_5 = exprs_5.dropna()

exprs_5.T.to_pickle('./Preprocessed_Data/Human/batch5_geno.p')
with open('./Preprocessed_Data/Human/batch5_geno.txt', 'w') as handle:
    exprs_5.T.to_csv(handle, sep='\t')
