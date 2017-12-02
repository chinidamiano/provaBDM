import GEOparse
import pandas as pd

gse6 = GEOparse.get_GEO(filepath="./Data/Human/GSE64567_family.soft.gz")

plats_6 = list(gse6.gpls.keys())[0]

# Annotation table

samples6 = gse6.phenotype_data[["characteristics_ch1.3.age","characteristics_ch1.4.bmi (kg/m2)", "characteristics_ch1.7.diastolic blood pressure (mm hg)", "characteristics_ch1.11.fasting plasma glucose (mg/dl)", "characteristics_ch1.12.fasting plasma insulin (μu/ml)", "characteristics_ch1.2.gender", "characteristics_ch1.9.high density lipoprotein cholesterol (mg/dl)", "characteristics_ch1.1.sample id", "characteristics_ch1.6.systolic blood pressure (mm hg)", "characteristics_ch1.0.tissue", "characteristics_ch1.8.total cholesterol (mg/dl)", "characteristics_ch1.10.triacylglycerol (mg/dl)", "characteristics_ch1.5.waist circumference (cm)"]]
samples6 = pd.DataFrame(samples6)
samples6.rename(columns={"characteristics_ch1.3.age": "age","characteristics_ch1.4.bmi (kg/m2)":"bmi (kg/m2)",\
"characteristics_ch1.7.diastolic blood pressure (mm hg)":"diastolic blood pressure (mm hg)",\
"characteristics_ch1.11.fasting plasma glucose (mg/dl)":"fasting plasma glucose (mg/dl)",\
"characteristics_ch1.12.fasting plasma insulin (μu/ml)":"fasting plasma insulin (Î¼u/ml)",\
"characteristics_ch1.2.gender":"gender", "characteristics_ch1.9.high density lipoprotein cholesterol (mg/dl)": "high density lipoprotein cholesterol (mg/dl)",\
"characteristics_ch1.1.sample id":"sample id", "characteristics_ch1.6.systolic blood pressure (mm hg)":"systolic blood pressure (mm hg)",\
"characteristics_ch1.0.tissue":"tissue", "characteristics_ch1.8.total cholesterol (mg/dl)" : "total cholesterol (mg/dl)",\
"characteristics_ch1.10.triacylglycerol (mg/dl)":"triacylglycerol (mg/dl)", "characteristics_ch1.5.waist circumference (cm)":"waist circumference (cm)"},\
inplace=True)

samples6["cbmi"] = samples6["bmi (kg/m2)"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )

samples6.to_pickle('./Preprocessed_Data/Human/batch6_pheno.p')
with open('./Preprocessed_Data/Human/batch6_pheno.txt', 'w') as handle:
    samples6.to_csv(handle, sep='\t')


samples6_ = gse6.pivot_samples('VALUE')[list(samples6.index)]

samples6_ann = samples6_.reset_index().merge(gse6.gpls['GPL10558'].table[["ID", "Symbol"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')

del samples6_ann["ID"]
samples6_ann['Symbol'] = samples6_ann['Symbol'].astype(str)

samples6_ann = samples6_ann[~samples6_ann['Symbol'].str.contains("///")]

samples6_ann = samples6_ann.dropna()

exprs_6 = samples6_ann.groupby('Symbol').mean()
del exprs_6.index.name

exprs_6 = exprs_6.T
exprs_6.to_pickle('./Preprocessed_Data/Human/batch6_geno.p')
with open('./Preprocessed_Data/Human/batch6_geno.txt', 'w') as handle:
    exprs_6.to_csv(handle, sep='\t')
