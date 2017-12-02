import GEOparse
import pandas as pd
import numpy as np
from functools import *
import re

gse1 = GEOparse.get_GEO(filepath="./Data/Human/GSE2508_family.soft.gz")

plats_1 = list(gse1.gpls.keys())

samples1 = gse1.phenotype_data[["platform_id", "title"]]
sample1 = samples1.groupby(["platform_id"]); sample1.groups
d = {}
for l in plats_1:
    ls = "".join(list(sample1.get_group(l)['title']))
    lf = re.findall("Lean F", ls)
    of = re.findall("Obese F", ls)
    lm = re.findall("Lean M", ls)
    om = re.findall("Obese M", ls)
    d[l] = {"LF": len(lf), "OF": len(of), "LM": len(lm), "OM": len(om)}

x = samples1.copy()
x["samples"] = x.index
x["title"] = x['title'].apply(lambda x: x[:-len(x.split()[-1])].strip()).to_frame('samples')
x['gender'] = x['title'].map(lambda x: x.split(' ')[1])
x['cbmi'] = x['title'].map(lambda x: x.split(' ')[0].lower())

grouped = x.groupby("title")
l = pd.DataFrame.from_dict(grouped.groups)

y = x[["title", "gender", "cbmi"]]

y = y.drop_duplicates("title")

y.index = y.title; y = y[["gender", "cbmi"]]

y.to_pickle('./Preprocessed_Data/Human/batch1_pheno.p')
with open('./Preprocessed_Data/Human/batch1_pheno.txt', 'w') as handle:
    y.to_csv(handle, sep='\t')

plats_1 = list(gse1.gpls.keys())

d = {}
samples1 = gse1.phenotype_data[["platform_id", "title"]]
sample1 = samples1.groupby(["platform_id"])

for plat in plats_1:
    d[plat] = gse1.pivot_samples('VALUE')[list(sample1.get_group(plat)[["title"]].index)]

d_ann = {}
for key in d.keys():
    d_ann[key] = d[key].reset_index().merge(gse1.gpls[key].table[["ID", "Gene Symbol"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
    d_ann[key].drop(['ID'], axis=1, inplace=True)


for key in d_ann.keys():
    d_ann[key].dropna(inplace=True)

for key in d_ann.keys():
    idx = d_ann[key]["Gene Symbol"].str.contains("///")
    idx = idx[idx==True].index
    d_ann[key] = d_ann[key].drop(list(idx), axis=0)
    d_ann[key] = d_ann[key].groupby("Gene Symbol").mean()

for key in d_ann.keys():
    d_ann[key] = np.log2(d_ann[key])

s = [d_ann[key] for key in d_ann.keys()]


df_final = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True, how='outer'), s)

df_c = {}

for column in l.columns:
    df_c[column] = df_final[list(l[column])]

for key in df_c.keys():
    df_c[key] = df_c[key].mean(axis=1)
    df_c[key] = df_c[key].to_frame()
    df_c[key] = df_c[key].rename(columns={0: key})

p = [df_c[key] for key in df_c.keys()]


df_f = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True, how='outer'), p)

del df_f.index.name

exprs_1 = df_f.T
exprs_1.to_pickle('./Preprocessed_Data/Human/batch1_geno.p')
with open('./Preprocessed_Data/Human/batch1_geno.txt', 'w') as handle:
    exprs_1.to_csv(handle, sep='\t')
