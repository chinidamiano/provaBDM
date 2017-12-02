import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy.stats import linregress,ks_2samp

geno1234 = pd.read_pickle("./output/1-Replication/batch1234_geno.p")
pheno1234 = pd.read_pickle("./output/1-Replication/batch1234_pheno.p")

geno5  = pd.read_pickle("./Preprocessed_Data/Human/batch5_geno.p")
pheno5 = pd.read_pickle("./Preprocessed_Data/Human/batch5_pheno.p")

geno6  = pd.read_pickle("./Preprocessed_Data/Human/batch6_geno.p")
pheno6 = pd.read_pickle("./Preprocessed_Data/Human/batch6_pheno.p")

geno7  = pd.read_pickle("./Preprocessed_Data/Human/batch7_geno.p")
pheno7 = pd.read_pickle("./Preprocessed_Data/Human/batch7_pheno.p")

geno8  = pd.read_pickle("./Preprocessed_Data/Human/batch8_geno.p")
pheno8 = pd.read_pickle("./Preprocessed_Data/Human/batch8_pheno.p")

def bmi_to_cbmi(x):
    assert type(x) is float or type(x) is np.float64
    if x<=25:
        return "lean"
    elif x<=30:
        return "overweight"
    else:
        return "obese"

pheno5["bmi"] = pheno5.bmi.astype(float)
pheno5["cbmi"] = pheno5.bmi.astype(float).apply(bmi_to_cbmi)
pheno6["bmi"] = pheno6["bmi (kg/m2)"].astype(float)
pheno6["cbmi"] = pheno6.bmi.apply(bmi_to_cbmi)
pheno7["bmi"] = pheno7.bmi.astype(float)
pheno7["cbmi"] = pheno7.bmi.apply(bmi_to_cbmi)
pheno6["BMI"] = pheno6.bmi
pheno6["FPI"] = pheno6["fasting plasma insulin (\xce\xbcu/ml)"].astype(float)
pheno6["FPG"] = pheno6["fasting plasma glucose (mg/dl)"].astype(float)

idx = pheno8.loc[pheno8.bmi==" Unk"].index
pheno8.drop(idx,axis=0,inplace=True)
geno8.drop(idx,axis=0,inplace=True)

pheno8["cbmi"]=0
for lab,df in pheno8.groupby("bmi"):
    if lab==" <25":
        pheno8.loc[df.index,"cbmi"]="lean"
    if lab==" 25-29.99":
        pheno8.loc[df.index,"cbmi"]="overweight"
    if lab==" 30+":
        pheno8.loc[df.index,"cbmi"]="obese"
pheno8["bmi"]=np.nan

def normalize_df(df):

    x = pow(2.,df)

    x = x.div(x.sum(axis=1), axis=0)

    return np.log2(x)

geno5 = normalize_df(geno5)
geno6 = normalize_df(geno6)
geno7 = normalize_df(geno7)
geno8 = normalize_df(geno8)

signature = pd.read_pickle("./output/1-Replication/signature.p")
signature.index = signature.index.get_level_values(0)

def compute_score(geno,signature):
    common_genes = np.intersect1d(geno.columns,signature.index)
    rawscore = geno[common_genes].dot(signature.loc[common_genes].coef)
    return rawscore - np.mean(rawscore)

pheno1234["score"] = compute_score(geno1234,signature)
pheno5["score"] = compute_score(geno5,signature)
pheno6["score"] = compute_score(geno6,signature)
pheno7["score"] = compute_score(geno7,signature)
pheno8["score"] = compute_score(geno8,signature)

### save batches for other notebooks
# geno5.to_pickle("./output/batches5-9_normalized/batch5_geno.p")
# pheno5.to_pickle("./output/batches5-9_normalized/batch5_pheno.p")
#
# geno6.to_pickle("./output/batches5-9_normalized/batch6_geno.p")
# pheno6.to_pickle("./output/batches5-9_normalized/batch6_pheno.p")
#
# geno7.to_pickle("./output/batches5-9_normalized/batch7_geno.p")
# pheno7.to_pickle("./output/batches5-9_normalized/batch7_pheno.p")
#
# geno8.to_pickle("./output/batches5-9_normalized/batch8_geno.p")
# pheno8.to_pickle("./output/batches5-9_normalized/batch8_pheno.p")


# ### Parameters for figures formatting
# This helps plotting all figures with the same parameters and color palettes:

labelsize=16  # fontsize of x,y labels
title_size=18  # fontsize of the title
ticks_size=14  # fontsize of numbers in x,y ticks
latex_size=18  # fontsize of latex-rendered pval, corr_coef annotations
xprop = 0.65   # relative x-position of latex-rendered pval, corr_coef annotations
yprop = 0.1    # relative y-position of latex-rendered pval, corr_coef annotations

colors_dict = {}

colors_dict["lean"] = sns.palettes.light_palette("green")[3]
colors_dict["overweight"] = sns.palettes.light_palette("orange")[3]
colors_dict["obese"] = sns.palettes.light_palette("red")[3]

colors_dict["NT"] = sns.palettes.light_palette("#f6ff42")[4]
colors_dict["TP"] = sns.palettes.light_palette("blue")[4]

plt.figure(figsize=(6,4))

sns.boxplot(
    data = pheno1234,x="cbmi",y="score",
    order = ["lean","obese"],
    palette=colors_dict,
    saturation=1,
    width=0.5*2/3
)

xlim = plt.xlim()
ylim = plt.ylim()
plt.ylabel("SCORE",size=labelsize)
plt.xlabel("")
plt.title("BATCHES 1 - 4",size=title_size)
plt.subplots_adjust(left=0.15)
plt.xticks([0,1],["Lean","Obese"])
plt.tick_params(labelsize=ticks_size)
plt.savefig("./output/1-Replication/figures/Figure1e.pdf")
#plt.show()

pheno12341 = pheno1234.loc[~(pheno1234.bmi.isnull())]
pheno = pheno12341

plt.figure(figsize=(6,4))
sns.regplot(pheno.bmi,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = plt.ylim()
sns.regplot(pheno.bmi.astype(float),pheno.score, fit_reg=True, scatter=False, x_bins=6, x_estimator=np.mean, color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("BMI",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 3",size=title_size)

plt.subplots_adjust(left=0.15,bottom=0.15)

plt.tick_params(labelsize=ticks_size)
if False:
    plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
             ylim[0]+yprop*(ylim[1]-ylim[0]),
             "$R=0.75$",
             size=latex_size,
             horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure2a.pdf")
#plt.show()

pheno = pheno5
print("p=",linregress(pheno.bmi,pheno.score).pvalue)
print("R=",linregress(pheno.bmi,pheno.score).rvalue)

plt.figure(figsize=(6,4))
sns.regplot(pheno.bmi,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = [-3,3]
sns.regplot(pheno.bmi,pheno.score, fit_reg=True,scatter=True,x_bins=6,x_estimator=np.mean, color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("BMI",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 5",size=title_size)

plt.subplots_adjust(bottom=0.15)

plt.tick_params(labelsize=ticks_size)
plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
         ylim[0]+yprop*(ylim[1]-ylim[0]),
         "$R=0.47$ \n$\hspace{0.3}p=7.3 \cdot 10^{-7}$",
         size=latex_size,
         horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure2b.pdf")
#plt.show()

pheno = pheno6
print(linregress(pheno.bmi,pheno.score).pvalue)
print(linregress(pheno.bmi,pheno.score).rvalue)

plt.figure(figsize=(6,4))
sns.regplot(pheno.bmi,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = plt.ylim()
sns.regplot(pheno.bmi,pheno.score,            fit_reg=True,scatter=True,x_bins=6,x_estimator=np.mean,            color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("BMI",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 6",size=title_size)

plt.subplots_adjust(bottom=0.15)

plt.tick_params(labelsize=ticks_size)
plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
         ylim[0]+yprop*(ylim[1]-ylim[0]),
         "$R=0.59$ \n$\hspace{0.3}p=2.7 \cdot 10^{-7}$",
         size=latex_size,
         horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure2c.pdf")
#plt.show()

pheno = pheno7
print("p=",linregress(pheno.bmi,pheno.score).pvalue)
print("R=",linregress(pheno.bmi,pheno.score).rvalue)

plt.figure(figsize=(6,4))
sns.regplot(pheno.bmi,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = [-1.7,1.5]
sns.regplot(pheno.bmi,pheno.score,            fit_reg=True,scatter=True,x_bins=6,x_estimator=np.mean,            color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("BMI",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 7",size=title_size)

plt.subplots_adjust(bottom=0.15)

plt.tick_params(labelsize=ticks_size)
plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
         ylim[0]+yprop*(ylim[1]-ylim[0]),
         "$R=0.27$ \n$\hspace{0.3}p=2.2 \cdot 10^{-2}$",
         size=latex_size,
         horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure2d.pdf")
#plt.show()

bmi_min = 25
bmi_max = 30
pheno66 = pheno6.loc[ pheno6.bmi < bmi_max ].loc[ pheno6.bmi > bmi_min ]

pheno = pheno66
print(linregress(pheno.FPG,pheno.score).pvalue)
print(linregress(pheno.FPG,pheno.score).rvalue)

plt.figure(figsize=(6,4))
sns.regplot(pheno.FPG,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = plt.ylim()
sns.regplot(pheno.FPG,pheno.score,            fit_reg=True,scatter=True,x_bins=4,x_estimator=np.mean,            color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("FPG",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 6\n$25 < \mathrm{BMI}<30$",size=title_size)

plt.subplots_adjust(top=0.85,left=0.15)

plt.tick_params(labelsize=ticks_size)
plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
         ylim[0]+yprop*(ylim[1]-ylim[0]),
         "$R=0.59$ \n$\hspace{0.3}p=2.7 \cdot 10^{-3}$",
         size=latex_size,
         horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure6d.pdf")
#plt.show()

pheno = pheno66
print(linregress(pheno.FPI,pheno.score).pvalue)
print(linregress(pheno.FPI,pheno.score).rvalue)

plt.figure(figsize=(6,4))
sns.regplot(pheno.FPI,pheno.score,fit_reg=False,color="black")
xlim = plt.xlim()
ylim = plt.ylim()
sns.regplot(pheno.FPI,pheno.score,            fit_reg=True,scatter=True,x_bins=4,x_estimator=np.mean,            color=sns.palettes.color_palette()[2])
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel("FPI",size=labelsize)
plt.ylabel("SCORE",size=labelsize)
plt.title("BATCH 6\n$25 < \mathrm{BMI}<30$",size=title_size)

plt.subplots_adjust(top=0.85,left=0.15)

plt.tick_params(labelsize=ticks_size)
plt.text(xlim[0]+xprop*(xlim[1]-xlim[0]),
         ylim[0]+yprop*(ylim[1]-ylim[0]),
         "$R=0.46$ \n$\hspace{0.3}p=2.9 \cdot 10^{-2}$",
         size=latex_size,
         horizontalalignment="left")

plt.savefig("./output/1-Replication/figures/Figure6b.pdf")
#plt.show()

plt.figure(figsize=(6,4))
plt.scatter(pheno6.BMI,pheno6.FPG,
            c=pheno6.score,
            cmap="Blues",
            edgecolors="black",
            alpha=0.8,
            s=30
           )
plt.ylim([72,120])
plt.xlabel("BMI",size=labelsize)
plt.ylabel("FPG",size=labelsize)
plt.title("BATCH 6",size=title_size)

plt.subplots_adjust(left=0.15)
plt.tick_params(labelsize=ticks_size)
cb = plt.colorbar(shrink=0.65)
cb.set_label("SCORE",fontsize=ticks_size)

plt.savefig("./output/1-Replication/figures/Figure6c.pdf")
#plt.show()

plt.figure(figsize=(6,4))
plt.scatter(pheno6.BMI,pheno6.FPI,
            c=(pheno6.score.values),
            cmap="Blues",
            edgecolors="black",
            alpha=0.8,
            s=30
           )
plt.ylim([-2,36])
plt.xlabel("BMI",size=labelsize)
plt.ylabel("FPI",size=labelsize)
plt.title("BATCH 6",size=title_size)

plt.subplots_adjust(left=0.15)
plt.tick_params(labelsize=ticks_size)
cb = plt.colorbar(shrink=0.65)
cb.set_label("SCORE",fontsize=ticks_size)

plt.savefig("./output/1-Replication/figures/Figure6a.pdf")
#plt.show()
