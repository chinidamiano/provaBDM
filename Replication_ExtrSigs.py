import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp
import copy
import sklearn

geno = pd.read_pickle("./output/1-Replication/batch1234_geno.p")
pheno = pd.read_pickle("./output/1-Replication/batch1234_pheno.p")

v = (PCA().fit(geno)).components_[0]
v_df = pd.DataFrame(columns=["coef","abs_coef"],index = geno.columns)
v_df["coef"] = v
v_df["abs_coef"] = np.abs(v)
v_df.sort_values(by="abs_coef",ascending=False,inplace=True)
v_df.drop("abs_coef",axis=1,inplace=True)

v_df.head()

N = geno.shape[1]
sigma = 1/np.sqrt(N)

idx = np.abs(v_df.coef)>5*sigma
print("genes beyond threshold:",idx.sum())

signature = v_df.loc[idx]

v_df["abscoef"] = v_df["coef"].apply(np.abs)

tmp = np.array([sorted(np.abs(np.random.normal(size=v_df.shape[0])/np.sqrt(v_df.shape[0]))) for x in range(100)])
m=tmp.mean(axis=0)
s=tmp.std(axis=0)

labelsize=16  # fontsize of x,y labels
title_size=18  # fontsize of the title
ticks_size=14  # fontsize of numbers in x,y ticks
latex_size=18  # fontsize of latex-rendered pval, corr_coef annotations

plt.figure(figsize=(6,4))

t = range(1,N+1)
Ng = 2000
mycolor = "#03589E"
M=20
K2=5

plt.plot((t[:58]),v_df["abscoef"].values[:58],color="black",linestyle="-",linewidth="2",label="SCORE")
plt.plot((t[58:]),v_df["abscoef"].values[58:],color="black",alpha=0.35,linestyle="-",linewidth="2")

for j in [-1,1]:
    for K in range(0,M,5):
        plt.fill_between(range(Ng),(m+j*(K)*s)[::-1][:Ng],(m+j*(K+5)*s)[::-1][:Ng],
                         color=mycolor,linewidth=0,alpha=0.7*(M-K)/M)

plt.plot((m)[::-1][:Ng],color=mycolor,linestyle="-",linewidth=2,label="RANDOM SCORE")

plt.axhline(K2/np.sqrt(v_df.shape[0]),linewidth=1,color=mycolor,linestyle="--",label="%d$\sigma$ threshold" %K2,)

plt.ylim([0,.14])
plt.xlim([0,200])
plt.xlabel("Rank-ordered genes",size=labelsize)
plt.ylabel("Absolute value of coefficient",size=labelsize)
plt.xticks([1,38,50,100,150,200],size=ticks_size)
plt.yticks(size=ticks_size)
plt.legend(fontsize=14)

plt.savefig("./output/1-Replication/figures/Figure1d.pdf")

signature.to_csv("./output/1-Replication/signature.csv")
signature.to_pickle("./output/1-Replication/signature.p")

np.savetxt("./output/1-Replication/background_genes.txt",v_df.index.get_level_values(0).values,fmt="%s")
