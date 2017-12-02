import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ks_2samp

import copy

def normalize_df(df):
    """
    input: pandas dataframe, values are log2 gene expression, genes in columns, samples in rows.
    output: same, with sum of genes equal 1 for all samples (in linear space).
    """
    x = pow(2.,df)

    x = x.div(x.sum(axis=1), axis=0)

    return np.log2(x)

gselist = ["batch1","batch2","batch3","batch4"]
dfs = {}
for gse in gselist:
    dfs[gse]=(normalize_df(pd.read_pickle("./Preprocessed_Data/Human/"+gse+"_geno.p")))
geno = pd.concat(dfs)
geno.index = geno.index.droplevel(0)

dfs = {}
for gse in gselist:
    tmp = pd.read_pickle("./Preprocessed_Data/Human/"+gse+"_pheno.p")
    tmp["batch"] = gse
    dfs[gse]=(tmp)
pheno = pd.concat(dfs)
pheno.index = pheno.index.droplevel(0)

geno.drop(geno.columns[geno.isnull().sum()!=0],axis=1,inplace=True)

print("batch \t\t lean \t overw \t obese \t total")
print("---")
for gse in gselist:
    l = ((pheno["batch"]==gse) & (pheno["cbmi"]=="lean")).sum()
    over = ((pheno["batch"]==gse) & (pheno["cbmi"]=="overweight")).sum()
    o = ((pheno["batch"]==gse) & (pheno["cbmi"]=="obese")).sum()
    print("%s \t %d \t %d \t %d \t %d" % (gse,l,over,o,l+over+o))
l =  (pheno["cbmi"]=="lean").sum()
over = ((pheno["cbmi"]=="overweight")).sum()
o = ((pheno["cbmi"]=="obese")).sum()
print("---")
print("%s \t %d \t %d \t %d \t %d" % ("ALL\t",l,over,o,l+over+o))

idx = pheno[(pheno["cbmi"]=="overweight")].index

geno.drop(idx,axis=0,inplace=True)
pheno.drop(idx,axis=0,inplace=True)

print("genes:",geno.shape[1])
print("samples (lean):",geno.loc[pheno.cbmi=="lean"].shape[0])
print("samples (obese):",geno.loc[pheno.cbmi=="obese"].shape[0])

ogeno = copy.deepcopy(geno)

def find_level(df,pheno,n=10):
    """
    returns 1./pvalue of distro of n first components being from different distros when split by cbmi
    args:
    df = dataframe with expression values
    pheno = dataframe with labels and metadata
    n = max numb of components
    """
    n=min(n,df.shape[0])

    res = []
    pca = PCA()
    trans = pca.fit_transform(df)
    df_trans = pd.DataFrame(data = trans[:,:n],index = df.index)

    for i in range(n):
        res.append([i,1/ks_2samp(
            df_trans.loc[(pheno.cbmi=="lean") | (pheno.cbmi=="overweight")].iloc[:,i],
            df_trans.loc[pheno.cbmi=="obese"].iloc[:,i]
        ).pvalue])
    return np.array(res).T


plt.figure(figsize=(17,3))
for i,batch in enumerate(gselist):
    plt.subplot(1,4,i+1)
    plt.title(batch,size=14)
    plt.xlabel("eigengene index",size=14)
    if i==0: plt.ylabel("$p$-value, KS lean/obese",size=14)
    y = np.log10(find_level(geno.loc[pheno.batch==batch],pheno)[1])
    plt.plot(range(1,y.shape[0]+1),y)
    yticks = plt.yticks()[0][::2]
    plt.yticks(yticks,["$10^{-%1.1f}$"%i for i in yticks],fontsize=14)
    plt.subplots_adjust(wspace=0.3)
plt.savefig("./output/1-Replication/figures/FigureP1.pdf")

pca = PCA()

for batch in gselist:
    effect_strength = find_level(geno.loc[pheno.batch==batch],pheno)[1]
    while np.argmax(effect_strength)!=0:
        tmp = pca.fit_transform(geno.loc[pheno.batch==batch])
        tmp[:,0]=0
        print("One principal component was set to zero on batch",batch)
        geno.loc[pheno.batch==batch] = pca.inverse_transform(tmp)
        effect_strength = find_level(geno.loc[pheno.batch==batch],pheno,n=min(10,geno.loc[pheno.batch==batch].shape[0]))[1]

plt.figure(figsize=(17,3))
for i,batch in enumerate(gselist):
    plt.subplot(1,4,i+1)
    plt.title(batch,size=14)
    plt.xlabel("eigengene index",size=14)
    if i==0: plt.ylabel("$p$-value, KS lean/obese",size=14)
    y = np.log10(find_level(geno.loc[pheno.batch==batch],pheno)[1])
    plt.plot(range(1,y.shape[0]+1),y)
    yticks = plt.yticks()[0][::2]
    plt.yticks(yticks,["$10^{-%1.1f}$"%i for i in yticks],fontsize=14)
    plt.subplots_adjust(wspace=0.3)
plt.savefig("./output/1-Replication/figures/FigureP2.pdf")

plt.title("ALL GSE's",size=14)
plt.plot(range(1,11),np.log10(find_level(geno,pheno)[1]))
plt.xlabel("eigegene index",size=14)
plt.ylabel("$p$-value (KS lean/obese)",size=14)
plt.xticks(fontsize=14)
yticks = plt.yticks()[0][::2]
plt.yticks(yticks,["$10^{-%d}$"%i for i in yticks],fontsize=16)
eigvect_k = np.argmax(find_level(geno,pheno)[1])
print("Picking eigengene number",eigvect_k+1)

pca = PCA()
trans = pca.fit_transform(geno)
trans[:,:eigvect_k]=0
n_geno = pd.DataFrame(index = geno.index,columns=geno.columns,data=pca.inverse_transform(trans))

K=7

pca_geno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(geno),
    index = geno.index).T.corr()

pca_ogeno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(ogeno),
    index = ogeno.index).T.corr()

pca_ngeno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(n_geno),
    index = ogeno.index).T.corr()

matrix_colors_dict = dict(zip(np.unique(pheno.batch),
                              np.array([sns.palettes.light_palette(x)[2] for x in sns.color_palette(n_colors=7)])[[3,4,5,6]]
                             ))
matrix_colors_dict["lean"] = sns.palettes.light_palette("green")[1]
matrix_colors_dict["obese"] = sns.palettes.light_palette("red")[1]

order = pheno.sort_values(by=["batch","cbmi"][::1]).index
batchcolors_bybatch = [matrix_colors_dict[pheno.loc[x,"batch"]] for x in order]
cbmicolors_bybatch = [matrix_colors_dict[pheno.loc[x,"cbmi"]] for x in order]

order = pheno.sort_values(by=["batch","cbmi"][::-1]).index
batchcolors_bycbmi = [matrix_colors_dict[pheno.loc[x,"batch"]] for x in order]
cbmicolors_bycbmi = [matrix_colors_dict[pheno.loc[x,"cbmi"]] for x in order]

batch_sizes = np.array([39,20,23,6])
cbmi_size = np.array([39,49])

order = pheno.sort_values(by=["batch","cbmi"][::1]).index
big_ax = sns.clustermap(pca_geno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bybatch,
               col_colors=batchcolors_bybatch,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]
nn=0
for b,n in enumerate(batch_sizes[:-1]):
    nn=nn+n
    ax.axhline(N-nn,color="black")
    ax.axvline(nn,color="black")
    ax.text(nn-n/2,N+1,"Batch"+str(b+1),rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
ax.text(N-3,N+1,"B4",rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
plt.savefig("./output/1-Replication/figures/Figure1a.pdf")

order = pheno.sort_values(by=["batch","cbmi"][::1]).index
big_ax = sns.clustermap(pca_ngeno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bybatch,
               col_colors=batchcolors_bybatch,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]
nn=0
for b,n in enumerate(batch_sizes[:-1]):
    nn=nn+n
    ax.axhline(N-nn,color="black")
    ax.axvline(nn,color="black")
    ax.text(nn-n/2,N+1,"Batch"+str(b+1),rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
ax.text(N-3,N+1,"B4",rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
plt.savefig("./output/1-Replication/figures/Figure1b.pdf")

order = pheno.sort_values(by=["batch","cbmi"][::-1]).index
big_ax = sns.clustermap(pca_ngeno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bycbmi,
               col_colors=batchcolors_bycbmi,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]

ax.axhline(N-39,color="black")
ax.axvline(39,color="black")
ax.text(-4.5,49/2,"OBESE",rotation=90,fontsize=13.5,color="black",verticalalignment="center")
ax.text(-4.5,N-39/2,"LEAN",rotation=90,fontsize=13.5,color="black",verticalalignment="center")
plt.savefig("./output/1-Replication/figures/Figure1c.pdf")

n_geno.to_pickle("./output/1-Replication/batch1234_geno.p")
pheno.to_pickle("./output/1-Replication/batch1234_pheno.p")
