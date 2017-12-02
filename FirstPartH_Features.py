import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import ks_2samp
import statsmodels
from statsmodels import *
import statsmodels.stats.multitest

from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV


import sklearn.tree as tr
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from scipy.stats import rankdata as rd
import operator

###############################################
#Kolmogorov-Smirnov (KS) + Bonferroni correction
###############################################

geno = pd.read_pickle("./output/2-First_part/batch123456_geno.p")
pheno = pd.read_pickle("./output/2-First_part/batch123456_pheno.p")

#print(geno.head(), pheno.head())

#print(geno.loc[pheno.cbmi=="lean"].T)

df_1 = geno.loc[pheno.cbmi=="lean"].T
df_2 = geno.loc[pheno.cbmi=="obese"].T

g_ks = geno.T.apply(lambda x: scipy.stats.ks_2samp(x[pheno.cbmi=="lean"],x[pheno.cbmi=="obese"]), axis=1)
g_ks = pd.DataFrame(g_ks)
g_ks = g_ks.rename(columns={0:'resul_ks'})
g_ks['p_value'] = g_ks.resul_ks.apply(lambda x: np.float(x[1]))
corrected_pvalues = statsmodels.stats.multitest.multipletests(pvals=g_ks['p_value'], method='bonferroni', alpha=0.00001)
g_ks['corrected_pvalue'] = corrected_pvalues[1]
g_ks['rejected'] = corrected_pvalues[0] #true for hypothesis that can be rejected for given alpha
g_ks = g_ks[g_ks.rejected==True]
#print(g_ks, g_ks.shape, type(g_ks))

ks_sigs = list(g_ks.index)
#print(ks_sigs)


###############################################
#LASSO
###############################################

y = pheno["cbmi"]

def label_to_int(x):
    if x == "lean":
        return 0
    else:
        return 1
y = y.apply(label_to_int).values
x = geno.as_matrix()

#Lasso feature selection
clf = LassoCV()

# Set a minimum threshold of 0.25
sfm = SelectFromModel(clf, threshold=0.05)
sfm.fit(x, y)
relevant_features_indices = sfm.get_support(indices=True)

r=sfm.transform(x)
n_features = sfm.transform(x).shape[1]

genes = list(geno.columns)
lasso_sigs = []

for index in relevant_features_indices:
    lasso_sigs.append(genes[index])


#print(lasso_sigs)

###############################################
#Random Forest (rf)
###############################################

y = pheno["cbmi"]
def label_to_int(x):
    if x == "lean":
        return 0
    else:
        return 1
y = y.apply(label_to_int).values
x = geno.as_matrix()

# Build a forest and compute the feature importances
runs = []
genes = list(geno.columns)
#ranks1 = {gene:[] for gene in genes}
ranks1 = {gene:0 for gene in genes}
for i in range(50,60):
    forest = ExtraTreesClassifier(n_estimators=100,
                                  random_state=i)
    
    forest.fit(x, y)

    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
     
    #t=rd(importances, method='average')
    #for l in range(0,len(t)):
    #    ranks1[genes[l]].append( t[l])
    
    for l in range(0,len(importances)):
        ranks1[genes[l]]+=importances[l]


sorted_genes1 = sorted(ranks1.items(), key=operator.itemgetter(1))[::-1]
signature = pd.DataFrame(sorted_genes1[:10])
rf_sigs = list(signature[0])
#print(rf_sigs)

###############################################
#Intersection
###############################################

intersection = set.union(set.intersection(set(ks_sigs), set(lasso_sigs)), set.intersection(set(ks_sigs), set(rf_sigs)), set.intersection(set(lasso_sigs), set(rf_sigs)))

their_sigs = pd.read_pickle("./output/2-First_part/theirH_signature.p")
their_sigs = set(their_sigs.index)

inter_our_their = set.intersection(intersection, their_sigs)

print(intersection)

with open("./output/2-First_part/ourH_signature.csv", "w") as handle:
    handle.write("\n".join(list(intersection)))
with open("./output/2-First_part/ourItheirH_signature.csv", "w") as handle:
    handle.write("\n".join(list(inter_our_their)))

