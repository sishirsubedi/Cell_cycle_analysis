import pandas as pd
import matplotlib.pylab as plt
import numpy as np


############## this section is before anova analysis

cos = pd.read_csv("2_deseqnorm_cos.csv")
cos = pd.read_csv("2_deseqnorm_rv.csv")
cos.head(2)

profile =[]

for ind,row in cos.iterrows():
    gene = row['rv']
    cos1 = row[1:16]
    cos2 = row[16:]


    for i in range(0,len(cos1)):
        temp =[gene]
        temp.extend(['cos1',str('tp_'+str(i+1)),cos1[i]])
        profile.append(temp)

    for i in range(0,len(cos2)):
        temp =[gene]
        temp.extend(['cos2',str('tp_'+str(i+1)),cos2[i]])
        profile.append(temp)

df_profile = pd.DataFrame(profile,columns=['rv','strain','tp','expr'])
df_profile.head(5)
df_profile.to_csv('2_deseqnorm_rv_foranova.csv',index=False)



################## this section is after anova analysis


cos = pd.read_csv("2_deseqnorm_cos.csv")
cos.head(2)

compare = pd.read_csv("compare_sig.csv")
compare.head(2)


protlen = pd.read_csv("H37Rv_protlen.csv")
print protlen.shape
protlen.head(4)



anova_sig_genes = [x for x in cos.rv.values if x not in compare.anova.values]
cos_all_anova_sig = cos[cos.rv.isin(anova_sig_genes)]
print cos_all_anova_sig.shape

cos_all_anova_sig.head(4)
protlen.head(4)


protlenhash = {}
for ind,row in protlen.iterrows():
    protlenhash[row[0]] = row[1]


plen =[]
for ind2, row2 in cos_all_anova_sig.iterrows():
    plen.append(protlenhash[row2[0]])


cos_all_anova_sig['plen'] = plen
cos_all_anova_sig_allvals = cos_all_anova_sig.iloc[:,1:]
cos_all_anova_sig['row_mean'] = cos_all_anova_sig_allvals.apply(lambda x : np.mean(x),axis=1)
cos_all_anova_sig.head(5)
cos_all_anova_sig['norm'] = cos_all_anova_sig.row_mean/cos_all_anova_sig.plen
cos_all_anova_sig.head(5)
plt.hist(cos_all_anova_sig.norm.values[cos_all_anova_sig.norm.values<2],bins=40)
plt.xlabel("Coverage estimate of significant genes", fontsize=16)
plt.ylabel("Number of genes", fontsize=16)

anova_nonsig_genes = [x for x in cos.rv.values if x in compare.anova.values]
cos_anova_nonsig = cos[cos.rv.isin(anova_nonsig_genes)]
print cos_anova_nonsig.shape

cos_anova_nonsig.head(4)
protlen.head(4)

plen = []
for ind2, row2 in cos_anova_nonsig.iterrows():
    plen.append(protlenhash[row2[0]])

cos_anova_nonsig['plen'] = plen
cos_anova_nonsig_allvals = cos_anova_nonsig.iloc[:, 1:]
cos_anova_nonsig['row_mean'] = cos_anova_nonsig_allvals.apply(lambda x: np.mean(x), axis=1)
cos_anova_nonsig.head(5)
cos_anova_nonsig['norm'] = cos_anova_nonsig.row_mean / cos_anova_nonsig.plen
cos_anova_nonsig.head(5)
plt.hist(cos_anova_nonsig.norm.values[cos_anova_nonsig.norm.values < 2], bins=40)
plt.xlabel("Coverage estimate of non significant genes", fontsize=16)
plt.ylabel("Number of genes", fontsize=16)





############## this section is before anova analysis

cos = pd.read_csv("2_deseqnorm_rv.csv")
cos.head(2)

profile =[]

for ind,row in cos.iterrows():
    gene = row['rv']
    cos1 = row[1:16]
    cos2 = row[16:]


    for i in range(0,len(cos1)):
        temp =[gene]
        temp.extend(['rv1',str('tp_'+str(i+1)),cos1[i]])
        profile.append(temp)

    for i in range(0,len(cos2)):
        temp =[gene]
        temp.extend(['rv2',str('tp_'+str(i+1)),cos2[i]])
        profile.append(temp)

df_profile = pd.DataFrame(profile,columns=['rv','strain','tp','expr'])
df_profile.head(5)
df_profile.to_csv('2_deseqnorm_rv_foranova.csv',index=False)



