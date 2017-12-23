
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.feature_selection import RFE
import scipy, scipy.stats
import statsmodels.formula.api as sm
import scipy.optimize
import math

import itertools

import collections
import csv



# ## cannot use panda data frame due to uneven columns for each row
# df_profile = pd.read_csv('output_tf_parents_coefs_cellcycle2.csv', delimiter=',',header=None)#, names=range(0,1))
# profile = [d1.tolist() for i1,d1 in df_profile.iterrows()]



f1 = open('output_tf_parents_coefs_cellcycle.csv')
profile=[]
for line in f1: profile.append([float(x) for x in line.strip().split(',')])



df_tfs = pd.read_csv('1_sig_tf.csv', delimiter=',',header=None)
df_tfs['id'] =df_tfs[0].str.cat("-"+ df_tfs[1])
tfs = [d1.tolist() for i1,d1 in df_tfs.iterrows()]
tf_name =[x for x in df_tfs['id'].values]



df_genes = pd.read_csv('1_mean_expression_cos_all4018_cellcycle.csv', delimiter=',',header=None)
df_genes['id'] =df_genes[0].str.cat("-"+ df_genes[1])
genes = [d2.tolist() for i2,d2 in df_genes.iterrows()]
gene_name = [x for x in df_genes['id'].values]



## draw  total regulators for each gene

regnum = [x[1] for x in profile]
plt.xlim(1,len(genes))
plt.ylim(5,20)
plt.plot(regnum,'ro-')
plt.xlabel("Cell cycle genes")
plt.ylabel("Total number of non-zero coefficients in the model")


network ={}
for i in range(0,len(profile)):
    line = profile[i]
    g =gene_name[i]
    parents =[]

    counter = 2

    for j in range(0,int(line[1])):
        tfac= tf_name[int(line[counter])-1]
        tfcoef= line[counter+1]

        ptup = (tfac,tfcoef)

        counter += 2

        parents.append(ptup)
    network[g] = parents


G = nx.Graph()

for gene in network:
    if len(network[gene]) >1:
        for neighbor in network[gene]:
            n,w = neighbor
            G.add_edges_from([(gene,n)],weight=round(w,2))




for gene in network:
    #G.add_node(gene)
    if gene == 'Rv2154c-ftsW' or gene == 'Rv2150c-ftsZ' or gene == 'Rv2145c-wag31' or gene == 'Rv3918c-parA':
        for neighbor in network[gene]:
            n,w = neighbor
            G.add_edges_from([(gene,n)],weight=round(w,2))



color_map=[]
for node in G:
    if node in gene_name:
        color_map.append('blue')
    else:
        color_map.append('red')

pos =  nx.spring_layout(G,scale=2)#,weight='weight')
pos = nx.circular_layout(G,scale=5)
#labels = nx.get_edge_attributes(G,'weight')
nx.draw(G,pos,with_labels=False,font_size=12,alpha=0.5,width=2,node_color=color_map)
#nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_size=10)
plt.draw()
plt.legend()
plt.show()



########### greedy algorithm - find master tf
tf_count={}
for gene in network:
    parents = network[gene]
    for tfactorpair in parents:
        tfactor = tfactorpair[0]
        if tfactor in tf_count:
            #tf_count[tfactor] += 1
            tf_count[tfactor] += abs(tfactorpair[1])
        else:
            tf_count[tfactor] = abs(tfactorpair[1])
            #tf_count[tfactor] = 1

tf_rank =[ key for key, value in sorted(tf_count.iteritems(), key=lambda (k,v): (v,k),reverse=True)]



######################## calculate best cut off #######################


result_rsq =[]

for k in range(1,26):
    print k

    rsq_table = []

    for egene in network:

        #egene = 'Rv0017c-rodA'

        y1 = df_genes[df_genes['id']==egene].transpose()
        y1 = y1.iloc[2:17,0]
        y1 = np.array(y1, dtype=float)


        tf_group_index = [x for x, d in df_tfs.iterrows() if d['id'] in tf_rank[0:k]]
        xmat = df_tfs.iloc[tf_group_index, 2:17].transpose()
        xmat_filt = np.array(xmat, dtype=float)
        model = linear_model.LinearRegression()

        # Train the model using the training sets
        model.fit( xmat_filt,y1)


        rsq_table.append(model.score(xmat_filt,y1))


    result_rsq.append(sum(rsq_table)/float(len(network)))


plt.ylim([0.5,1.1])
plt.xlim([1,24])
plt.plot(result_rsq,'bo-')#,label='Rsquare')
plt.xlabel("Number of Transcription Factors in the model")
plt.ylabel(" Average Coefficient of Determination (R^2)")
plt.legend()
plt.axhline(0.90,0,25)


## we need top 14:
print tf_rank[0:7]


###########################################################




tf_pass = tf_rank[0:7]

for gene in network:
    parents = network[gene]
    parents = [(tfactor,val) for (tfactor,val) in parents if tfactor in tf_pass]
    network[gene] = parents



### new linear model with only 7 tfs


interaction =[]

tf_group_index = [x for x, d in df_tfs.iterrows() if d['id'] in tf_pass]
xmat = df_tfs.iloc[tf_group_index, 2:17].transpose()
xmat_filt = np.array(xmat, dtype=float)


for egene in network:

    #egene = 'Rv0017c-rodA'

    y1 = df_genes[df_genes['id']==egene].transpose()
    y1 = y1.iloc[2:17,0]
    y1 = np.array(y1, dtype=float)



    model = linear_model.LinearRegression()

    model.fit( xmat_filt,y1)

    interaction.append(model.coef_)


df_int = pd.DataFrame(interaction,columns=tf_pass)
df_int.head(2)

tf_reg_num =[]
for tf in tf_pass:
    pos = df_int.loc[:,tf][df_int.loc[:,tf]>=0.1].count()
    neg = df_int.loc[:,tf][df_int.loc[:,tf]<-0.1].count()
    print pos,neg
    tf_reg_num.append([tf,pos,neg])

for i in tf_reg_num:
    print i[0],',', i[1],',', i[2]





######draw graph for entire cell cycle group


G = nx.Graph()

for gene in network:
    if len(network[gene]) >1:
        for neighbor in network[gene]:
            n,w = neighbor
            G.add_edges_from([(gene,n)],weight=round(w,2))




for gene in network:
    #G.add_node(gene)
    if gene == 'Rv2154c-ftsW' or gene == 'Rv2150c-ftsZ' or gene == 'Rv2145c-wag31' or gene == 'Rv3918c-parA':
        for neighbor in network[gene]:
            n,w = neighbor
            G.add_edges_from([(gene,n)],weight=round(w,2))



color_map=[]
for node in G:
    if node in gene_name:
        color_map.append('blue')
    else:
        color_map.append('red')

pos =  nx.spring_layout(G,scale=2)#,weight='weight')
pos = nx.circular_layout(G,scale=5)
#labels = nx.get_edge_attributes(G,'weight')
nx.draw(G,pos,with_labels=True,font_size=15,alpha=0.5,width=2,node_color=color_map)
#nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_size=10)
plt.draw()
plt.legend()
plt.show()


for items in network:
    print items, len(network[items])




G = nx.Graph()


for gene in network:
    #G.add_node(gene)
    if gene == 'Rv2145c-wag31' :
        for neighbor in network[gene]:
            n,w = neighbor
            G.add_edges_from([(gene,n)],weight=round(w,2))



color_map=[]
for node in G:
    if node in gene_name:
        color_map.append('blue')
    else:
        color_map.append('red')

pos =  nx.spring_layout(G,scale=2)#,weight='weight')
pos = nx.circular_layout(G,scale=5)
labels = nx.get_edge_attributes(G,'weight')
nx.draw(G,pos,with_labels=True,font_size=15,alpha=0.5,width=2,node_color=color_map)
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_size=10)
plt.draw()
plt.legend()
plt.show()









############################## DRAW AS BIPARTITE #############################################

result = []
for gene in network:
    if gene == 'Rv2154c-ftsW' or gene == 'Rv2150c-ftsZ' or gene == 'Rv2145c-wag31' or gene == 'Rv3918c-parA':
        for tf in network[gene]:
            print gene, tf[0], round(tf[1],2)
            temp = [gene, tf[0], round(tf[1],2)]
            result.append(temp)

    # for tf in network[gene]:
    #     print gene, tf[0], round(tf[1], 2)
    #     temp = [gene, tf[0], round(tf[1], 2)]
    #     result.append(temp)

resultdf = pd.DataFrame(result,columns=['gene','tf','weight'])

df= resultdf
B = nx.Graph()
B.add_nodes_from(df['gene'], bipartite=0)
B.add_nodes_from(df['tf'], bipartite=1)
B.add_weighted_edges_from([(row['gene'], row['tf'], row['weight']) for idx, row in df.iterrows()], weight='weight')

# print(B.edges(data=True))

pos = {node:[0, i] for i,node in enumerate(df['tf'])}
pos.update({node:[1, i] for i,node in enumerate(df['gene'])})


color_map=[]
for node in B:
    if node in gene_name:
        color_map.append('blue')
    else:
        color_map.append('red')


for p in pos:  # raise text positions
    pos[p][1] += 1.0

labels = nx.get_edge_attributes(B,'weight')
nx.draw(B,pos,with_labels=True,font_size=25,alpha=0.5,width=2,node_color=color_map)
nx.draw_networkx_edge_labels(B,pos,edge_labels=labels,font_size=12, label_pos=0.467)

plt.show()






##########################################################

aic_table = []
rss_table = []
ll_table = []

G = len(network)
TP = len(df_genes.iloc[1,2:17])
n = G * TP
for k in range(1,16):

    print k

    rss = 0.0

    for egene in network:

        #egene = 'Rv0017c-rodA'

        y1 = df_genes[df_genes['id']==egene].transpose()
        y1 = y1.iloc[2:17,0]
        y1 = np.array(y1, dtype=float)

        ### try to permute all combination -- takes too long
        #tf_group_index = list(itertools.permutations(range(0, len(network[egene])), 4))


        ### to test best k tf using recursive feature selection algorithm from the entire network

        tf_group_index = [x for x,d in df_tfs.iterrows() if d['id'] in  tf_rank]
        xmat = df_tfs.iloc[tf_group_index,2:17].transpose()
        #xmat = np.array(xmat, dtype=float)
        reg = linear_model.LinearRegression()
        rfe = RFE(reg, k)
        rfe.fit(xmat,y1)
        xmat_filt = xmat.iloc[:,[x for x,d in enumerate(rfe.ranking_) if d ==1]]
        xmat_filt = np.array(xmat_filt, dtype=float)
        model = linear_model.LinearRegression()
        model.fit(xmat_filt,y1)

        yhat = model.predict(xmat_filt)
        yhat2 = xmat_filt.dot(model.coef_)


        ### to test k tf according to tf ranking

        # tf_group = tf_rank[0:k]
        # tf_group_index = [x for x, d in df_tfs.iterrows() if d['id'] in tf_group]
        # xmat = df_tfs.iloc[tf_group_index, 2:17].transpose()
        # xmat = np.array(xmat, dtype=float)
        # model = linear_model.LinearRegression()
        # model.fit(xmat,y1)
        # #yhat = model.predict(xmat_filt)
        # yhat2 = xmat.dot(model.coef_)


        rss += sum((yhat2-y1)**2)


    ll = -0.5*(n+n*np.log(2*np.pi)+( n * np.log(rss/n)))

    aic = (-2.0 * ll) + (2.0 * (G* float(k + 1)))

    aic_table.append(aic)

    rss_table.append(rss)

    ll_table.append(ll)


plt.xlim([1,15])
plt.plot(aic_table,'ro-',label='aic')
#plt.plot(result_aicc,'bo-',label='aicc')
plt.plot(rss_table,'go-',label='rss')
plt.plot(ll_table,'go-',label='ll')
plt.legend()




plt.ylim([10,35])
plt.xlim([1,12])
plt.plot(result_rsq,'bo-',label='rsq')
plt.legend()












###########################################################


#
# tf_dist = []
# for items in tf_count:
#     tf_dist.append(tf_count[items])
# plt.plot(tf_dist)
# plt.show()
#
#
# tf_dist_unq= np.unique(tf_dist)
# tf_dist_unq
# threshold = tf_dist_unq[len(tf_dist_unq)-10]
# tf_count = {i:tf_count[i] for i in tf_count if tf_count[i] >threshold}
# tf_pass = [i for i in tf_count]




#
# y = df_genes.iloc[:,2:]
# ymat = np.array(y)
# gene_name
#
# x = df_tfs.iloc[:,2:]
# xmat= np.array(x)
# tf_name
#
#
# ssr_table = []
#
# for num in range(1,len(tf_rank[0:20])):
#
#     x1 = []
#
#     tf_totest = tf_rank[0:num]
#
#     for tf in tf_totest:
#         tf_ind = tf_name.index(tf)
#         x1.append(xmat[tf_ind])
#
#     x1 = np.matrix(x1).transpose()
#
#     ssr =0.0
#
#     for egene in range(0,len(profile)):
#
#         y1 = ymat[egene]
#
#
#         reg = linear_model.LinearRegression()
#
#         reg.fit(x1,y1)
#         yhat = reg.predict(x1)
#
#         ssr += sum([pow((a-b),2) for a,b in zip(y1,yhat)])
#
#         print ssr
#
#     ssr_table.append(ssr)
#
#
# plt.plot(ssr_table,'ro')
#

################## likelihood



#   define a function to calculate the log likelihood
# def calcLogLikelihood(yhat, y1):
#     n = len(y1)
#     error = y1-yhat
#     sigma = np.std(error)
#     f = ((1.0/(2.0*math.pi*sigma*sigma))**(n/2)) * np.exp(-1*((np.dot(error.T,error))/(2*sigma*sigma)))
#     return np.log(f)
#
#
# result_aicc = []
# result_aic =[]
# result_bic =[]
# result_llike=[]
# n=15.0
# for k in range(1,14):
#     print k
#
#     llike_table =[]
#
#     for egene in network:
#
#
#         yrow = df_genes[df_genes['id']==egene].values
#         y1 = yrow[:,2:17].transpose()
#
#
#         #tf_group_index = list(itertools.permutations(range(0, len(network[egene])), 4))
#
#         tf_group = tf_rank[0:k]
#         tf_group_index = [x for x,d in df_tfs.iterrows() if d['id'] in  tf_group]
#
#
#         xmat = np.array(df_tfs.iloc[tf_group_index,2:17]).transpose()
#
#         reg = linear_model.LinearRegression()
#         reg.fit(xmat, y1)
#         yhat = reg.predict(xmat)
#
#         llike = calcLogLikelihood(yhat.flatten(), y1.flatten())
#
#         llike_table.append(llike)
#
#
#     bic = (np.log(n)*k) - 2. * sum(llike_table)
#
#     aic = (2. * k) - 2. * sum(llike_table)
#
#     aicc = aic + ( (2.*k*(k+1) ) / ( n-k - 1 ))
#
#     result_aicc.append(aicc)
#     result_bic.append(bic)
#     result_aic.append(aic)
#     result_llike.append(sum(llike_table))
#
#
# plt.xlim([1,12])
#
# plt.plot(result_aic,'ro-',label='aic')
# plt.plot(result_aicc,'bo-',label='aicc')
# plt.plot(result_bic,'go-',label='bic')
# plt.plot(result_llike,'yo-',label='logl')
# plt.legend()
#
#
#




#   define a function to calculate the log likelihood
# def calcLogLikelihood(yhat, y1):
#     n = len(y1)
#     error = y1-yhat
#     sigma = np.std(error)
#     f = ((1.0/(2.0*math.pi*sigma*sigma))**(n/2)) * np.exp(-1*((np.dot(error.T,error))/(2*sigma*sigma)))
#     return np.log(f)
#
# def calcRsq(yhat, y1):
#     return 1- (sum((yhat-y1)**2) / sum((y1-np.mean(y1))**2))

#

egene = 'Rv0017c-rodA'

y1 = df_genes[df_genes['id'] == egene].transpose()
y1 = y1.iloc[2:17, 0]
y1 = np.array(y1, dtype=float)


aic =[]
for k in range(3,9):

    tf_group = tf_rank[0:k]
    tf_group_index = [x for x, d in df_tfs.iterrows() if d['id'] in tf_group]
    xmat = df_tfs.iloc[tf_group_index, 2:17].transpose()

    reg = linear_model.LinearRegression()
    rfe = RFE(reg, k)
    rfe.fit(xmat, y1)
    xmat_filt = xmat.iloc[:, [x for x, d in enumerate(rfe.ranking_) if d == 1]]
    xmat_filt = np.array(xmat_filt, dtype=float)
    model = sm.OLS(y1, xmat_filt).fit()


    aic.append(model.aic)

#plt.xlim([3,7])
plt.plot(aic,'ro-',label='aic')


