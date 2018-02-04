import pandas as pd
import numpy as np
import statsmodels.api as sm
import  matplotlib.pylab as plt
from sklearn.cluster import KMeans,AgglomerativeClustering

# df_tfs = pd.read_csv('tfs_peak7_.csv',header=0)
#
# df_genes = pd.read_csv('genes_peak7_.csv',header=0)
#
# xdata = df_tfs.loc[:,'3hrs':'55hrs']
# #xdata = df_tfs.iloc[:,3:18]
# xdata = xdata.T
# xdata.columns = df_tfs['ORF ID']
# print xdata.head(1)
# plt.plot( xdata['Rv0038'])
# print xdata.shape
#
# df_genes.head(1)
# ydata = df_genes.loc[df_genes.loc[df_genes['gene'] == 'ftsW'].index.values.astype(int)[0],'3hrs':'55hrs']
# ydata = ydata.T
# print ydata.shape
#
# plt.plot( xdata['Rv0038'],'ro-')
# plt.plot(ydata,'bo-')


############# run one polynomial for rv0038 ###########
#
#
# x = xdata['Rv0038']
# y = np.array(ydata,dtype=float)
# x2 = x ** 2.0
# df = pd.concat([x, x2], axis=1)
# df.columns = ['x','x2']
# poly_2 = sm.OLS(y,df.values).fit()
#
# plt.plot( xdata['Rv0038'], 'ro-',label='original x')
# plt.plot(y,'bo-',label ='original y ')
# plt.plot( poly_2.predict(df.values), 'go-',label = 'poly2' + '/rsq:' + str(poly_2.rsquared_adj))
#
# x = xdata['Rv0038']
# y = np.array(ydata,dtype=float)
# x2 = x ** 2.0
# x3 = x ** 3.0
# df = pd.concat([x, x2,x3], axis=1)
# df.columns = ['x','x2','x3']
# poly_3 = sm.OLS(y,df.values).fit()
# plt.plot( poly_3.predict(df.values), 'mo-',label = 'poly3' + '/rsq:'  + str(poly_3.rsquared_adj))
#
# x = xdata['Rv0038']
# y = np.array(ydata,dtype=float)
# x2 = x ** 2.0
# x3 = x ** 3.0
# x4 = x ** 4.0
# df = pd.concat([x, x2,x3,x4], axis=1)
# df.columns = ['x','x2','x3','x4']
# poly_4 = sm.OLS(y,df.values).fit()
# plt.plot( poly_4.predict(df.values), 'yo-',label = 'poly4'+ '/rsq:' + str(poly_4.rsquared_adj))
#
#
#
#
# x = xdata['Rv0038']
# y = np.array(ydata,dtype=float)
# x2 = x ** 2.0
# x3 = x ** 3.0
# x4 = x ** 4.0
# x5 = x ** 5.0
# df = pd.concat([x, x2,x3,x4,x5], axis=1)
# df.columns = ['x','x2','x3','x4','x5']
# poly_5 = sm.OLS(y,df.values).fit()
# plt.plot( poly_5.predict(df.values), 'co-',label = 'poly5'+ '/rsq:'  + str(poly_5.rsquared_adj))
# plt.legend()

###############################################################################


### cluster tfs



df_tfs = pd.read_csv('tfs_peak7_.csv',header=0)

df_genes = pd.read_csv('genes_peak7_.csv',header=0)

xdata = df_tfs.loc[:,'3hrs':'55hrs']

xdata.index = df_tfs.loc[:,'ORF ID']



################ check number of cluster using kMeans error

# cluster_range = range( 1, 20 )
# cluster_errors = []
#
# for num_clusters in cluster_range:
#   clusters = KMeans( num_clusters )
#   clusters.fit( xdata )
#   cluster_errors.append( clusters.inertia_ )
# clusters_df = pd.DataFrame( { "num_clusters":cluster_range, "cluster_errors": cluster_errors } )
# clusters_df[0:10]
#
# plt.figure(figsize=(12,6))
# plt.plot( clusters_df.num_clusters, clusters_df.cluster_errors, marker = "o" )

##############################################################

clust_num = 6
#kmeans = KMeans(n_clusters=clust_num, random_state=0).fit(xdata)
hc = AgglomerativeClustering(n_clusters=clust_num).fit(xdata)
xdata['clust'] = hc.labels_
xdata.head(2)


######## plot clusters

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')

c=0
for ax in axs.flatten():
        xdata_1 = xdata.loc[xdata['clust'] == c]
        for x in range(0, len(xdata_1)):
            ax.plot(xdata_1.iloc[x, 0:15],'k-')
            ax.set_title(c+1)
        ax.plot(xdata_1.iloc[:, 0:15].mean(), 'r-')

        c += 1



print xdata.loc[xdata['clust']==4,:]


########################################



def rbf(m,v,x):
    sum =0.0
    for n in x:
        sum += np.exp(- ((m-n) ** 2)/(2*v))
    return np.exp(- ((m-sum/len(x)) ** 2)/(2*v))

xhat = []
for tps in range(0,15):
    t1_mean = xdata.groupby('clust').mean().iloc[:,tps]

    t1_var = xdata.groupby('clust').var().iloc[:,tps]

    x = xdata.iloc[:,tps]

    temp = []
    for i in range(0,6):
        temp.append(rbf(t1_mean[i],t1_var[i],x))

    xhat.append(temp)

df_xhat = pd.DataFrame(xhat)
print df_xhat
df_xhat.insert(6, 6, [1]*15)
df_xhat.columns = [ 1, 2, 3, 4, 5, 6,0]
df_xhat = df_xhat[[0, 1, 2, 3, 4, 5, 6]]
print df_xhat


plt.plot(df_xhat.iloc[:,1])
plt.plot(df_xhat.iloc[:,2])
plt.plot(df_xhat.iloc[:,3])
plt.plot(df_xhat.iloc[:,4])
plt.plot(df_xhat.iloc[:,5])
plt.plot(df_xhat.iloc[:,6])




############# check for ftsW
df_genes.head(1)
ydata = df_genes.loc[df_genes.loc[df_genes['gene'] == 'ftsW'].index.values.astype(int)[0],'3hrs':'55hrs']
ydata = ydata.T
print ydata
plt.plot(ydata)

model = sm.OLS(np.array(ydata,dtype=float),df_xhat.values).fit()
model.summary()
plt.plot(ydata)
plt.plot(model.predict(df_xhat.values))


####################### run for all ###################################



df_genes.head(1)
ydata = df_genes.loc[:,'3hrs':'55hrs']
ydata = ydata.T

yhat =[]
ycoef = []
for ind in range(0,536):
    model = sm.OLS(ydata.iloc[:,ind].values,df_xhat.values).fit()
    model.summary()
    yhat.append(model.predict(df_xhat.values))
    temp =[]
    for x,y in zip(model.params,model.pvalues):
        if y<0.05:
            temp.append(x)
        else:
            temp.append(0.0)
    ycoef.append(temp)



yhat2 = pd.DataFrame(yhat)
yhat2 = yhat2.T
print yhat2.shape

plt.plot(ydata.iloc[:,5],'ro-')
plt.plot(yhat2.iloc[:,5],'bo-')

ycoef2 = pd.DataFrame(ycoef)
ycoef2 = ycoef2.T
print ycoef2.shape
ycoef2.columns = df_genes.iloc[:,0] +'/' + df_genes.iloc[:,1]
ycoef2.to_csv("rbf_result.csv")


################################################################




