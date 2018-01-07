# This algorithm is modified version of
# Cheng, Y., & Church, G. M.(2000, August). Biclustering of expression data. In Ismb (Vol. 8, pp. 93-103).
# here, biclustering is done only allowing rows(genes) to cluster together in fixed time-points(columns)
# good explanation of the paper
#http://www.kemaleren.com/cheng-and-church.html


import numpy as np
import pandas as pd

#df_cos= pd.read_csv('2_test.csv', delimiter=',', header=0)
df_cos= pd.read_csv('1_meanexpression_cos.csv', delimiter=',', header=0)
cos1 = []
for row in df_cos.iterrows():
    index, data = row
    cos1.append(data.tolist())

cos1m = df_cos.as_matrix()
cos1mdata = cos1m[:,4:19]
print cos1mdata.shape
#cos1mdata = df_cos.iloc[:,4:19]


# for i in range(0,10):
#     print cos1mdata[i,:]

def colwise_mean(cosm, rown, coln):
    #total_col = cosm.shape[1]-coln
    #print total_col
    sum =0
    for j in range(0,coln+1):
        sum += cosm[rown,j]
    return sum/(coln+1)

def rowwise_mean(cosm, rown, coln):
    sum =0
    for i in range(0,rown+1):
        sum += cosm[i,coln]
    return sum/(rown+1)


def biclust_mean(cosm, rown, coln):
    sum =0.0
    for i in range(0,rown+1):
        for j in range(0, coln +1):
            sum += cosm[i,j]
    return sum/float((rown+1)*(coln+1))


def variance(cosm, rown, coln):
    sum =0.0
    edge = biclust_mean(cosm, rown, coln)
    for i in range(0,rown+1):
        for j in range(0, coln +1):
            sum +=  (cosm[i,j] - edge)**2
    return sum/float((rown+1)*(coln+1))


def residual(cosm, i, j,I,J):
    rowmean = rowwise_mean(cosm, i,J)
    colmean = colwise_mean(cosm,I,j)
    rowcolmean = biclust_mean(cosm,I,J)
    temp =  pow((cosm[i,j]-rowmean-colmean+rowcolmean),2)
    return temp


def abs_residual(cosm, i, j,I,J):
    rowmean = rowwise_mean(cosm, i,J)
    colmean = colwise_mean(cosm,I,j)
    rowcolmean = biclust_mean(cosm,I,J)
    temp =  abs(cosm[i,j]-rowmean-colmean+rowcolmean)
    return temp

def remove_rows(cosm,I,J):
    temp = []
    val = 0.0
    for i in range(0,cosm.shape[0]):
        for j in range(0, cosm.shape[1]):
            val += residual(cosm,i,j,I,J)
            #val += abs_residual(cosm, i, j, I, J)
        val = val/float(cosm.shape[1]+1)
        temp.append(val)
    maxindex = temp.index(max(temp))
    return maxindex

def remove_rows_wth(cosm,I,J,th):
    temp = []
    delind =[]
    rval = 0.0
    for i in range(0,cosm.shape[0]):
        for j in range(0, cosm.shape[1]):
            rval += residual(cosm,i,j,I,J) # travel across ith row
            #val += abs_residual(cosm, i, j, I, J)
        rval = rval/float(cosm.shape[1]+1)
        temp.append(rval)
    for di in range(0,len(temp)):
        if temp[di] >= th:
            delind.append(di)
    return delind

def update_matrix(mat,index,rc):
    mat = np.array(mat)
    if rc =='r':
        mat = np.delete(mat,(index),axis=0)
    else:
        mat = np.delete(mat, (index), axis=1)
    mat = np.matrix(mat)
    return mat

def check_columns(mat):
    lowest = 100
    temp = 100
    for x in range(4, mat.shape[1]):
        mse =0.0
        for i in range(0, mat.shape[0]):
            for j in range(0, x):
                mse += residual(mat, i, j, mat.shape[0] - 1, x)
                #mse += abs_residual(mat, i, j, mat.shape[0] - 1, x)
        mse = mse/float((mat.shape[0])*(x))
        print x, mse
        if mse<=temp:
            temp = mse
            lowest = x
    return lowest

#index = (0,0)

delta = 0.5
alpha = 1.0
mse =0.0
for i in range(0,cos1mdata.shape[0]):
    for j in range(0, cos1mdata.shape[1]):
        mse += residual(cos1mdata, i, j,cos1mdata.shape[0]-1,cos1mdata.shape[1]-1)
        #mse += abs_residual(cos1mdata, i, j, cos1mdata.shape[0] - 1, cos1mdata.shape[1] - 1)
mse = mse/float((cos1mdata.shape[0])*(cos1mdata.shape[1]))
print "first mse", mse

while  mse > delta :

    # #remove single row one at a time
    # ind = remove_rows(cos1mdata, cos1mdata.shape[0] - 1, cos1mdata.shape[1] - 1)
    # print 'ind',ind
    # cos1mdata = np.array(cos1mdata)
    # cos1mdata = np.delete(cos1mdata, ind, axis=0)
    # cos1mdata = np.matrix(cos1mdata)

    oldmse = mse
    alpha -= 0.01
    inds = remove_rows_wth(cos1mdata,cos1mdata.shape[0]-1,cos1mdata.shape[1]-1,alpha*mse)

    if len(inds) >=1 :
        print 'group del inds', inds
        cos1mdata = np.array(cos1mdata)
        cos1mdata = np.delete(cos1mdata,inds,axis=0)
        cos1mdata = np.matrix(cos1mdata)

        for i in range(0, cos1mdata.shape[0]):
            for j in range(0, cos1mdata.shape[1]):
                mse += residual(cos1mdata, i, j, cos1mdata.shape[0] - 1, cos1mdata.shape[1] - 1)
        mse = mse / float((cos1mdata.shape[0]) * (cos1mdata.shape[1]))

        print "new mse", mse
        print cos1mdata.shape
    else:
        print "did not find row cluster lower that delta"
        break

cut = check_columns(cos1mdata)
print "cut column at position:", cut

dfout= pd.DataFrame(cos1mdata)
dfout.to_csv("output.csv")
