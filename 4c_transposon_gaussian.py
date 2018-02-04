
########################################### 3:GPy package - test on LFC data #####################################
########################  https://sheffieldml.github.io/GPy/
import pandas as pd
import GPy, numpy as np
import math
from matplotlib import pyplot as plt
import random



f1 = open('avium_combined_read_counts_TTR.txt','r')

dat =[]
for line in f1:
    temp =[]
    data = line.split(' ')
    for l in data:
        if l != '':
            temp.append(l)
    dat.append(temp)

dat2 =[]
for lines in dat:
    temp=[]
    temp.append(int(lines[0]))
    temp.append(float(lines[1]))
    dat2.append(temp)

ta_sites = [x[0] for x in dat2]
counts = [x[1]for x in dat2]
nzmean = np.mean([i for i in counts if i>0.0])
counts = counts - nzmean



###################################################
#plt.plot(tasites,counts)



X = ta_sites
X = np.atleast_2d(X).T
x_pred = X


y = counts
Y = np.atleast_2d(y).T




kernel = GPy.kern.RBF(input_dim=1, variance=10.0, lengthscale=500.0)


############# no overlap batch of 1250 sites

gp_pred = []

arr_length= len(Y)
batch = 1250
rounds = int(arr_length/batch)
start=0
for i in range(rounds):

    print i

    X2 = X[start:start+batch]
    Y2 = Y[start:start+batch]
    x_pred2 = X2
    start = start + batch

    m = GPy.models.GPRegression(X2, Y2, kernel, noise_var=1.0)

    mu=m.predict(x_pred2)

    mu = mu[0]

    for item in mu:
        gp_pred.append(item)

    if len(X)<start +batch:
        X2 = X[start:]
        Y2 = Y[start:]
        x_pred2 = X2

        m = GPy.models.GPRegression(X2, Y2, kernel, noise_var=1.0)

        mu = m.predict(x_pred2)

        mu = mu[0]

        for item in mu:
            gp_pred.append(item)


gp_pred = np.array(gp_pred)

gp_pred_noolap = np.array(gp_pred)

plt.plot(Y,'bo');
plt.plot(gp_pred,'r-')

plt.plot(Y[5000:5500],'bo')
plt.plot(gp_pred[5000:5500],'r-')




################# with overlap batch of 1000 sites with 100 overlap 10% overlap




gp_pred = []

arr_length= len(Y)
batch = 1000
overlap = 100
rounds = int(arr_length/(batch-overlap) )
start=0
for i in range(rounds):

    print i

    X2 = X[start:start+batch]
    Y2 = Y[start:start+batch]
    x_pred2 = X2
    start = start + batch - overlap

    m = GPy.models.GPRegression(X2, Y2, kernel, noise_var=1.0)

    mu=m.predict(x_pred2)

    mu = mu[0]

    for item in mu:
        gp_pred.append(item)

    if len(X)<start +batch-overlap:
        X2 = X[start:]
        Y2 = Y[start:]
        x_pred2 = X2

        m = GPy.models.GPRegression(X2, Y2, kernel, noise_var=1.0)

        mu = m.predict(x_pred2)

        mu = mu[0]

        for item in mu:
            gp_pred.append(item)


gp_pred = np.array(gp_pred)


gp_pred_final =[]
unique = 900
overlap = 100


for i in range(0,len(gp_pred),1000):
    if i == 0:
        temp = gp_pred[i:i+unique]

        olap1 = gp_pred[i+unique: i+unique+overlap]

        olap2 = gp_pred[i + unique + overlap: i + unique + overlap+overlap]

        temp2 = [np.mean(x) for x in zip(olap1, olap2)]

        for num in temp: gp_pred_final.append(num)

        for num2 in temp2: gp_pred_final.append(num2)
    else:
        temp = gp_pred[i+overlap:i+unique]

        olap1 = gp_pred[i+unique: i+unique+overlap]

        olap2 = gp_pred[i + unique + overlap: i + unique + overlap+overlap]

        temp2 = [np.mean(x) for x in zip(olap1, olap2)]

        for num in temp: gp_pred_final.append(num)

        for num2 in temp2: gp_pred_final.append(num2)



gp_pred_final = np.array(gp_pred_final)


plt.plot(Y,'b-')
plt.plot(gp_pred_final,'g-')



# comparsion overlap vs no ol

li = 9000
hi = 10000

plt.plot(Y[li:hi],'o')
plt.plot(gp_pred_noolap[li:hi],'r-')
plt.plot(gp_pred_final[li:hi],'g-')


