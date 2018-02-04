import pandas as pd
import GPy, numpy as np
from matplotlib import pyplot as plt
from scipy import stats

## input data


df_cos1 = pd.read_csv('4_cos1_logmean.csv', delimiter=',',header=0)
temp_cos1 = []
for row in df_cos1.iterrows():
    index, data = row
    temp_cos1.append(data.tolist())

df_cos2 = pd.read_csv('4_cos2_logmean.csv', delimiter=',',header=0)
temp_cos2 = []
for row in df_cos2.iterrows():
    index, data = row
    temp_cos2.append(data.tolist())

df3 = pd.read_csv('cos_sinfit_adjtp_datapoints.csv', delimiter=',',header=0)
temp_sinfit = []
for row in df3.iterrows():
    index, data = row
    temp_sinfit.append(data.tolist())
genenames_sin = [x[0] for x in temp_sinfit]


############## GP Fit

cos_smooth = []

kernel = GPy.kern.RBF(input_dim=1, variance=1.0, lengthscale=20.0)

X = [3.0, 6.5, 9.0, 12.0, 18.5, 21.0, 27.0, 31.0, 33.0, 36.0, 39.5, 42.0, 45.5, 52.0, 55.0, \
     3.0, 6.5, 9.0, 12.0, 18.5, 21.0, 27.0, 31.0, 33.0, 36.0, 39.5, 42.0, 45.5, 52.0, 55.0]

X = np.atleast_2d(X).T
x_pred = X[0:15]

for rv in range(0,len(temp_cos1)):
    ycos1= temp_cos1[rv]
    ycos2 = temp_cos2[rv]
    y = ycos1[1:] + ycos2[1:]
    Y = np.atleast_2d(y).T
    m = GPy.models.GPRegression(X, Y,kernel,noise_var=0.1)
    mu,C=m.predict(x_pred, full_cov=True)
    cos_smooth.append(mu.tolist())


    #figure
    m.plot(legend=None)
    ymin =  min(Y)
    ymax =  max(Y)
    plt.ylim([ymin-0.25 , ymax+0.25]);plt.xlim(0,56)
    #plt.plot(x_pred,mu,'b-',linewidth=2.5,label="cosGP")
    #plt.plot(x_pred,mu,'ro',label ="GP")
    plt.plot(x_pred,ycos1[1:],'c--', linewidth=1.0,label="cos1")
    plt.plot(x_pred,ycos2[1:],'y--', linewidth=1.0,label="cos2")

    gene=ycos1[0];y_sin=[]
    if gene in genenames_sin:
        geneid_sin = genenames_sin.index(gene)
        y_sin = temp_sinfit[geneid_sin]

    if len(y_sin) > 2:
        plt.plot(x_pred, y_sin[1:], 'm--', linewidth=1.0, label='Sinfit')
        plt.plot(x_pred, y_sin[1:], 'mo')
        lgd = plt.legend(["GP-Mean", "Cos1", "Cos2", "Sinfit"], loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12})
        plt.plot(x_pred, ycos1[1:], 'co')
        plt.plot(x_pred, ycos2[1:], 'yo')
        plt.xlabel('Time-points(hr)', fontsize=12)
        plt.ylabel('Relative Expression (mean:1)', fontsize=12)
        plt.title('Gaussian Process (GP) regression - ' +str(gene) , fontsize=12)
    else:
        plt.plot(x_pred, ycos1[1:], 'co')
        plt.plot(x_pred, ycos2[1:], 'yo')
        lgd = plt.legend(["GP-Mean", "Cos1", "Cos2"], loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12})
        plt.xlabel('Time-points(hr)', fontsize=12)
        plt.ylabel('Relative Expression (mean:1)', fontsize=12)
        plt.title('Gaussian Process (GP) regression - ' +str(gene) , fontsize=12)

    plt.savefig(str(gene) + ".png", format="PNG",bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close('all') # to remove figure from memory