import sys
import pandas as pd
import random
import csv
import numpy as np
import scipy.stats as ss
import math
from sklearn import linear_model
from copy import copy


random.seed(1)



df_profile = pd.read_csv('2704_geneswith_peakResult_cos_GPSmooth.csv', delimiter=',',header=0)
profile = []
for row in df_profile.iterrows():
    index, data = row
    profile.append(data.tolist())

siggenes = [x[0] for x in profile]
profilemat = df_profile.values
profilematdata = profilemat[:,7:22]

ematrix = np.array([[0.0 for x in range(0, 15)] for y in range (0,len(siggenes))])

for x in range (0,len(siggenes)):
    for y in range(0, 15):
        ematrix[x,y] = float(profilematdata[x,y])
print " ematrix before centering :- max:", ematrix.max(), "min: ", ematrix.min()

ematrix2 = copy(ematrix)
ematrix2t = ematrix2.T
for y in range(0,2704):
    mean = ematrix2t[:,y].mean()
    std  = math.sqrt(ematrix2t[:,y].var())
    for x in range(0,15):
        ematrix2t[x, y] = (ematrix2t[x,y] - mean )/std

for y in range(0,15):
    mean = ematrix[:,y].mean()
    std  = math.sqrt(ematrix[:,y].var())
    for x in range(0,2704):
        ematrix[x, y] = (ematrix[x,y] - mean )/std

print "~~~~ematrix ~~~~~"
print "shape:", ematrix.shape
print "max:" , ematrix.max(), "min: " , ematrix.min()
print "mean:" , ematrix.mean(), "var: " , ematrix.var()
print "shape2:", ematrix2t.shape
print "max2:" , ematrix2t.max(), "min2: " , ematrix2t.min()
print "mean2:" , ematrix2t.mean(), "var2: " , ematrix2t.var()


f = open("ematrix_observed.csv",'wt')
for items in ematrix:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()



##################################################################

df_paper = pd.read_csv('2_OUTPUT_TFwithTargetPval_AND_Regulation_nomin_filtered.csv', delimiter=',',header=None)
paper = []
for row in df_paper.iterrows():
    index, data = row
    paper.append(data.tolist())

print "Paper record: ", len(paper)

sigtf=[]
for i in range(0,len(paper)):
    current = paper[i]
    if current[0] in siggenes and current[0] not in sigtf:
        sigtf.append(current[0])

print "sig genes : ", len(siggenes)
print "sig tfs: ", len(sigtf)


########## initialize amatrix from paper


amatrix= np.array([[0.0 for x in range(0, len(sigtf))] for y in range (0,len(siggenes))])



for k in range (0, len(paper)):
    interaction = paper[k]
    tf = interaction[0]
    gene= interaction[1]
    strength = float(interaction[3])
    if gene in siggenes and tf in sigtf:
        geneindex = siggenes.index(gene)
        tfindex = sigtf.index(tf)
        amatrix[geneindex,tfindex] = strength

#
# for y in range(0,132):
#     mean = amatrix[:,y].mean()
#     std  = math.sqrt(amatrix[:,y].var())
#     for x in range(0,2704):
#         amatrix[x, y] = (amatrix[x,y] - mean )/std


print "~~~~amatrix ~~~~~"
print "shape:", amatrix.shape
print "max:", amatrix.max(), "min: ", amatrix.min()
print "mean:" , amatrix.mean(), "var: " , amatrix.var()




f = open("amatrix_paper.csv",'wt')
for items in amatrix:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()

print "done"


############################################


df_profile_gpsmooth = pd.read_csv('7_final_sig_2704_cos_smooth.csv', delimiter=',',header=0)
smooth = []
for row in df_profile_gpsmooth.iterrows():
    index, data = row
    smooth.append(data.tolist())
smoothgenes =[x[0] for x in smooth]

tfdata=[]
for i in range(0,len(sigtf)):
    tf= sigtf[i]
    tfindex = smoothgenes.index(tf)
    tfprofile = smooth[tfindex]
    tfdata.append(tfprofile[1:])

pmatrix= np.array([[0.0 for x in range(0, 15)] for y in range (0,132)])


tfdata = np.array(tfdata)
print "tfdata shape", tfdata.shape
for x in range (0,132):
    for y in range(0, 15):
        if tfdata[x,y]<0.0:
            pmatrix[x, y]=0.0
        else:
            pmatrix[x,y] = float(tfdata[x,y])


print "~~~~pmatrix ~~~~~"
print "shape:", pmatrix.shape
print "max:", pmatrix.max(), "min: ", pmatrix.min()
print "mean:", pmatrix.mean(), "var: ", pmatrix.var()


f = open("pmatrix_paper_raw.csv",'wt')
for items in pmatrix:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()

print "done"

#exit()
# ########################


print " ~~~~~~~ ~~~~~~~ initial amatrix and pmatrix ~~~~~~~ ~~~~~~"

print "amatrix zero : ", str( 356928 - np.count_nonzero(amatrix))
print "pmatrix zero : ", str(1980 - np.count_nonzero(pmatrix))

print " ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~"


##########################  CREATE pie matrix ###############################################



def sigmoid(x,k):
    return 1 /(1+(math.e**((-x)*k)))


pie_matrix= np.array([[0.0 for x in range(0, len(sigtf))] for y in range (0,len(siggenes))])



for x in range(0, len(siggenes)):
    for y in range (0,len(sigtf)):
        #pie_matrix[x, y] = sigmoid(amatrix[x,y],1)
        pie_matrix[x, y] = abs(math.tanh(amatrix[x,y]))


print "~~~~piematrix ~~~~~"
print "shape:", pie_matrix.shape
print "max:", pie_matrix.max(), "min: ", pie_matrix.min()



pie_matrix_t = pie_matrix.T



# shooting algorithm


def solve_shooting(y,X,lmbda1, lmbda2,gene):

    #initialize w with solution of ridge regression
    # beta0 = np.dot(X.T,X) + lmbda1*np.matrix(np.identity(132))
    # beta = np.dot (beta0.I , np.dot(X.T,y))

    #initialize beta as zeros
    beta = [[ 0.0 for x0 in range(0, 1)] for y0 in range(0,132)]
    beta =np.array(beta)

    found = False
    tolerance = 1e-1

    while (not found):
        betaold = beta # one gene interaction with all tfs
        for j in range(0,132): # each interaction with TF one at a time

            # get the ith column
            xi =X[:,j] # take one TF state in all time points
            zi = np.dot(xi.T,xi)
            yhat = np.array(np.dot(X, beta))
            yhat = np.array(yhat).flatten()
            yhatwj = np.dot(xi,beta.item(j))
            # calculate the residual without ith w
            yi = (y- yhat)+ yhatwj
            # calculate xi'*yi and see where it falls # unit value describing effect of ith w
            w_j =  np.dot(xi.T,yi)


            #get pieij from pie_matrix
            pie_col =pie_matrix_t[:,gene]
            pie_ij = pie_col[j]
            pie_log =0.0
            if pie_ij != 0.0:
                pie_log = -1 * math.log(pie_ij,2)

            gamma = lmbda1 * pie_log

            #print "wj", w_j, "gamma", gamma,"pieij", pie_ij, "pielog", pie_log, "beta", beta[j]

            if (w_j < -gamma):
                beta[j] = ((w_j + gamma)) / (1 + lmbda2)
            elif (w_j > gamma):
                beta[j] = ((w_j - gamma)) / (1 + lmbda2)
            else:
                beta[j] =0.0



        diff= (np.absolute(beta - betaold)).max()

        if diff <= tolerance :
            found = True

    return beta

# ################################################################################

th = sys.maxint
k=0
a0 = copy(amatrix.T)
p0 = copy(pmatrix)
log = []
while k<5:


    k +=1

    log.append(" ----- optimizing p and a matrix ------ "+ "CYCLE : "+ str(k))


    pn = [[0.0 for x0 in range(0, 1)] for y0 in range (0,132)]


    x1 = copy(amatrix)
    for ci in range(0, 132):
        if np.count_nonzero(x1[:,ci]) != 0:
            x1[:,ci] = ss.zscore(x1[:,ci])


    for i1 in range(0,15):

        y1= copy(ematrix[:,i1])
        #model1 = linear_model.LassoLars(alpha=0.0001,fit_intercept=False,positive=True)
        
        model1 = linear_model.Lars(fit_intercept=False)
        model1.fit(x1,y1)
        temp1 = model1.coef_
        #print "lar coef", temp1
        pn = np.c_[pn, temp1]

    pn = pn[:,1:]
    Pdiffs = abs(pn-pmatrix)
    pmatrix = pn

    col = 0
    for i in range(0, 15):
        if (np.count_nonzero(pmatrix[:,i])== 0):
            col = col + 1
    log.append(" CYCLE: "+ str(k)+ "  ZEROS columns: pmatrix: "+ str(col))


    an = [[0.0 for x2 in range(0, 1)] for y2 in range (0,132)]

    x3 = copy(pmatrix.T)

    for i2 in range(0, 132):
        if np.count_nonzero(x3[:, i2]) != 0:
            x3[:, i2] = ss.zscore(x3[:, i2])

    for i3 in range(0,2704):

        y3 = copy(ematrix2t[:,i3])

        # with prior
        #temp2 = solve_shooting(y3, x3, 10.0, 10.0, i)
        
        
        #model2 = linear_model.LassoLars(alpha=0.0001, fit_intercept=False)
        
        model2 = linear_model.Lars(fit_intercept=False)
        model2.fit(x3, y3)
        temp2 = model2.coef_
        #print temp2
        an = np.c_[an, temp2]

    an= an[:,1:]
    an2 = an.T

    Adiffs = abs(an2 - amatrix)
    amatrix = an2

    col = 0
    for i in range(0, 132):
        if np.count_nonzero(amatrix[:,i])== 0:
            col = col + 1
    log.append(" CYCLE: "+ str(k)+ "  ZEROS columns: amatrix: "+ str(col))

    log.append("Adiff : " + str(Adiffs.max()))
    log.append("Pdiff : " + str(Pdiffs.max()))
    log.append("amatrix zero : "+ str(356928 - np.count_nonzero(amatrix)))
    log.append("pmatrix zero : "+ str(1980 - np.count_nonzero(pmatrix)))

    j = sigtf.index("Rv0023")

    log.append(pmatrix[j, 0:15])

    th = min(Adiffs.max(),Pdiffs.max())

    log.append("current min : " +str(th))

for records in log:
    print records



####################################################
#         Analysis of Result
####################################################




ematrix_hat = np.dot(amatrix, pmatrix)
# for y in range(0,15):
#     mean = ematrix_hat[:,y].mean()
#     std  = math.sqrt(ematrix_hat[:,y].var())
#     for x in range(0,2704):
#         ematrix_hat[x, y] = (ematrix_hat[x,y] - mean )/std

f = open("ematrix_hat.csv",'wt')
for items in ematrix_hat:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()



# for y in range(0,132):
#     if np.count_nonzero(amatrix[:, y]) != 0:
#         mean = amatrix[:,y].mean()
#         std  = math.sqrt(amatrix[:,y].var())
#         for x in range(0,2704):
#             amatrix[x, y] = (amatrix[x,y] - mean )/std
f = open("amatrix_hat.csv",'wt')
for items in amatrix:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()



f = open("pmatrix_hat.csv",'wt')
for items in pmatrix:
     writer = csv.writer(f)
     writer.writerow(items)
f.close()


print "optimization complete...printing result.."

profileforgraph = []
for i in range(0,len(sigtf)):
    tf = sigtf[i]
    for j in range(0,len(siggenes)):
        target = siggenes[j]
        value = amatrix[j,i]
        if abs(value) > 0.0:
            temp=[]
            temp.append(tf)
            temp.append(target)
            temp.append(value)
            profileforgraph.append(temp)

df = pd.DataFrame(profileforgraph, columns=['tf', 'target','w'])
result = np.array(df)

resultcompare =[]
for i in range(0,len(paper)):
    current = paper[i]
    pval = current[2]
    if pval <0.01:
        tf = current[0]
        target =current[1]
        for j in range (0,len(result)):
            currentresult = result[j]
            if tf == currentresult[0] and target == currentresult[1]:
                temp = []
                temp.append(tf)
                temp.append(target)
                temp.append(pval)
                temp.append(current[3])
                temp.append(currentresult[0])
                temp.append(currentresult[1])
                temp.append(currentresult[2])
                resultcompare.append(temp)
                print temp

df_paper_r = pd.read_csv('paperresult_pvallessthan0_01_2.csv', delimiter=',',header=None)
paper_r = []
for row in df_paper_r.iterrows():
    index, data = row
    paper_r.append(data.tolist())

paper_r_tf = [x[0] for x in paper_r]




print "\n\n"
print "~~~~~~~~~~~~~~~~~Result~~~~~~~~~~~~~~\n"

tffreq = [x[0] for x in resultcompare]
tffreqtable = {x: tffreq.count(x) for x in tffreq}

print "Total orig TFs:132  Model TFs: ", len(tffreqtable), "\n"
print "     o-tf  o-tar      m-tf    m-match   %m \n"

for w in sorted(tffreqtable, key=tffreqtable.get, reverse=True):

     tfi = paper_r_tf.index(w)

     original = paper_r[tfi]
     originaltf = original[0]
     originalinteractions = original[1]

     modeltf = w
     modelmatched = tffreqtable[w]

     percentmatch = float(modelmatched)/float(originalinteractions) * 100

     percentmatch = float("{0:.2f}".format(percentmatch))

     print "\t", originaltf , originalinteractions," --- ",  modeltf,  modelmatched,"---",  percentmatch ,"%"



