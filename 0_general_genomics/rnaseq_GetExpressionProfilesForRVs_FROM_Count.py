import numpy as np
import pylab as plt
import pandas as pd


df_cos= pd.read_csv('test.csv', delimiter=',', header=None)
temp_cos1 = []
for row in df_cos.iterrows():
    index, data = row
    temp_cos1.append(data.tolist())
print temp_cos1[0]

df_target= pd.read_csv('test2.csv', delimiter=',', header=None)
temp_2 = []
for row in df_target.iterrows():
    index, data = row
    temp_2.append(data.tolist())

temp_2 = [val for sublist in temp_2 for val in sublist]

print temp_2

sig_profile =[]
genenames = [x[0] for x in temp_cos1]
for i in range(0,len(temp_2)):
    rvid = temp_2[i]
    print rvid
    geneid = genenames.index(rvid)
    profile= temp_cos1[geneid]
    temp1 = profile[1:16]
    temp2 = profile[16:31]
    print temp1
    print temp2
    temp = [(x+y)/2. for x,y in zip(temp1,temp2)]
    print temp
    temp.append(profile[0])
    sig_profile.append(temp)

dfout0 = pd.DataFrame(sig_profile )
dfout0.to_csv('TF_profiles.csv')