import numpy as np
import pylab as plt
import pandas as pd


df_cos= pd.read_csv('PeakAssignmentGeneProfile_FINAL.csv', delimiter=',', header=0)
temp_cos1 = []
for row in df_cos.iterrows():
    index, data = row
    temp_cos1.append(data.tolist())
print temp_cos1[0]

df_target= pd.read_csv('TEST4_clara_ftsZ.csv', delimiter=',', header=None)
temp_2 = []
for row in df_target.iterrows():
    index, data = row
    temp_2.append(data.tolist())

temp_2 = [val for sublist in temp_2 for val in sublist]

print temp_2
sig_profile =[]


genenames = [x[0] for x in temp_cos1]

for i in range(0,len(temp_2)):
    rv_name = temp_2[i]
    try:
        geneid = genenames.index(rv_name)
        profile= temp_cos1[geneid]
        sig_profile.append(profile)
    except Exception:
        print"Error - gene not found in profile", rv_name
        #pass
dfout0 = pd.DataFrame(sig_profile )
dfout0.to_csv('clara_FTSZ_profiles.csv')