import sys
import pandas as pd

df_cos1 = pd.read_csv('1_meanexpression_cos_sig.csv', delimiter=',',header=0)
temp_cos1 = []
for row in df_cos1.iterrows():
    index, data = row
    temp_cos1.append(data.tolist())

of1 = open('temp_html','w')

for i in range(0,len(temp_cos1)):
    profile = temp_cos1[i]
    des=profile[2]
    of1.write(
        "<TD>" + profile[0] + "<TD>" + profile[1] + "<TD>" + \
       "<A HREF =" + "images3/" +  str(profile[0]) + ".svg" + ">" + str(profile[0]) + "</A><TR>" + '\n')

of1.close()



