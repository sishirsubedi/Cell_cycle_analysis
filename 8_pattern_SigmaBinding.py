
import pandas as pd



df_sigmap = pd.read_csv('sigmapping.csv', delimiter=',',header=None)
sigmap = []
for row in df_sigmap.iterrows():
    index, data = row
    sigmap.append(data.tolist())


df_prottable = pd.read_csv('H37Rv4_prottable.csv', delimiter=',',header=None)
prottable = []
for row in df_prottable.iterrows():
    index, data = row
    prottable.append(data.tolist())

result = []

for i in range(0,len(prottable)):

    gene= prottable[i]
    gene_sign = gene[0]



    for j in range(0,len(sigmap)):

        sigma= sigmap[j]
        sigma_position = int(sigma[3])

        if sigma[2]==' location ' and gene_sign == '+':

            gene_position = int(gene[1])

            if sigma_position < gene_position-26 and sigma_position > gene_position - 45:
                temp=[]
                for items in sigma:
                    temp.append(items)
                for items in gene:
                    temp.append(items)
                result.append(temp)

        elif sigma[2]==' rlocation ' and gene_sign =='-':

            gene_position = int(gene[1])

            if sigma_position > gene_position+26 and sigma_position < gene_position + 45:
                temp=[]
                for items in sigma:
                    temp.append(items)
                for items in gene:
                    temp.append(items)
                result.append(temp)



for items in result:
    print items



