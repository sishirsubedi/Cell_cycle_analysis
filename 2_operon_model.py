# this algorithm models and compares operons from our data vs published literature



import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import random

df_cos= pd.read_csv('1_mean_expression_cos_all4018_withgenename.csv', delimiter=',', header=None)
profile = [d1.tolist() for i1,d1 in df_cos.iterrows()]
genenames =[x[0] for x in profile]


df_prottable = pd.read_csv('H37Rv4_prottable.csv', delimiter=',',header=None)
prottable =[d2.tolist() for i2,d2 in df_prottable.iterrows()]

# pos_prot = [x for x in prottable if x[0]=='+']
# neg_prot = [x for x in prottable if x[0]=='-']



def calculate_cor_score(datamat):
    cormat = np.corrcoef(datamat)
    cormat = np.matrix(cormat)
    cormat_mod = np.triu(cormat, 1)
    score = np.sum(cormat_mod) / np.count_nonzero(cormat_mod)
    return score






df_operons= open('shell_operons.csv', 'r')
shell_operons = []
for line in df_operons:
    word = line.strip().split(",")
    shell_operons.append(list(word))

temp=[]
for op in shell_operons:
    if op not in temp:
        temp.append(op)

shell_operons = temp

pos_profile =[]
neg_profile =[]
for i in range(0,len(profile)):
    gene = profile[i]

    pinfo = [x for x in prottable if gene[0]==x[0]].pop()

    if pinfo[1] == '+':
        gene.extend(pinfo)
        pos_profile.append(gene)
    else:
        gene.extend(pinfo)
        neg_profile.append(gene)


#for items in neg_profile:print items


threshold = 0.4
dis_th =100


while threshold < 0.5:

    threshold += 0.1

    operons =[]

    i=0


    while i < len(pos_profile):

        gene1 = pos_profile[i]
        rvnum = gene1[0]#+ '/'+gene1[1]
        genedata = gene1[2:17]

        if i == len(pos_profile )-1:
            operons.append([rvnum])
            i += 1

        else:
            gene2 = pos_profile[i+1]
            rvnum2 = gene2[0] #+ '/'+gene2[1]
            genedata2 = gene2[2:17]

            datamat = []
            score = 0.0
            datamat.append(genedata)
            datamat.append(genedata2)
            score = calculate_cor_score(datamat)

            if gene2[19] <= gene1[20] + dis_th and  score >= threshold:

                temp =[]
                temp.append(rvnum)
                temp.append(rvnum2)
                i += 1

                if i == len(pos_profile)-1:
                    operons.append(temp)

                else:

                    while score >= threshold and i< len(pos_profile)-1:

                            oldgene = pos_profile[i]

                            i += 1
                            newgene = pos_profile[i]
                            newrv = newgene[0]#+ '/'+newgene[1]
                            newgenename = newgene[1]
                            newgenedata = newgene[2:17]

                            datamat.append(newgenedata)

                            score = calculate_cor_score(datamat)

                            if newgene[19] <= oldgene[20] + dis_th and score >= threshold and i< len(pos_profile)-1:

                                temp.append(newrv)

                            else:
                                operons.append(temp)
                                break
            else:
                operons.append([rvnum])
                i += 1



    i = len(neg_profile)-1



    while i >= 0:

        gene1 = neg_profile[i]
        rvnum = gene1[0]#+ '/'+gene1[1]
        genedata = gene1[2:17]


        if i == 0:
            operons.append([rvnum])
            i -= 1

        else:
            gene2 = neg_profile[i-1]
            rvnum2 = gene2[0]# + '/'+gene2[1]
            genedata2 = gene2[2:17]



            datamat = []
            score = 0.0
            datamat.append(genedata)
            datamat.append(genedata2)
            score = calculate_cor_score(datamat)

            if gene2[20] >= gene1[19] - dis_th and score >= threshold:

                temp =[]

                temp.append(rvnum)

                temp.append(rvnum2)

                i = i-1

                if i ==0 :
                    operons.append(temp)
                    i -= 1

                else:
                    while score >= threshold and i> 0:

                            oldgene = neg_profile[i]

                            i -= 1

                            newgene = neg_profile[i]
                            newrv = newgene[0]#+ '/'+newgene[1]
                            newgenename = newgene[1]
                            newgenedata = newgene[2:17]

                            datamat.append(newgenedata)

                            score = calculate_cor_score(datamat)

                            if newgene[20] >= oldgene[19] - dis_th and score >= threshold and i>0:

                                temp.append(newrv)

                            else:
                                operons.append(temp)
                                break


            else :
                operons.append([rvnum])
                i -= 1

    #
    # msum = 0
    # for soperon in shell_operons:
    #     if len(soperon) == 1:
    #         continue
    #     else:
    #         for newoperon in operons:
    #             if len(soperon) == len(newoperon) and len(soperon)== len([x for x in newoperon if x in soperon]):
    #                 #print soperon, newoperon
    #                 msum += 1
    #                 break


    msum = 0
    for i in range(0,len(genenames)-1):
        operon_pair = [genenames[i],genenames[i+1]]
        if operon_pair in shell_operons and operon_pair in operons:
            msum += 1
        else:
            completed = False
            for soperon in shell_operons:
                if operon_pair[0] in soperon and operon_pair[1] in soperon:
                    for noperon in operons:
                        if operon_pair[0] in noperon and operon_pair[1] in noperon and soperon[0]==noperon[0]:
                            msum += 1
                            completed = True
                            print soperon , noperon
                            break
                if completed:
                    break



        # else:
        #     for j in range(0,len(shell_operons)-1):
        #         if operon_pair[0] in shell_operons[j] and operon_pair[1] in shell_operons[j+1]:
        #             for k in range(0,len(operons)-1):
        #                 if operon_pair[0] in operons[k] and operon_pair[1] in operons[k + 1]:
        #                     msum += 1
        #





    print threshold, len(operons), msum



dfoperons = pd.DataFrame(operons)
dfoperons.to_csv('cormodel_operons.csv')



#
# operons_2more = [x for x in operons if len(x)>1]
#
# shell_operons_2more = [x for x in shell_operons if len(x)>1]
#

oplen = []
max = [0,0]
for items in operons:
    n = len(items)
    oplen.append(n)
    if n>=8:
        print items, len(items)
    if n>max[0]:
        max[0] = len(items)
        max[1] = items


print max



ol  = pd.value_counts(oplen)
x = ol.values
y = np.array(ol.index)
plt.bar(y[:],x[:] )
plt.xlim(1,12)



shell_oplen = []
max = [0,0]
for items in shell_operons:
    n = len(items)
    shell_oplen.append(n)
    if n>1:
        print items, len(items)
    if n>max[0]:
        max[0] = len(items)
        max[1] = items




shell_ol  = pd.value_counts(shell_oplen)
shell_x = shell_ol.values
shell_y = np.array(shell_ol.index)
plt.bar(shell_y[:],shell_x[:])
plt.xlim(1,12)


oldf = pd.DataFrame(ol)
oldf.to_csv('ol.csv')
shell_oldf = pd.DataFrame(shell_ol)
shell_oldf.to_csv('shell_ol.csv')
