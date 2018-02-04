
import pandas as pd
import numpy as np
import csv


df_prottable = pd.read_csv('protwith_tss_1968_shelloperons.csv', delimiter=',',header=None)
prottable = []
for row in df_prottable.iterrows():
    index, data = row
    prottable.append(data.tolist())


def rcomplement(seq):
    temp =""
    for i in range(len(seq)-1,-1,-1):
        if seq[i] =='A': temp+='T'
        elif seq[i] == 'T': temp+= 'A'
        elif seq[i] == 'C': temp+= 'G'
        elif seq[i] == 'G': temp+= 'C'
    return temp


infile1 = open("H37Rv.fna", 'r')
genome =""
for line in infile1:
    words = line.strip()
    if words[0]==">":
        continue
    genome +=words
infile1.close()



class Node(object):
    def __init__(self, data):
        nts=[]
        self.group = 0
        for i in range(0,len(data)):
            if data[i]=='/':
                if data[i+1] not in nts:
                    nts.append(data[i+1])
                continue
            if data[i]=='-':
                self.group = int(data[i+1:i+3])
                break
            if data[i] not in nts:
                nts.append(data[i])

        self.nucleotides = nts
        self.next_node = None

GO ="go"
END = "end"

class FSM(object):
    def __init__(self, seq):
        self.sequence = seq
        self.head = None
        for i in range(0,len(self.sequence)):
            if i ==0:
                newnode = Node(self.sequence[i])
                self.head = newnode
                newnode.next_node = None
            else:
                newnode = Node(self.sequence[i])
                tempnode = self.head
                while tempnode.next_node != None:
                    tempnode= tempnode.next_node
                tempnode.next_node = newnode
                newnode.next_node = None

    def check_seq(self,testseq):
        temp =testseq
        currentnode = self.head
        state = GO
        position = 0
        mismatch = 0
        while(state !=END and currentnode!=None):
            #print currentnode.group
            for i in range(0,currentnode.group):
                if testseq[i] in currentnode.nucleotides:
                    position += 1
                    #print testseq[i], currentnode.nucleotides
                else:
                    mismatch += 1
                    position += 1
                    if mismatch >1:
                        state = END
                        break
            testseq = testseq[position:]
            position = 0
            currentnode = currentnode.next_node

        if state == GO :return True,temp[0:6]
        else: return False,-1



# significant sigmas A-B-C-D-E-F-G-H-J-K

def count_mm(A,B):
    mm = 0
    for i in range(len(A)):
        if A[i]!=B[i]: mm += 1
    return mm

def find_best(window,patt):
    n,m, = len(window),len(patt)
    besti,bestmm = -1,-1
    for i in range(n-m):
        mm = count_mm(window[i:i+m],patt)
        if besti==-1 or mm<bestmm: besti,bestmm = i,mm
    return besti,bestmm



siga_1=["T-2","G-1","C-1","G-1","A-1"]
siga_1_fsm = FSM(siga_1)

siga_2=["T-1","A-1","A/T/G/C-3","T-1"]
siga_2_fsm = FSM(siga_2)



sigb_1=["A/T/G/C-1","G-1","T-1","G-2"]
sigb_1_fsm = FSM(sigb_1)

sigb_2=["A/T/G/C-2","G-1","A/T/G/C-2","G-1"]
sigb_2_fsm = FSM(sigb_2)


sigc_1=["C/G-3", "A-2","T-1"]
sigc_1_fsm = FSM(sigc_1)

sigc_2=["C-1", "G-1", "T-1","C/G-3"]
sigc_2_fsm = FSM(sigc_2)




sigd_0=["G-1","T-1","A-2","C-1","G-1"]
sigd_0_fsm = FSM(sigd_0)

sigd_1=["A-1","G-1","A-3","G-1"]
sigd_1_fsm = FSM(sigd_1)

sigd_2=["C-1","G-1","T-2","A-2"]
sigd_2_fsm = FSM(sigd_2)



sige_1=["G-2","A-2","C-1","C/T-1"]
sige_1_fsm = FSM(sige_1)

sige_2=["G","T-2"]
sige_2_fsm = FSM(sige_2)



sigf_1=["G-2","A/T-2","T-1"]
sigf_1_fsm = FSM(sigf_1)

sigf_2=["G-3","T-1","A-1","C/T-1"]
sigf_2_fsm = FSM(sigf_2)



sigg_1=["G-1","C-1","G-1","A/T/G/C-1","G-1","T-1"]
sigg_1_fsm = FSM(sigg_1)


sigg_2=["C-1","G-1","A-1","A/T/G/C-1","C-1","A-1"]
sigg_2_fsm = FSM(sigg_2)


sigh_1=["G-2","A-2","C/T-1","A-1"]
sigh_1_fsm = FSM(sigh_1)

sigh_2=["G","T-2"]
sigh_2_fsm = FSM(sigh_2)


sigj_1=["G-1","T-1","C-1","A-1","C-1","A-1"]
sigj_1_fsm = FSM(sigj_1)


sigj_2=["C-1","G-1","T-1","C-2","T-1"]
sigj_2_fsm = FSM(sigj_2)



sigk_1=["C-2","A-1","T-1","C-2"]
sigk_1_fsm = FSM(sigk_1)

sigk_2=["C-2","G-1","A-2","T-1"]
sigk_2_fsm = FSM(sigk_2)


sigmas ={}
sigmas["sigA"] =[]
sigmas["sigB"] =[]
sigmas["sigC"] =[]
sigmas["sigD"] =[]
sigmas["sigE"] =[]
sigmas["sigF"] =[]
sigmas["sigG"] =[]
sigmas["sigH"] =[]
sigmas["sigJ"] =[]
sigmas["sigK"] =[]

for j in range(0,len(prottable)):
    profile = prottable[j]
    rvnum = profile[5]
    geneid = profile[4]
    tss = profile[6]
    sign = profile[0]

    if sign == '+':
        tss -= 1
        UTR,leader =  genome[tss-50:tss], genome[tss:tss+10]
    else:
        UTR, leader = rcomplement(genome[tss:tss+50]), rcomplement(genome[tss-10:tss])



    for i in range(0, len(UTR) -25):
        f1, seq1 = siga_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = siga_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigA"]:
                        sigmas["sigA"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigA"



    # for i in range(0, len(UTR) - 6):
    #     f2, seq2 = siga_2_fsm.check_seq(UTR[i:])
    #     if f2 and (i-len(UTR)) >=-15 :
    #             if rvnum not in sigmas["sigA"]:
    #                 sigmas["sigA"].append(rvnum)
    #                 print UTR, "\t", leader, "\t", rvnum[0:10], "\t",  "\t",i-len(UTR),seq2, "\t", "sigA"


    for i in range(0, len(UTR) -25):
        f1, seq1 = sigb_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigb_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigB"]:
                        sigmas["sigB"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigB"


    for i in range(0, len(UTR) -25):
        f1, seq1 = sigc_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigc_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigC"]:
                        sigmas["sigC"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigC"


    for i in range(0, len(UTR) - 25):
        f2, seq2 = sigd_0_fsm.check_seq(UTR[i:])
        if f2 and (i-len(UTR)) >=-15 :
                if rvnum not in sigmas["sigD"]:
                    sigmas["sigD"].append(rvnum)
                    print UTR, "\t", leader, "\t", rvnum[0:10], "\t",  "\t",i-len(UTR),seq2, "\t", "sigD"



    for i in range(0, len(UTR) -25):
        f1, seq1 = sigd_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigd_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigD"]:
                        sigmas["sigD"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigD"



    for i in range(0, len(UTR) -25):
        f1, seq1 = sigc_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 3):
                f2, seq2 = sige_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigE"]:
                        sigmas["sigE"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigE"


    for i in range(0, len(UTR) -25):
        f1, seq1 = sigf_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigf_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigF"]:
                        sigmas["sigF"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigF"

    for i in range(0, len(UTR) -25):
        f1, seq1 = sigg_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigg_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigG"]:
                        sigmas["sigG"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigG"


    for i in range(0, len(UTR) -25):
        f1, seq1 = sigh_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 3):
                f2, seq2 = sigh_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigH"]:
                        sigmas["sigH"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigH"

    for i in range(0, len(UTR) -25):
        f1, seq1 = sigj_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigj_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigJ"]:
                        sigmas["sigJ"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigJ"


    for i in range(0, len(UTR) -25):
        f1, seq1 = sigk_1_fsm.check_seq(UTR[i:])
        if f1 and (i-len(UTR)) >=-40:
            UTRm = UTR[i+21:]
            for j in range(0, len(UTRm) - 6):
                f2, seq2 = sigk_2_fsm.check_seq(UTRm[j:])
                if f2 and ((i+21+j)-len(UTR)) >=-15 :
                    if rvnum not in sigmas["sigK"]:
                        sigmas["sigK"].append(rvnum)
                        print UTR, "\t", leader, "\t", rvnum[0:10], "\t", i-len(UTR),seq1, "\t",(i+21+j)-len(UTR),seq2, "\t", "sigK"



for items in sigmas:
    print items, len(sigmas[items])#, sigmas[items]


sigmafactors = ["sigA","sigB","sigC","sigD","sigE","sigF","sigG","sigH","sigJ","sigK"]
operons = [x[5] for x in prottable]
sigmat = [[0 for x in range(0,len(sigmafactors))] for y in range(0,len(prottable))]

for i in range(0,len(operons)):
    operon = operons[i]
    for j in range(0,len(sigmafactors)):
        sigfs = sigmafactors[j]
        if operon in sigmas[sigfs]:
            sigmat[i][j] = 1


f = open("OUTPUT_sigmafactor_binding_operon_matrix.csv",'wt')
writer = csv.writer(f)
temp = ["operons"] + sigmafactors
writer.writerow(temp)
for i in range(0,len(prottable)):
    temp =[]
    temp.append(operons[i])
    for j in range(0,len(sigmat[i])):temp.append(sigmat[i][j])
    writer.writerow(temp)
f.close()