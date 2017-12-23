
import pandas as pd


###  PEAK FINDING RULES:
# 0. end rule : two ends on both sides 0-1 and 14-15 are not considered for peak finding
# 1. neighbour rule : peak should be greater than two nearest neighbor on both sides
# 3. low expression rule : difference between global minimum and the highest peak (global max)
#       must be less than threshold(0.4)
# --th--0 peaks:
# 0.4   - 953
#  0.5   - 1108
# 0.6   - 1258
# 0.7    - 1413
# 0.8   - 1580
# 4. peak difference rule : the height of second highest peak must be at least 1/2 of first highest peak

def getindex(arr):
    th1=2.0 # for second peak as least as 1/2 high as first one
    th2=0.5 # 0.4 is to remove lower fluctuating profiles
    indx=[]



# add points to indx which is higher than two neighbors on both sides
    for i in range(2,13):
        if(arr[i]>arr[i-1] and arr[i]>arr[i+1] and \
          arr[i] > arr[i-2] and arr[i] > arr[i+2]):
            indx.append(i)

#max is the max peak from indx and min is lowest point in entire array
# clear index if difference between lowest point in array and the highest peak is <0.4
    max = -100.0
    min = 100.0
    #print indx
    for i in range(0,15):
        if arr[i] < min:
            min = arr[i]
    for i in range(0,len(indx)):
        if arr[indx[i]] > max:
            max = arr[indx[i]]

    if max-min<th2:
        indx=[]

#in case of two peaks - make sure both peaks have at least half of each others height
#if not then delete second peak from the indx array
    if len(indx)>1:
        if arr[indx[0]] > arr[indx[1]]:
            peak1 =  arr[indx[0]]
            peak2 =  arr[indx[1]]
            if(peak2 < peak1/th1):
                del indx[1]
        else:
            peak1 = arr[indx[1]]
            peak2 = arr[indx[0]]
            if (peak2 < peak1 / th1):
                del indx[0]

    return indx


df_cos= pd.read_csv('final_cos_smooth.csv', delimiter=',', header=0)

#df_cos= pd.read_csv('final_rv_smooth.csv', delimiter=',', header=0)


temp_cos1 = []
for row in df_cos.iterrows():
    index, data = row
    temp_cos1.append(data.tolist())

numofpeaks=[]
zeropeaks=[]
onepeaks =[]
twopeaks =[]
threepeaks =[]
fourpeaks =[]
divisomepeaks =[]
for genes in temp_cos1:
    ygene=genes[1:]
    peaks = getindex(ygene)
    length=len(peaks)
    if length not in numofpeaks:
        numofpeaks.append(length)
    if length == 0:
        genes.append(0)
        genes.append("#")
        genes.append("#")
        zeropeaks.append(genes)
    if length == 1:
        genes.append(1)
        genes.append(peaks[0]+1)
        genes.append("#")
        onepeaks.append(genes)
    elif length ==2 :
        ### idea of primary and secondary peak ####
        ppeak =0
        speak =0
        if (ygene[peaks[0]] )>(ygene[peaks[1]]):
            ppeak, speak = (peaks[0]), (peaks[1])
        else:
            ppeak, speak= (peaks[1]), (peaks[0])
        genes.append(2)
        genes.append(ppeak+1)  # +1 to adjust for zero array
        genes.append(speak+1)
        twopeaks.append(genes)
    elif length ==3:
        genes.append(3)
        genes.append(peaks[0]+1)  # +1 to adjust for zero array
        genes.append(peaks[1]+1)
        #genes.append(peaks[2]+1)
        threepeaks.append(genes)
    elif length == 4:
        fourpeaks.append(genes)

    # if len([i for i, j in zip(divisometp, peaks) if i == j])>=1:
    #     print [i for i, j in zip(divisometp, peaks) if i == j]
    if length ==1:
        for i in peaks:
            if i in range(11,12):
                if genes not in divisomepeaks:
                    divisomepeaks.append(genes)

print "peak diversity " , numofpeaks
print "genes with 0 peak - ", len(zeropeaks)
print "genes with 1 peak - ", len(onepeaks)
print "genes with  2 peaks - ", len(twopeaks)
print "genes with 3 peaks - ", len(threepeaks)
print "genes with 4 peaks - ", len(fourpeaks)
print "genes with one peak at timepoint 12", len(divisomepeaks)

#fourpeaks = [val for sublist in fourpeaks for val in sublist]

for items in onepeaks:
    zeropeaks.append(items)

for items in twopeaks:
    zeropeaks.append(items)

for items in threepeaks:
    zeropeaks.append(items)


dfout0 = pd.DataFrame(zeropeaks )

dfout0.columns = ["ORF ID", "3hrs" , "6.5hrs" , "9hrs" , "12hrs" , "18.5hrs" , "21hrs" , "27hrs" , \
                  "31hrs" , "33hrs" , "36hrs" , "39.5hrs" , "42hrs" , "45.5hrs" , "52hrs" , "55hrs",\
                  'Number of Peaks',"Primary Peak", "Secondary Peak"]

#dfout0.to_csv('COS_peakResult.csv',index=False)


#dfout0.to_csv('RV_peakResult.csv',index=False)