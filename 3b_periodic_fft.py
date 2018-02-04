import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from scipy.fftpack import fft,ifft

xdata = np.array([0.,1.,2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.])


df3 = pd.read_csv('cos_smooth_divbymean_adjtp_eventp.csv', delimiter=',',header=0)



temp3 = []
for row in df3.iterrows():
    index, data = row
    temp3.append(data.tolist())

fft_out =[]
for i in range(0,len(temp3)):
    gene = temp3[i]
    yf = np.fft.fft(gene[1:])
    freq = np.fft.fftfreq(15, 1)
    pfs = sp.where(freq > 0)  # Select postive frequencies
    freqs = freq[pfs]
    power = abs(yf)[pfs] ** 2  # abs(yf) is getting magnitude value of complex numbers
    rvname =gene[0]
    temp=[]
    temp.append(rvname)
    # for element in freqs:
    #     temp.append(element)
    for element in power:
        temp.append(element)

    fft_out.append(temp)

dfout0 = pd.DataFrame(fft_out )
dfout0.to_csv('FFT_profiles.csv')