"""
Created on Sun Jan 30 21:55:27 2022

@author: bartek
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy.signal import butter, lfilter, freqz
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy import signal
from scipy.integrate import trapz
from scipy.interpolate import interp1d

def timedomain(rr):

    hr = 60000/rr
    
    df2['Mean RR (ms)'] = np.mean(rr)
    df2['STD RR/SDNN (ms)'] = np.std(rr)
    df2['Mean HR (beats/min)'] = np.mean(hr)
    df2['STD HR (beats/min)'] = np.std(hr)
    df2['RMSSD (ms)'] = np.sqrt(np.mean(np.square(np.diff(rr))))
    df2['pNNxx (%)'] = 100 * np.sum((np.abs(np.diff(rr)) > 50)*1) / len(rr)
    df2['Min HR (beats/min)'] = np.min(hr)
    df2['Max HR (beats/min)'] = np.max(hr)

def frequency_domain(rri, fs):
    # Estimate the spectral density using Welch's method
    
    '''
    Segement found frequencies in the bands 
     - Very Low Frequency (VLF): 0-0.04Hz 
     - Low Frequency (LF): 0.04-0.15Hz 
     - High Frequency (HF): 0.15-0.4Hz
    '''
    cond_vlf = (fxx >= 0) & (fxx < 0.04)
    cond_lf = (fxx >= 0.04) & (fxx < 0.15)
    cond_hf = (fxx >= 0.15) & (fxx < 0.4)
    
    # calculate power in each band by integrating the spectral density 
    vlf = trapz(pxx[cond_vlf], fxx[cond_vlf])
    lf = trapz(pxx[cond_lf], fxx[cond_lf])
    hf = trapz(pxx[cond_hf], fxx[cond_hf])
    
    # sum these up to get total power
    total_power = vlf + lf + hf

    # find which frequency has the most power in each band
    peak_vlf = fxx[cond_vlf][np.argmax(pxx[cond_vlf])]
    peak_lf = fxx[cond_lf][np.argmax(pxx[cond_lf])]
    peak_hf = fxx[cond_hf][np.argmax(pxx[cond_hf])]

    # fraction of lf and hf
    lf_nu = 100 * lf / (lf + hf)
    hf_nu = 100 * hf / (lf + hf)
    
    df2['Power VLF (ms2)'] = vlf
    df2['Power LF (ms2)'] = lf
    df2['Power HF (ms2)'] = hf   
    df2['Power Total (ms2)'] = total_power

    df2['LF/HF'] = (lf/hf)
    df2['Peak VLF (Hz)'] = peak_vlf
    df2['Peak LF (Hz)'] = peak_lf
    df2['Peak HF (Hz)'] = peak_hf

    df2['Fraction LF (nu)'] = lf_nu
    df2['Fraction HF (nu)'] = hf_nu
    
#%%

df = pd.read_csv(r'/Volumes/ExFAT-EMTEC/VNS/VNS-sampled/006/VNS_ON.csv') # drop file

# convert from ms to s
df['RRI'] = df['RRI']*1000

# Binning by SBP
BPbin= []

for x in df.SBP:
    BPbin.append(int((x - min(df.SBP))/3))
    
df['bin'] = BPbin
df = df[:(len(df))-2] #removes trailing incomplete cardiac cycle
print(df)

groupedRR = df['RRI'].groupby(df['bin'])
RRarray = groupedRR.mean() 

groupedSBP = df['SBP'].groupby(df['bin'])
SBParray = np.asarray(groupedSBP.mean())
SBParray2 = groupedSBP.mean()

bin_weight = groupedSBP.size()/df.index.max()
df2 = df[['RRI','SBP']].mean()
#df2.drop(['bin'],axis=0)

slope, intercept, r_value, p_value, std_err = linregress(SBParray, RRarray)
df2['BRS slope'] = slope
df2['R^2'] = r_value**2

sns.regplot(x = SBParray, y = RRarray)
plt.title('BRS')
plt.xlabel('Brachial SBP (mmHg)')
plt.ylabel('RRI')
plt.xlim(SBParray.min()-1,SBParray.max()+1)
plt.ylim(RRarray.min()-50,RRarray.max()+50)
#plt.text(100,1240, df2['R^2'])
plt.text(SBParray.min()+1 ,RRarray.max()+35, r'$Slope: %s,\ R^2=%s$' % (np.round(slope, decimals =2), np.round(df2['R^2'], decimals = 2)))

#time axis in s
x_ecg = np.cumsum(df['RRI'])/1000

#fit function to the dataset
f_ecg = interp1d(x_ecg, df['RRI'], kind='cubic', fill_value= 'extrapolate')

#sample rate for interpolation
fs = 4
steps = 1 / fs

#sample using the interpolation function
xx_ecg = np.arange(0, np.max(x_ecg), steps)
rr_interpolated_ecg = f_ecg(xx_ecg)

plt.figure(figsize=(20,8))

plt.subplot(211)
plt.title('rr-intervals')
plt.plot(x_ecg, df['RRI'], color='k', markerfacecolor='#A999D1',marker='o', markersize=3)
plt.xlabel('Time (s)')
plt.ylabel('rr-interval (ms)')

plt.subplot(212)
plt.title('rr-intervals (cubic interpolation)')
plt.plot(xx_ecg, rr_interpolated_ecg, color='r')
plt.xlabel('Time (s)')
plt.ylabel('RR-interval (ms)')
plt.show()

fxx, pxx = signal.welch(x=rr_interpolated_ecg, fs=fs, nperseg=256)

#fit a function for plotting bands
powerspectrum_f = interp1d(fxx, pxx, kind='cubic', fill_value= 'extrapolate')

plt.figure(figsize=(15,6))
#plt.plot(fxx,pxx,color='k',linewidth=0.5)
plt.title("FFT Spectrum (Welch's periodogram)")

# setup frequency bands for plotting
x_VLF = np.linspace(0, 0.04, 100)
x_LF = np.linspace(0.04, 0.15, 100)
x_HF = np.linspace(0.15, 0.4, 100)

plt.gca().fill_between(x_VLF, powerspectrum_f(x_VLF), alpha=0.2, color="#F5866F", label="VLF")
plt.gca().fill_between(x_LF, powerspectrum_f(x_LF), alpha=0.2, color="#51A6D8", label="LF")
plt.gca().fill_between(x_HF, powerspectrum_f(x_HF), alpha=0.2, color="#17BC0F", label="HF")

plt.gca().set_xlim(0, 0.5)
plt.gca().set_ylim(0)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Density")

plt.legend()
plt.show()

timedomain(df['RRI'])
frequency_domain(df['RRI'],fs)

#df2 = df2.drop(['bin'],axis=0)

df2.to_clipboard(excel=True)
df2







