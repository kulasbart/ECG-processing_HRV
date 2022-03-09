import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy.signal import butter, lfilter, freqz
from scipy.stats import linregress

path = r'/volumes/ExFAT-EMTEC/VNS/VNS-BRSanalysis/001/ihg_OFF.csv'
df = pd.read_csv(path)
#df.head()

# Binning by SBP

BPbin= []

for x in df.brachial_SBP:
    BPbin.append(int((x - min(df.brachial_SBP))/3))
    
df['bin'] = BPbin

df = df[:(len(df))-2] #removes trailing incomplete cardiac cycle

groupedRR = df['RRI'].groupby(df['bin'])
RRarray = groupedRR.mean() 

groupedSBP = df['brachial_SBP'].groupby(df['bin'])
SBParray = np.asarray(groupedSBP.mean())
SBParray2 = groupedSBP.mean()

bin_weight = groupedSBP.size()/df.index.max()
df2 = df.mean()

# Linear Regression and time domain metrics

slope, intercept, r_value, p_value, std_err = linregress(SBParray, RRarray)
df2['BRS slope'] = slope
df2['R^2'] = r_value**2
df2['Mean SBP'] = df['brachial_SBP'].mean()
df2['STD SBP'] = np.std(df['brachial_SBP'])
df2['Mean RRI'] = df['RRI'].mean()
df2['Mean HR (bpm)'] = 60000/df['RRI'].mean()
df2['STD HR (bpm)'] = np.std(df['RRI'])
df2['HRV RMSSD'] = np.sqrt(np.mean(np.square(np.diff(df['RRI']))))
df2['NNxx'] = np.sum(np.abs(np.diff(df['RRI'])) > 50)*1
df2['P value'] = p_value



# plots plots plots

sns.regplot(x = SBParray, y = RRarray)
plt.title('BRS')
plt.xlabel('Brachial SBP (mmHg)')
plt.ylabel('RRI')
plt.yticks()
plt.xticks()
plt.xlim(SBParray.min()-1,SBParray.max()+1)
plt.ylim(RRarray.min()-50,RRarray.max()+50)
plt.text(SBParray.min()+1 ,RRarray.max()+35, r'$Slope: %s,\ R^2=%s$' % (np.round(slope, decimals =2), np.round(df2['R^2'], decimals = 2)))

df2 = df2.drop(['bin'])




