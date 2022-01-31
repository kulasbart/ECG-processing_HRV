# ECG Signal Processing & HRV analysis

Guide to the basics of working with ECG signals, R-R peak detection, fitting functions

Also provided are metrics for time-series ECG data, including heart rate variability periodogram analysis using welch's method 

# BRS-regression 

Linear model plotting systolic blood pressure (SBP) and r-r intervals (RRI) by binning ALL beats in a given sample. Includes the following time and frequency domain metrics for HRV:

STD RR/SDNN (ms)
STD HR (beats/min)
RMSSD (ms)
pNNxx (%)

Power of the following frequency bands:
Very Low Frequency (VLF): 0 - 0.4 Hz
Low Frequency (LF): 0.04 - 0.15 HZ
High Frequency (HF): 0.15 - 0.4 Hz
Total Power (VLF + LF + HF)
LF/HF ratio
Peak (Hz) in VLF, LF, HF bands
