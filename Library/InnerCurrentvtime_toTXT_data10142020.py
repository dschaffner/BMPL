# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 13:08:52 2021

@author: Josh0
"""

#This is a code to get a current .txt file for inner coil flux measurements

import numpy as np
import matplotlib.pylab as plt
from scipy.signal import savgol_filter


location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2020\\Flux_10142020\\'
scopename = '10142020pico1\\'

filename1 = '20201014-0001'
data1 = np.loadtxt(location+scopename+filename1+'.txt',skiprows=3,unpack=True)

#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_presentationplots.py')


"""
Data structure follows:
    0           1               2               3               4
    time(ms)    emf_outside     emf_middle      emf_gunEdge     coil_current(V) - conversion below


As of now, the coil current is in volts, but we need to be in amps.

Working backwards, 333mV : 20A, then 5A : 3000A

So for the first data point of 0.00195319V = 1.95319mV that corresponds to:
    1.95319mV * (20/333) * (3000/5) = 1.95319mV * 0.06006A * 600 = 70.378308A
    
    Conversion factor is 0.06006 * 600 = 36.036
"""

time_ms1 = data1[0]
time_s = time_ms1*1e-3

#Shift start of discharge to 0sec
time_ms1 = time_ms1+1.5281281

coil_current_in_V1 = data1[4]
coil_current_in_mV1 = coil_current_in_V1*1e3
coil_current1 = coil_current_in_mV1*36.036

smoothed_coil_current1 = savgol_filter(coil_current1, 1051, 1)
smoothed_coil_current_kA1 = smoothed_coil_current1*1e-3
#smoothed_coil_current4 = savgol_filter(coil_current4, 551, 1)
#smoothed_coil_current_kA4 = smoothed_coil_current4*1e-3

time_ms_window1 = time_ms1[31000:]
smoothed_coil_current_window1 = smoothed_coil_current_kA1[31000:]

fig, (ax1) = plt.subplots(1)
ax1.plot(time_ms1, smoothed_coil_current_kA1)
#ax1.plot(time_ms1, smoothed_coil_current1, color='red')
#ax1.axhline(0, 0, 1)
#ax1.set_title('Smoothed Inner Coil Current - Picoscope Capture')
ax1.set_ylabel('Current (kA)')
ax1.set_xlabel('Time (ms)')

#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\smoothed_innerCoil_current.pdf', dpi=600)


data_to_be_writtenTxt = np.column_stack([time_ms_window1, smoothed_coil_current_window1])
CurrentandTime_0to45_path = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\CurrentandTime_0to45.txt'
np.savetxt(CurrentandTime_0to45_path, data_to_be_writtenTxt, fmt=['%10.8f','%10.8f'])
