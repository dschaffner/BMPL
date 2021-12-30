# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 13:03:28 2021

@author: Josh0
"""

# Takes the 08062021 oscilloscope current(V) vs time data and writes to .txt file

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')


FarOuter = pd.read_csv('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Oscilloscope Traces\\08052021_Far_Outer_200VDischarge_500ms.csv', usecols=(3,4))
time = FarOuter['ms'].values
voltage = FarOuter['V'].values

time_ms = time*1e3
current = voltage*1e3

smoothed_current = savgol_filter(current, 11, 1)

fig, (ax1) = plt.subplots(1)
ax1.plot(time_ms, smoothed_current, color='blue', label='200V Far Discharge')
ax1.legend(loc='best', frameon=False)
ax1.set_ylabel('Current (A)')
ax1.set_xlabel('Time (ms)')

datatobewrittenTxt = np.column_stack([time_ms, smoothed_current])
Far200V_CurrentandTime_path = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\Far200V_CurrentandTime.txt'
np.savetxt(Far200V_CurrentandTime_path, datatobewrittenTxt, fmt=['%10.8f', '%10.8f'])