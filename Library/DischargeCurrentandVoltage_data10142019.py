# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 12:47:44 2021

@author: Josh0
"""

import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import numpy as np
import indexfinderfuncs as iff
from scipy.signal import savgol_filter

"""data:
    discharge:
        dis_I
        dis_I_raw
        dis_V
        dis_V_raw
        
    time:
        time_s
        time_us
"""
plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')
#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_fullpageplots.py')

# Load the HDF5 file:
data_directory_location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\'
datafilename = '2kV_0p5ms_10msGas_20shots_DisHV_10142019.h5'
data = load_hdf5(data_directory_location + datafilename, verbose = True)

# Calling time from our data
time_us = data['time']['time_us'][:]
time_s = data['time']['time_s'][:]

# Calling the discharge current and discharge voltages:
dis_current = data['discharge']['dis_I'][:]
dis_current = dis_current*1e-3
raw_dis_current = data['discharge']['dis_I_raw'][:]
dis_voltage = data['discharge']['dis_V'][:]/1e3
raw_dis_voltage = data['discharge']['dis_V_raw'][:]

start_time = 50
end_time = 150
time_range = [start_time, end_time]
start_time_index = iff.tindex_min(start_time, time_us)
end_time_index = iff.tindex_min(end_time, time_us)

# Shot 4 was a dud, so its excluded from the data calculation.
avg_dis_current = np.zeros(25004)
for s in range (0, 21):
    if s == 3:
        continue
    avg_dis_current = avg_dis_current + np.array((dis_current[s,:])/20)

avg_dis_voltage = np.zeros(25004)
for s in range(0, 21):
    if s == 3:
        continue
    avg_dis_voltage = avg_dis_voltage + np.array((dis_voltage[s,:])/20)
    

max_I = np.zeros(25004)
min_I = np.zeros(25004)
# Array excluding the problematic shot:
x1 = np.array(range(0, 3))
x2 = np.array(range(4, 21))
x = np.append(x1, x2)

for i in range(25004):
    max_I[i] = max(dis_current[x, i])
    min_I[i] = min(dis_current[x, i])

max_V = np.zeros(25004)
min_V = np.zeros(25004)
for i in range(25004):
    max_V[i] = max(dis_voltage[x, i])
    min_V[i] = min(dis_voltage[x, i])

zero_array_for_AvgLines = np.zeros(end_time_index - start_time_index)
mean_discurrent_array = np.zeros(21)
mean_disvoltage_array = np.zeros(21)

for s in range (0, 21):
     mean_discurrent_array[s] = np.mean(dis_current[s,start_time_index:end_time_index])
     mean_disvoltage_array[s] = np.mean(dis_voltage[s, start_time_index:end_time_index])
     
# Now exclude shot 4:
     ## Is there an easier or fancier way to do this not brute force?
dis_current_append1 = mean_discurrent_array[0:3]
dis_current_append2 = mean_discurrent_array[4:]
mean_discurrent_array = np.append(dis_current_append1, dis_current_append2)

dis_voltage_append1 = mean_disvoltage_array[0:3]
dis_voltage_append2 = mean_disvoltage_array[4:]
mean_disvoltage_array = np.append(dis_voltage_append1, dis_voltage_append2)

smoothed_dis_current = savgol_filter(avg_dis_current, 11, 3)
smoothed_max_current = savgol_filter(max_I, 15, 3)
smoothed_min_current = savgol_filter(min_I, 15, 3)

current_line = zero_array_for_AvgLines + np.mean(mean_discurrent_array)
voltage_line = zero_array_for_AvgLines + np.mean(mean_disvoltage_array)

# Here is the masking process: omitting any values less than 1
threshold_I = 1
min_I = np.ma.array(smoothed_min_current)
max_I = np.ma.array(smoothed_max_current)
avg_I = np.ma.array(smoothed_dis_current)
max_V = np.ma.array(max_V)
min_V = np.ma.array(min_V)
avg_V = np.ma.array(avg_dis_voltage)

masked_minI = np.ma.masked_where(min_I < threshold_I, min_I)
masked_maxI = np.ma.masked_less(max_I, 0)
masked_avgI = np.ma.masked_less(avg_I, 0)
masked_maxV = np.ma.masked_outside(max_V, 0, 2)
masked_minV = np.ma.masked_outside(min_V, 0, 2)
masked_avgV = np.ma.masked_outside(avg_V, 0, 2)

fig1, (ax1, ax2) = plt.subplots(2)
ax1.plot(time_us, masked_maxI, linewidth=0.5, color='blue', alpha=0.4)
ax1.plot(time_us, masked_minI, linewidth=0.5, color='blue', alpha=0.4)
ax1.plot(time_us, masked_avgI, linewidth=1, color='black')
ax1.plot(time_us[start_time_index:end_time_index], current_line, linewidth=2, linestyle='--', alpha=1, color='red', label= r'Mean Current = 47.4kA')
ax1.fill_between(time_us, masked_minI, masked_maxI, alpha=0.3, color='blue')
ax1.set_ylabel(r'Discharge Current (kA)')
ax1.set_xlabel(r'Time ($\mu$s)')
ax1.text(0.07, 0.92,'(a)',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
ax1.legend(loc='best', frameon=False)

ax2.plot(time_us, masked_maxV, linewidth=0.5, color='blue', alpha=0.4)
ax2.plot(time_us, masked_minV, linewidth=0.5, color='blue', alpha=0.4)
ax2.plot(time_us, masked_avgV, linewidth=1, color='black')
ax2.plot(time_us[start_time_index:end_time_index], voltage_line, linewidth=2, linestyle='--', alpha=1, color='red', label=r'Mean Voltage = 0.485kV')
ax2.fill_between(time_us, masked_minV, masked_maxV, alpha=0.3, color='blue')
ax2.set_ylabel(r'Discharge Voltage (kV)')
ax2.set_xlabel(r'Time ($\mu$s)')
ax2.text(0.07, 0.92,'(b)',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
ax2.legend(loc='best', frameon=False)















