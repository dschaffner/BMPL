# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 23:28:27 2020

@author: Josh0
"""

# Inner coil flux data from 10142020, and injected helicity as a function of flux

import numpy as np
import scipy.integrate as sp
import matplotlib.pylab as plt


location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2020\\Flux_10142020\\'
scopename = '10142020pico1\\'
filename1 = '20201014-0001'
#filename2 = '20201014-0002'
#filename3 = '20201014-0003'
#filename4 = '20201014-0004'
data1 = np.loadtxt(location+scopename+filename1+'.txt',skiprows=3,unpack=True)
#data2 = np.loadtxt(location+scopename+filename2+'.txt',skiprows=3,unpack=True)
#data3 = np.loadtxt(location+scopename+filename3+'.txt',skiprows=3,unpack=True)
#data4 = np.loadtxt(location+scopename+filename4+'.txt',skiprows=3,unpack=True)

#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_presentationplots.py')

# data[0] = time(ms), data[1] = emf_outside, data[2] = emf_middle, data[3] = emf_gunEdge, data[4] = discharge current
time_ms1 = data1[0]
time_s1 = time_ms1*1e-3
time_ms = time_ms1[1:]+1.4
emf_outside1 = data1[1]
emf_middle1 = data1[2]
emf_gunEdge1 = data1[3]

#time_ms2 = data2[0]
#time_s2 = time_ms2*1e-3
#emf_middle2 = data2[2]

#time_ms3 = data3[0]
#time_s3 = time_ms3*1e-3
#emf_middle3 = data3[2]

#time_ms4 = data4[0]
#time_s4 = time_ms4*1e-3
#emf_middle4 = data3[2]

probe_diameter = 0.09652 #m (3.8"probe)
probe_area = np.pi*((probe_diameter/2)**2)


flux_middle1 = sp.cumtrapz(emf_middle1, time_ms1)
flux_gunEdge1 = sp.cumtrapz(emf_gunEdge1, time_ms1)
flux_outside1 = sp.cumtrapz(emf_outside1, time_ms1)

#flux_middle2 = sp.cumtrapz(emf_middle2, time_ms2)
#flux_middle3 = sp.cumtrapz(emf_middle3, time_ms3)
#flux_middle4 = sp.cumtrapz(emf_middle4, time_ms4)

# Remember the bit about locating the index of a single value in a data set via enumerate:
for i, j in enumerate(time_ms1):
    if j >= 20 and j <= 20.0001:
        print(i)
        # This gives us the index where time_ms1 = 20ms.

fig1, (ax1, ax2) = plt.subplots(2)
ax1.plot(time_ms, flux_middle1, color='blue', label='Flux Middle')
ax1.plot(time_ms, flux_gunEdge1*-1, color='red', label='Flux Gun Edge')
ax1.plot(time_ms, flux_outside1*-1, color='green', label='Flux Outside')
ax1.plot(1.5, 2.009, 'ro', markersize=5, color='black', label='2 mWb at 1.5ms delay')
ax1.set_ylabel(r'Flux (mWb)')
ax1.set_xlabel(r'Delay Time (ms)')
ax1.axvline(1.5, ymin=0, ymax=0.36, linestyle='--', color='black', linewidth=0.7)
ax1.axhline(2.009, xmin=0, xmax=0.143, linestyle='--', color='black', linewidth=0.7)
ax1.tick_params(axis='x')
ax1.tick_params(axis='y')
ax1.legend(loc='best', frameon=False)
    
# I don't need magnetic field as well because it only scales the flux.
# Going to look into the injected helicity, defined as integral(flux*voltage)
# Average voltage from 50 to 150us from 10142019 is 484.94 V
# Average voltage from 50 to 200us from 10142019 is 417.766 V

int_helicity_50_150us1 = flux_middle1*484.94
int_helicity_50_200us1 = flux_middle1*417.766 # Units of V*Wb

inj_helicity_50_1501 = sp.cumtrapz(int_helicity_50_150us1, time_ms1[1:]*1e-3)
inj_helicity_50_2001 = sp.cumtrapz(int_helicity_50_200us1, time_ms1[1:]*1e-3)

ax2.plot(time_ms[1:]+1.3, inj_helicity_50_1501, color='red', label='Injected Helicity Rate')
ax2.plot(1.5, 1.1, 'ro', markersize=5, color='black', label='1.216 Wb$\cdot$m$^{-2}$ at 1.5ms delay')
ax2.set_ylabel(r'Inj. Helicity Rate (Wb$\cdot$m$^{-2}$)')
ax2.set_xlabel(r'Delay Time (ms)')
ax2.tick_params(axis='x')
ax2.tick_params(axis='y')
ax2.legend(loc='best', frameon=False)



#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\9222021_redo_flux_middleCoil.eps', dpi=600)


