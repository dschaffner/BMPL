# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 22:53:59 2021

@author: Josh0
"""

# Inner coil flux comparisons at Mid, GunEdge, Outside Bdot probe locations with COMSOL simulations of the same locations.

import numpy as np
import scipy.integrate as sp
import matplotlib.pylab as plt


location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2020\\Flux_10142020\\'
scopename = '10142020pico1\\'
filename1 = '20201014-0001'
data1 = np.loadtxt(location+scopename+filename1+'.txt', skiprows=3, unpack=True)

comsol_location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\COMSOL_SimulationData\\Data\\08182021\\'

comsol_gunedge = '08182021_COM_GE_norm_0p625cm'
comsol_gunedge1cm = '08182021_COM_GE_norm_1cm'
comsol_gunedge1p25cm = '08182021_COM_GE_norm_1p25cm'

comsol_middle = '08182021_COM_Middle_zcomp_0p625cm'
comsol_middle1cm = '08182021_COM_Middle_zcomp_1cm'
comsol_middle1p25cm = '08182021_COM_Middle_zcomp_1p25cm'
#comsol_currentPhi1cm = '09222021_COM_Mid_Current_densityPhi'

comsol_outside = '08182021_COM_Outside_zrnorm_0p625cm'
comsol_outside1cm = '08182021_COM_Outside_norm_1cm'
comsol_outside1p25cm = '08182021_COM_Outside_norm_1p25cm'


com_fluxGE = np.loadtxt(comsol_location+comsol_gunedge+'.txt', skiprows=5, unpack=True)
com_fluxGE1cm = np.loadtxt(comsol_location+comsol_gunedge1cm+'.txt', skiprows=5, unpack=True)
com_fluxGE1p25cm = np.loadtxt(comsol_location+comsol_gunedge1p25cm+'.txt', skiprows=5, unpack=True)

com_fluxMid = np.loadtxt(comsol_location+comsol_middle+'.txt', skiprows=5, unpack=True)
com_fluxMid1cm = np.loadtxt(comsol_location+comsol_middle1cm+'.txt', skiprows=5, unpack=True)
com_fluxMid1p25cm = np.loadtxt(comsol_location+comsol_middle1p25cm+'.txt', skiprows=5, unpack=True)
#com_currentMid1cm = np.loadtxt(comsol_location+comsol_currentPhi1cm+'.txt', skiprows=5, unpack=True)

com_fluxOut = np.loadtxt(comsol_location+comsol_outside+'.txt', skiprows=5, unpack=True)
com_fluxOut1cm = np.loadtxt(comsol_location+comsol_outside1cm+'.txt', skiprows=5, unpack=True)
com_fluxOut1p25cm = np.loadtxt(comsol_location+comsol_outside1p25cm+'.txt', skiprows=5, unpack=True)

plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_tallplots.py')
#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')

# data[0] = time(ms), data[1] = emf_outside, data[2] = emf_middle, data[3] = emf_gunEdge, data[4] = discharge current
time_ms1 = data1[0]
time_s1 = time_ms1*1e-3
emf_outside1 = data1[1]
emf_middle1 = data1[2]
emf_gunEdge1 = data1[3]


probe_diameter = 0.09652 #m (3.8"probe)
probe_area = np.pi*((probe_diameter/2)**2)

comsol_time = com_fluxGE[0]
comsol_fluxGE = com_fluxGE[1]
comsol_fluxGE1cm = com_fluxGE1cm[1]
comsol_fluxGE1p25cm = com_fluxGE1p25cm[1]

comsol_fluxMid = com_fluxMid[1]
comsol_fluxMid1cm = com_fluxMid1cm[1]
comsol_fluxMid1p25cm = com_fluxMid1p25cm[1]
#comsol_currentMid1cm = com_currentMid1cm[1]

comsol_fluxOut = com_fluxOut[3]
comsol_fluxOut1cm = com_fluxOut1cm[1]
comsol_fluxOut1p25cm = com_fluxOut1p25cm[1]

flux_middle1 = sp.cumtrapz(emf_middle1, time_ms1)
flux_gunEdge1 = sp.cumtrapz(emf_gunEdge1, time_ms1)
flux_outside1 = sp.cumtrapz(emf_outside1, time_ms1)

#fig1, (ax1) = plt.subplots(1)
fig1, (ax1, ax2, ax3) = plt.subplots(3)
#fig2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True)

#Flux Middle 0.625cm Inner Electrode
ax1.plot(time_ms1[1:210716]+1.4, flux_middle1[:210715], color='red', label='Inner Coil - Middle')
ax1.plot(comsol_time, comsol_fluxMid1cm, color='black', label='1cm COM Middle')
ax1.set_ylabel(r'Magnetic Flux (mWb)')
ax1.set_xlabel(r'Time (ms)')
ax1.legend(loc='best', frameon=False)
#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\08182021\\08182021_1p25cm_middle.pdf', dpi=600)

ax2.plot(time_ms1[1:210716]+1.4, flux_gunEdge1[:210715]*-1, color='red', label='GunEdge')
ax2.plot(comsol_time, comsol_fluxGE1cm, color='black', label='1cm COM GunEdge')
ax2.set_ylabel(r'Magnetic Flux (mWb)')
ax2.set_xlabel(r'Time (ms)')
ax2.legend(loc='best', frameon=False)

ax3.plot(time_ms1[1:210716]+1.4, flux_outside1[:210715]*-1, color='red', label='Outside')
ax3.plot(comsol_time, comsol_fluxOut1cm, color='black', label='1cm COM Outside')
ax3.set_xlabel(r'Delay Time (ms)')
ax3.set_ylabel(r'Flux (mWb)')
ax3.legend(loc='best', frameon=False)

#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\08182021\\08182021_1p25cm_comparison.pdf', dpi=600)






