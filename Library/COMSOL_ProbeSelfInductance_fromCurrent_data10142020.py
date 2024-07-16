# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:17:25 2021

@author: Josh0
"""

# Calculating the self inductance of the probe due to inner coil magnetic field (focusing on middle probe), and subtracting this from COMSOL's flux arrays and comparing to experiment.

import numpy as np
import scipy.integrate as sp
import matplotlib.pylab as plt


location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2020\\Flux_10142020\\'
scopename = '10142020pico1\\'
filename1 = '20201014-0001'
data1 = np.loadtxt(location+scopename+filename1+'.txt', skiprows=3, unpack=True)

comsol_location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\COMSOL_SimulationData\\Data\\11012021\\'

comsol_middle1cm = '11012021_com_middle1cm_0to20ms'
comsol_currentPhi_Mid1cm = '11012021_com_CurrentMiddle_phi_0to20ms'

comsol_gunedge1cm = '11012021_com_GE_1cm_norm_r_z_0to20ms'
comsol_currentPhi_GE1cm = '11012021_com_CurrentGE_phi_0to20ms'


com_fluxMid1cm = np.loadtxt(comsol_location+comsol_middle1cm+'.txt', skiprows=5, unpack=True)
com_currentMid1cm = np.loadtxt(comsol_location+comsol_currentPhi_Mid1cm+'.txt', skiprows=5, unpack=True)
com_fluxGE1cm = np.loadtxt(comsol_location+comsol_gunedge1cm+'.txt', skiprows=5, unpack=True)
com_currentGE1cm = np.loadtxt(comsol_location+comsol_currentPhi_GE1cm+'.txt', skiprows=5, unpack=True)

#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_presentationplots.py')

# data[0] = time(ms), data[1] = emf_outside, data[2] = emf_middle, data[3] = emf_gunEdge, data[4] = discharge current
time_ms1 = data1[0]
time_s1 = time_ms1*1e-3
emf_outside1 = data1[1]
emf_middle1 = data1[2]
emf_gunEdge1 = data1[3]*-1

probe_diameter = 0.09652 #m (3.8"probe)
probe_area = np.pi*((probe_diameter/2)**2)

comsol_time = com_currentMid1cm[0]

comsol_fluxMid1cm = com_fluxMid1cm[1]
comsol_currentMid1cm = com_currentMid1cm[1]
comsol_dIdt_1cm = np.diff(comsol_currentMid1cm)
comsol_self_emf = comsol_dIdt_1cm*0.000000696 #Inductance of probe in henries (calculated manually)


flux_Middle_SelfEMF = sp.trapz(comsol_self_emf)

comsol_fluxGE1cm = com_fluxGE1cm[1]
comsol_currentGE1cm = com_currentGE1cm[1]
comsol_GE_dIdt_1cm = np.diff(comsol_currentGE1cm)
comsol_GE_self_emf = comsol_GE_dIdt_1cm*0.000000696

flux_GE_SelfEMF = sp.trapz(comsol_GE_self_emf)

#comsol_fluxOut1cm = com_fluxOut1cm[1]
#EXPERIMENTAL DATA:
flux_middle1 = sp.cumtrapz(emf_middle1, time_ms1)
flux_gunEdge1 = sp.cumtrapz(emf_gunEdge1, time_ms1)
flux_outside1 = sp.cumtrapz(emf_outside1, time_ms1)
element_array = np.array(range(44645, 210715, 8929))
shortened_flux_middle = np.array(flux_middle1[element_array])
shortened_flux_GE = np.array(flux_gunEdge1[element_array])


COM_flux_mid_with_selfInductance = comsol_fluxMid1cm[1:20] + flux_Middle_SelfEMF
COM_flux_GE_with_selfInductance = comsol_fluxGE1cm[1:20] - flux_GE_SelfEMF

fig1, (ax1, ax2) = plt.subplots(2)
#fig1, (ax1, ax2, ax3) = plt.subplots(3)
#fig2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True)


ax1.plot(time_ms1[1:210716]+1.4, flux_middle1[:210715], color='blue', label='Flux Middle')
ax1.plot(comsol_time[:19]+0.75, COM_flux_mid_with_selfInductance, color='black', label='Flux from COMSOL current')
ax1.plot(comsol_time, comsol_fluxMid1cm, color='red', label='Flux_NoSelfEMF', alpha=0.5)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Flux (mWb)')
ax1.legend(loc='best', frameon=False)

# GunEdge flux doesn't change as much from the self inductance compared with the Flux Middle.  Why is this??
ax2.plot(time_ms1[1:210716]+1.4, flux_gunEdge1[:210715], color='blue', label='Flux GunEdge')
ax2.plot(comsol_time[:19]+0.75, COM_flux_GE_with_selfInductance, color='black', label='Flux from COM. current')
ax2.plot(comsol_time, comsol_fluxGE1cm, color='red', label='Flux_NoSelfEMF')
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Flux (mWb)')
ax2.legend(loc='best', frameon=False)

#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\10222021_flux_comparisons_withSelfInductance.pdf', dpi=600)
#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\11012021_flux_comparisonMiddle_withSelfInductance_Poster.png', dpi=600)






