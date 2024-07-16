# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:41:33 2021

@author: Josh0
"""

# This code is written to compare physical Far/Near Outer Coil data with LTSpice Simulations (Confirming Manufacturer Inductances)

import numpy as np
import matplotlib.pylab as plt
import pandas as pd

plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')

LTSpice_filepath = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\LTSpice Simulations\\'
LTSpice_filename = '08052021_Far_Outer_LTSpice_4p27mH_80mF'
LTSpice_filename_6mH = '08052021_Far_Outer_LTSpice_6mH_80mF'
LTSpice_filename_7mH_85mF = '08052021_Far_Outer_LTSpice_85mF_7mH'
LTSpice_filename_7mH_100mF = '08052021_Far_Outer_LTSpice_7mH_100mF'
LTSpice_filename_6p5mH_90mF = '08052021_Far_Outer_LTSpice_6p5mH_90mF'
LTSpice_filename_6mH_85mF = '08052021_Far_Outer_LTSpice_6mH_85mF'
LTSpice_filename_6mH_87p5mF = '08052021_Far_Outer_LTSpice_6mH_87p5mF'
LTSpice_filename_5p5mH_83p5mF = '08052021_Far_Outer_LTSpice_5p5mH_83p5mF'
LTSpice_filename_5p5mH_85mF = '08052021_Far_Outer_LTSpice_5p5mH_85mF'
LTSpice_filename_5p5mH_90mF = '08052021_Far_Outer_LTSpice_5p5mH_90mF'
LTSpice_filename_5p5mH_88mF = '08052021_Far_Outer_LTSpice_5p5mH_88mF'

LTspice_data1 = np.loadtxt(LTSpice_filepath+LTSpice_filename+'.txt', skiprows=1, unpack=True)
LTspice_data2 = np.loadtxt(LTSpice_filepath+LTSpice_filename_5p5mH_88mF+'.txt', skiprows=1, unpack=True)

LTspice_time1 = LTspice_data1[0]
LTspice_current1 = LTspice_data1[1]
LTspice_time2 = LTspice_data2[0]
LTspice_current2 = LTspice_data2[1]


Oscilloscope_data = pd.read_csv("C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Oscilloscope Traces\\08052021_Far_Outer_200VDischarge_500ms.csv", usecols=(3, 4))
oscilloscope_time = Oscilloscope_data['ms'].values
oscilloscope_time_shift = oscilloscope_time + 0.01

oscilloscope_voltage = Oscilloscope_data['V'].values
oscilloscope_current = oscilloscope_voltage*1e3


fig, (ax1) = plt.subplots(1, sharex=True, sharey=True)

#ax1.plot(LTspice_time1, LTspice_current1, color='blue', label='LTSpice - L=4.27mH, C=80mF (Manufacturer Inductance)')
ax1.plot(LTspice_time2, LTspice_current2, color='blue', label='LTSpice L=5.5mH, C=88mF')
ax1.plot(oscilloscope_time, oscilloscope_current, color='red', label='Oscilloscope Trace', alpha=0.4)
ax1.set_ylabel(r'Current (A)')
ax1.set_xlabel(r'Time (s)')
ax1.legend(loc='best', frameon=False)

plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\Spice_OscillComparison_Far200V_5p5mH_88mF.pdf', dpi=600)


