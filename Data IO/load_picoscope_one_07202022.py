# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:11:59 2022

@author: Josh0
"""


import numpy as np
import scipy.io as spio
import scipy.integrate as sp
import matplotlib.pylab as plt
from scipy import signal
import sys


###################################################
""" The picoscope-probe structure 7202022 """
###################################################
"""          A           B           C           D
pico1       HV          DC          5R          5T
pico2       5Z          7R          7T          7Z
pico3       19R         19T         19Z         21R
pico4       21T         21Z         33R         33T
pico5       33Z         35R         35T         35Z

N --- Density
DC --- Discharge Current
HV --- High Voltage
Total Shots = 61, 51 good
Nozzle Strength = 225V
"""
###################################################


def load_picoscope_one_07202022(shot_number, maxrange=1, scopenum=1, time_range=[-6.0, 194.0], location='', plot=False):
    
    def butter_highpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff/nyq
        b, a = signal.butter(order, normal_cutoff, btype='highpass', analog=False)
        return b, a
    
    def butter_highpass_filter(data, cutoff, fs, order=5):
        b, a = butter_highpass(cutoff, fs, order=order)
        y = signal.filtfilt(b, a, data)
        return y
    
  #if (type(scopenum) == int):
#        if scopenum == 1:
#            scopename = 'Pico1\\'
#        elif scopenum == 2:
#           scopename = 'Pico2\\'
#       elif scopenum == 3:
#           scopename = 'Pico3\\'
#       elif scopenum == 4:
#           scopename = 'Pico4\\'
#       else:
#           scopename = 'Pico5\\'
#    else:
#        print(f'scopenum is not an int, {scopenum}')
#        sys.exit()
    
    
    probe_dia = 0.003175    #m (1/8'' probe)
    probe_dia = 0.00158755  #m (1/16'' probe)
    r_probe_area = np.pi*(probe_dia/2)**2
    startintg_index = 0 #3000
    meancutoff = 1000

    location = 'C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\07202022\\Pico1\\'
    filename = '20220720-0001 ('
    print(location + filename + str(shot_number) +').mat')
    data = spio.loadmat(location+filename+str(shot_number)+').mat', appendmat=True)
    dataraw = data
    
    HighV_raw = dataraw['A']
    DisCurrent_raw = dataraw['B']
    Bdot1raw = dataraw['C']
    Bdot2raw = dataraw['D']
    
    DisCurrent_raw[np.where(DisCurrent_raw==1)]=5.0
    DisCurrent_raw[np.where(DisCurrent_raw==-1)]=-5.0
    HighV_raw[np.where(HighV_raw==1)]=5.0
    HighV_raw[np.where(HighV_raw==-1)]=-5.0
    
    Rogowski_gain = 10000.0
    Rogowski_dir = 1.0
    Rogowski_factor = 2.0
    HV_gain = 1000.0
    HV_dir = -1.0
    
    neginfs = np.isneginf(DisCurrent_raw)
    DisCurrent_raw[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(DisCurrent_raw)
    DisCurrent_raw[np.where(posinfs)] = maxrange
    DisCurrent = DisCurrent_raw - np.mean(DisCurrent_raw[0:meancutoff])
    DisCurrent = DisCurrent * Rogowski_gain * Rogowski_dir * Rogowski_factor
    DisCurrent = DisCurrent.T
    
    neginfs = np.isneginf(HighV_raw)
    HighV_raw[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(HighV_raw)
    HighV_raw[np.where(posinfs)] = maxrange
    HighV = HighV_raw - np.mean(HighV_raw[0:meancutoff])
    HighV = HighV * HV_gain * HV_dir
    HighV = HighV.T
    
    start_time = dataraw['Tstart']
    time_interval = dataraw['Tinterval']
    iteration_array = np.array(range(1, 25001, 1))
    time_interval_array = time_interval * iteration_array
    time_s = start_time + time_interval_array
    time_us = time_s*1e6
    timeB_s = time_s[:,1:]
    timeB_us = time_us[:,1:]
        
    Bdot1 = data['C']
    neginfs = np.isneginf(Bdot1)
    Bdot1[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot1)
    Bdot1[np.where(posinfs)] = maxrange
    Bdot1 = Bdot1-np.mean(data['C'][0:meancutoff])
    
    Bdot2 = data['D']
    neginfs = np.isneginf(Bdot2)
    Bdot2[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot2)
    Bdot2[np.where(posinfs)] = maxrange
    Bdot2 = Bdot2-np.mean(data['D'][0:meancutoff])

    B1_probe = np.array(Bdot1/r_probe_area).T
    B2_probe = np.array(Bdot2/r_probe_area).T
    B1 = sp.cumtrapz(B1_probe, time_s)*1e4 #Gauss
    B2 = sp.cumtrapz(B2_probe, time_s)*1e4 #Gauss
    
    B1filt = butter_highpass_filter(B1,5e4,125e6,order=3)
    B2filt = butter_highpass_filter(B2,5e4,125e6,order=3)
    
    if plot:
        fig, (ax1, ax2, ax3) = plt.subplots(3)
        ax1.plot(time_us[0], DisCurrent[0], label='Current')
        ax1.legend(loc='best')
        ax2.plot(time_us[0], HighV[0], label='HV')
        ax2.legend(loc='best')
        ax3.plot(timeB_us[0], B1[0], label='5R')
        ax3.plot(timeB_us[0], B2[0], label='5T')
        ax3.legend(loc='best')
    
    return time_s, time_us, timeB_s, timeB_us, DisCurrent_raw, HighV_raw, DisCurrent, HighV, Bdot1raw, Bdot2raw, Bdot1, Bdot2, B1, B2, B1filt, B2filt




