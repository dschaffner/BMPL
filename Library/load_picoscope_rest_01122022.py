# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 10:26:52 2022

@author: Josh0
"""

import numpy as np
import scipy.io as spio
import scipy.integrate as sp
import matplotlib.pylab as plt
from scipy import signal
import sys


def load_picoscope_rest_01122022(shot_number, maxrange=1, scopenum=5, time_range=[-6.0, 194.0], location='', plot=False):
    
    def butter_highpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff/nyq
        b, a = signal.butter(order, normal_cutoff, btype='highpass', analog=False)
        return b, a
    
    def butter_highpass_filter(data, cutoff, fs, order=5):
        b, a = butter_highpass(cutoff, fs, order=order)
        y = signal.filtfilt(b, a, data)
        return y
    
    
    if (type(scopenum) == int):
        if scopenum == 2:
           scopename = 'Pico2\\'
        elif scopenum == 3:
           scopename = 'Pico3\\'
        elif scopenum == 4:
           scopename = 'Pico4\\'
        else:
           scopename = 'Pico5\\'
    else:
        print(f'scopenum is not an int, {scopenum}')
        sys.exit()
        
    probe_dia = 0.003175    #m (1/8'' probe)
    probe_dia = 0.00158755  #m (1/16'' probe)
    ##hole_sep = 0.001016     #m (1/16''probe)  ## Aparently unused variable
    r_probe_area = np.pi*(probe_dia/2)**2
    #tz_probe_area = probe_dia*hole_sep  ## Aparently unused variable
    startintg_index = 0 #3000
    meancutoff = 1000
    location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\01122022\\'
    filename = '20220112-0001 ('
    print(location + scopename + filename + str(shot_number)+').mat')
    data = spio.loadmat(location + scopename + filename + str(shot_number) + ').mat', appendmat=True)
    dataraw = data
    
    Bdot3raw = dataraw['A']
    Bdot4raw = dataraw['B']
    Bdot5raw = dataraw['C']
    Bdot6raw = dataraw['D']
    
    start_time = dataraw['Tstart']
    time_interval = dataraw['Tinterval']
    iteration_array = np.array(range(1, 25001, 1))
    time_interval_array = time_interval * iteration_array
    time_s = start_time + time_interval_array
    time_us = time_s*1e6
    timeB_s = time_s[:,1:]
    timeB_us = time_us[:,1:]
    
    Bdot3 = data['A'] - np.mean(data['A'][0:meancutoff])
    neginfs = np.isneginf(Bdot3)
    Bdot3[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot3)
    Bdot3[np.where(posinfs)] = maxrange
    
    Bdot4 = data['B'] - np.mean(data['B'][0:meancutoff])
    neginfs = np.isneginf(Bdot4)
    Bdot4[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot4)
    Bdot4[np.where(posinfs)] = maxrange
    
    Bdot5 = data['C'] - np.mean(data['C'][0:meancutoff])
    neginfs = np.isneginf(Bdot5)
    Bdot5[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot5)
    Bdot5[np.where(posinfs)] = maxrange
    
    Bdot6 = data['D'] - np.mean(data['D'][0:meancutoff])
    neginfs = np.isneginf(Bdot6)
    Bdot6[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot6)
    Bdot6[np.where(posinfs)] = maxrange
    
    B3_probe = np.array(Bdot3/r_probe_area).T
    B4_probe = np.array(Bdot4/r_probe_area).T
    B5_probe = np.array(Bdot5/r_probe_area).T
    B6_probe = np.array(Bdot6/r_probe_area).T
    
    B3 = sp.cumtrapz(B3_probe, time_us)*1e4 #Gauss
    B4 = sp.cumtrapz(B4_probe, time_us)*1e4 #Gauss
    B5 = sp.cumtrapz(B5_probe, time_us)*1e4 #Gauss
    B6 = sp.cumtrapz(B6_probe, time_us)*1e4 #Gauss
    
    B3filt = butter_highpass_filter(B3, 5e4, 125e6, order=3)
    B4filt = butter_highpass_filter(B4, 5e4, 125e6, order=3)
    B5filt = butter_highpass_filter(B5, 5e4, 125e6, order=3)
    B6filt = butter_highpass_filter(B6, 5e4, 125e6, order=3)
    
    if plot:
        fig, ax1 = plt.subplots(1)
        ax1.plot(timeB_us[0], B3filt[0])
        ax1.plot(timeB_us[0], B4filt[0])
        ax1.plot(timeB_us[0], B5filt[0])
        ax1.plot(timeB_us[0], B6filt[0])
    
    return time_s, time_us, timeB_s, timeB_us, Bdot3raw, Bdot4raw, Bdot5raw, Bdot6raw, Bdot3, Bdot4, Bdot5, Bdot6, B3, B4, B5, B6, B3filt, B4filt, B5filt, B6filt
