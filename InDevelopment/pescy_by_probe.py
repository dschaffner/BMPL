# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
"""
import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import os
from load_hdf5 import load_hdf5
import spectrum_wwind as spec
import indexfinderfuncs as iff
import get_corr as gc
import scipy.integrate as sp
from scipy.interpolate import interp1d

from calc_PE_SC import PE, CH, PE_dist, PE_calc_only
from collections import Counter
from math import factorial
#%%cell 1


from scipy.signal import butter, sosfiltfilt, sosfreqz

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos
#%%
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y

import smooth as sm

def iter_smooth(array,loops=6,window_len=3):
    for l in np.arange(loops):
        array = sm.smooth(array,window_len=window_len)
    return array

def cross_spectrum(array1,array2,interval):
    n = array1.shape[0]
    print ('n = ',n)
    factor = 2.0/(interval)
    cross_spec = np.conj(array1)*(array2)*factor
    cross_phase = np.angle(cross_spec)
    arr1_autospec = np.conj(array1)*(array1)*factor
    arr1_autospec = np.real(arr1_autospec)
    arr2_autospec = np.conj(array2)*(array2)*factor
    arr2_autospec = np.real(arr2_autospec)
    cross_coh = ((np.abs(cross_spec))**2)/(arr1_autospec*arr2_autospec)

    return cross_spec,cross_phase,arr1_autospec,arr2_autospec,cross_coh

# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
datafilename='Dataset_01122022.h5'
data=load_hdf5(directory+datafilename,verbose=True)


time_s = data['time']['time_s']
timeB_s = time_s[1:]
time_us = data['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m
dt=timeB_s[1]-timeB_s[0]

numshots=94
direction_list=['r','t','z']
probelist=['probe5','probe7','probe19','probe21','probe33','probe35']
directions = len(direction_list)
numprobes = 1#len(probelist)

#determine FFT size and generate an output array
#fsize=int((data['mag_probe']['positions']['probe5']['r']['bdot'][0,start_time_index:end_time_index].shape[0]))#/2)+1
#avebspec_frombdot = np.zeros([numprobes,directions,fsize])
#avebspec_direct = np.zeros([numprobes,directions,fsize])
#avebmagspec = np.zeros([numprobes,fsize])
#spec1 = np.zeros([numshots,numprobes,directions,fsize])
##spec2 = np.zeros([numshots,numprobes,directions,fsize])
#crossphase_fromb = np.zeros([numshots,numprobes,directions,fsize])
#kfrom_crossphase = np.zeros([numshots,numprobes,directions,fsize])

delay_array = np.arange(1,5000)
num_delays = len(delay_array)
embeddelay = 5
nfac = factorial(embeddelay)

PEs = np.zeros([num_delays])
SCs = np.zeros([num_delays])

for loop_delay in np.arange(len(delay_array)):
    if (loop_delay%100)==0: print( 'On Delay ',delay_array[loop_delay])
    permstore_counter = []
    permstore_counter = Counter(permstore_counter)
    tot_perms = 0
    for shot in np.arange(numshots):
        print ('On Shot: ',shot)
        for probe_index, probe in enumerate(probelist):
            #print ('On Probe: ',probe)
            for direction_index, direction in enumerate(direction_list):
                #print('On Direction: ',direction)
                data1=data['mag_probe']['positions'][probe][direction]['b'][shot,:]
                arr,nperms = PE_dist(data1,embeddelay,delay=delay_array[loop_delay])
                permstore_counter = permstore_counter+arr
                tot_perms = tot_perms+nperms
    PE_tot,PE_tot_Se = PE_calc_only(permstore_counter,tot_perms)
    C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*np.log2(nfac))
             /((1 + 1./nfac)*np.log2(nfac+1) - 2*np.log2(2*nfac) 
            + np.log2(nfac))*(PE_tot/np.log2(nfac)))
    PEs[loop_delay]=PE_tot/np.log2(nfac)
    SCs[loop_delay]=C



