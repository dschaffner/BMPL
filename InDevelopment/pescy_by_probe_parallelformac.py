# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
"""
import time
start=time.time()
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


from PESCy_functions import *
from collections import Counter
from math import factorial


# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
directory='/Users/dschaffner/Dropbox/Data/BMPL/BMX/2022/01122022/'
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

numshots=10#94
direction_list=['r','t','z']
probelist=['probe5']#,'probe7','probe19','probe21','probe33','probe35']
directions = len(direction_list)
numprobes = 1#len(probelist)

delay_array = np.arange(1,10)
num_delays = len(delay_array)
embeddelay = 5
nfac = factorial(embeddelay)

PEs = np.zeros([num_delays])
SCs = np.zeros([num_delays])


#Need to convert hdf5 file into an array in order to pickle it for multiprocessing
data=np.array(data)

def pescy(data,loop_delay,embeddelay,numshots,probelist,direction_list):
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
                arr,nperms,embed = constructPatternCount(data1,embeddelay,delay=loop_delay)
                permstore_counter = permstore_counter+arr
                tot_perms = tot_perms+nperms
    PE_tot,PE_tot_Se = calcS_fromPatternCount(permstore_counter,tot_perms,embed)
    C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*np.log2(nfac))
             /((1 + 1./nfac)*np.log2(nfac+1) - 2*np.log2(2*nfac) 
            + np.log2(nfac))*(PE_tot/np.log2(nfac)))
    PE=PE_tot/np.log2(nfac)
    SC=C
    return PE,SC


from joblib import Parallel, delayed
import multiprocessing
#what are your inputs, and what operation do you want to perform
#on each input. For example
inputs = range(10)
def processInput(i):
    return i*i
num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=1)(delayed(pescy)(data,delay_array[loop_delay],embeddelay,numshots,probelist,direction_list) for loop_delay in np.arange(len(delay_array)))

p=np.array(results)
PEs2=p[:,0]
SCs2=p[:,1]

for loop_delay in np.arange(len(delay_array)):
    if (loop_delay%100)==0: print( 'On Delay ',delay_array[loop_delay])
    PEs[loop_delay],SCs[loop_delay]=pescy(data,delay_array[loop_delay],embeddelay,numshots,probelist,direction_list)    
"""
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
                arr,nperms,embed = constructPatternCount(data1,embeddelay,delay=delay_array[loop_delay])
                permstore_counter = permstore_counter+arr
                tot_perms = tot_perms+nperms
    PE_tot,PE_tot_Se = calcS_fromPatternCount(permstore_counter,tot_perms,embed)
    C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*np.log2(nfac))
             /((1 + 1./nfac)*np.log2(nfac+1) - 2*np.log2(2*nfac) 
            + np.log2(nfac))*(PE_tot/np.log2(nfac)))
    PEs[loop_delay]=PE_tot/np.log2(nfac)
    SCs[loop_delay]=C
"""

end=time.time()
totaltime=end-start
print('Time to Run \n'+str(totaltime))