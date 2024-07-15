
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
from scipy.signal import butter, sosfiltfilt, sosfreqz


# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06052024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset_06052024.h5'
data=load_hdf5(directory+datafilename,verbose=True)

#magnetic probe time
time_s = data['time']['time_s']
timeB_s = time_s[1:]
time_us = data['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 100
analysis_end_time = 140
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m

#wire current time
#time_s_wc=data['time_wc']['time_s'][:]

#probe layout
#probe 0 - port 15 in front of target wire
#probe 1 - port 17 right against target wire
#probe 2 - probe 19 past target wire
#probe 3 - probe 21 past target wire

#current scan
#current 0 = 0A
#current 1 = 500A
#current 2 = 750A
#current 3 = 1000A
#current 4 = 1250A
#current 5 = 1500A
#current 6 = 1750A
#current 7 = 2000A
#current 8 = 2250A
#current 9 = 2500A

#vf probe located at port 17

#There were 10 shots per current


#magnetic data is organized in arrays of rad_pos, current, probe_number, shot_number, samples
#Example below extracts magnetic field data for r-direction, first (and only) radial position, current 8, probe 15 and 17, shot 5 (of current 8)

directions = ['r','t','z']
arr1=data['mag_probe'][directions[0]]['b'][0,8,0,5,:]
arr2=data['mag_probe'][directions[0]]['b'][0,8,1,5,:]
plt.figure(1)
plt.plot(timeB_s,arr1)
plt.plot(timeB_s,arr2)
