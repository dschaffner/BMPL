
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
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\03262024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset_03262024.h5'
data=load_hdf5(directory+datafilename,verbose=True)

#magnetic probe time
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

#wire current time
time_s_wc=data['time_wc']['time_s'][:]

#The first probe has two loops in the same glass. I've listed this as probe 1a and probe1b
#probe 1a, probe_number=0 -> port 13 
#probe 1b, probe_number=1 -> port 13
#probe 2,  probe_number=2 -> port 15
#probe 3,  probe_number=3 -> port 19
#probe 4,  probe_number=4 -> port 21 (but only r and t directions recorded)


#magnetic data is organized in arrays of rad_pos, probe_number, shot_number, samples
#Example below extracts magnetic field data for r-direction, first (and only) radial position, probe 13, shot 70

directions = ['r','t','z']
arr=data['mag_probe'][directions[0]]['b'][0,0,70,:]
plt.figure(1)
plt.plot(timeB_s,arr)



#loading in the wire current data
plt.figure(2)
wire=data['wirecurrent']['wirecurrent'][0,40,:]
#plt.plot(time_s_wc,wire)
b=np.arange(2500)
Bfield_wire = (4*np.pi*1e-7*b/(2*np.pi*0.0254))*10000
Bfield_wire2 = (4*np.pi*1e-7*b/(2*np.pi*0.007))*10000
plt.plot(b,Bfield_wire)
plt.plot(b,Bfield_wire2)