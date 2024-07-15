
import scipy.io as spio
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff



# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06192024\\'
# ***** CHANGE Pathway!! *****************************************************

#datafilename= 'Dataset_06192024.h5'
datafilename= 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'
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
#probe 0 - port 1
#probe 1 - port 9 (theta direction bad)
#probe 2 - port 11
#probe 3 - port 13
#probe 4 - port 15
#probe 5 - port 19


#efield probe located at port 13

#There were 10 shots per current


#magnetic data is organized in arrays of rad_pos, current, probe_number, shot_number, samples
#Example below extracts magnetic field data for r-direction, first (and only) radial position, current 8, probe 15 and 17, shot 5 (of current 8)

#directions = ['r','t','z']
arr1=data['mag_probe']['r']['b'][0,5,:]
arr2=data['mag_probe']['t']['b'][0,5,:]
plt.figure(1)
plt.plot(timeB_us,arr1)
plt.plot(timeB_us,arr2)
