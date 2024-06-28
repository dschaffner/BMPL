
import scipy.io as spio
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff



# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06262024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_nocurrent_20shots.h5'
datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_4kA_20shots.h5'
datafilename= 'Dataset_06262024_1kV_1stuff_wire_nocurrent_20shots.h5'
datafilename= 'Dataset_06262024_1kV_1stuff_wire_4kA_20shots.h5'
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

"""
5 probes:
    Triplet Bdot at Port 1, 3, 5
    Squid Probe at Ports 27 and 29 (only one loop in z direction)
    Electric field probe at opposite side port 2 (postive efield points away from gun)
    
    Bdots at 2V/div on scope
"""



#bdot triplet probes
#organized by direction, b or bdot, then port position, shot number
arr1=data['mag_probe']['t']['b'][0,5,:]#probe in port 1 t-direction, shot 5
arr2=data['mag_probe']['r']['b'][1,12,:]#probe in port 3 r-direction shot 12
plt.figure(1)
plt.plot(timeB_us[:],arr1)
plt.plot(timeB_us[:],arr2)

#bdot squid probes
#organized by b or bdot, then probe position (port 27 or 29), then loop number (1-6), then shot number
arr3=data['mag_probe']['squid']['b'][0,4,10,:]#squid probe in port 27,loop 4, shot 10
arr4=data['mag_probe']['squid']['b'][1,3,18,:]#squid probe in port 29, loop 3, shot 18
plt.figure(2)
plt.plot(timeB_us,arr3)
plt.plot(timeB_us,arr4)
