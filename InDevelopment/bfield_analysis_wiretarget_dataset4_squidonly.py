

import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff
import numpy as np

# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\Dataset 4 Summer 2024 Wire Target\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_0kAwire_40shots.h5'
data1=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_2kAwire_40shots.h5'#without disk
data2=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_4kAwire_40shots.h5'#without disk
data3=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_6kAwire_40shots.h5'#without disk
data4=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_0kAwire_40shots.h5'#without disk
data5=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_8kAwire_40shots.h5'#without disk
data6=load_hdf5(directory+datafilename,verbose=True)

#magnetic probe time
time_s = data1['time']['time_s']
timeB_s = time_s[1:]
time_us = data1['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m


probes = 2
loops = 6
numshots=40

aveb_squid_1kV_0kA = np.zeros([probes,loops,25003])
aveb_squid_1kV_2kA = np.zeros([probes,loops,25003])
aveb_squid_1kV_4kA = np.zeros([probes,loops,25003])
aveb_squid_1kV_6kA = np.zeros([probes,loops,25003])
#spec_squid_frombdot_wodisk = np.zeros([probes,loops,numshots,fsize])

for probe in np.arange(probes):
    for loop in np.arange(loops):
        for shot in np.arange(numshots):
            datasquid1=data1['mag_probe']['squid']['b'][probe,loop,shot,:]
            aveb_squid_1kV_0kA[probe,loop,:]=aveb_squid_1kV_0kA[probe,loop,:]+np.sqrt(datasquid1**2)
            
            datasquid2=data2['mag_probe']['squid']['b'][probe,loop,shot,:]
            aveb_squid_1kV_2kA[probe,loop,:]=aveb_squid_1kV_2kA[probe,loop,:]+np.sqrt(datasquid2**2)
            
            datasquid3=data3['mag_probe']['squid']['b'][probe,loop,shot,:]
            aveb_squid_1kV_4kA[probe,loop,:]=aveb_squid_1kV_4kA[probe,loop,:]+np.sqrt(datasquid3**2)
            
            datasquid4=data4['mag_probe']['squid']['b'][probe,loop,shot,:]
            aveb_squid_1kV_6kA[probe,loop,:]=aveb_squid_1kV_6kA[probe,loop,:]+np.sqrt(datasquid4**2)
            
aveb_squid_1kV_0kA=aveb_squid_1kV_0kA/40.0
aveb_squid_1kV_2kA=aveb_squid_1kV_2kA/40.0
aveb_squid_1kV_4kA=aveb_squid_1kV_4kA/40.0
aveb_squid_1kV_6kA=aveb_squid_1kV_6kA/40.0

pickloop=0
pickprobe=0
plt.plot(timeB_us,aveb_squid_1kV_0kA[pickprobe,pickloop], label='0kA')
plt.plot(timeB_us,aveb_squid_1kV_2kA[pickprobe,pickloop], label='2kA')
plt.plot(timeB_us,aveb_squid_1kV_4kA[pickprobe,pickloop], label='4kA')
plt.plot(timeB_us,aveb_squid_1kV_6kA[pickprobe,pickloop], label='6kA')
plt.legend()
