
import scipy.io as spio
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff
import numpy as np
import spectrum_wwind as spec

#load dataset 1
#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06192024\\'
#datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
#data1=load_hdf5(directory+datafilename,verbose=True)

#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06202024\\'
#datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
#data2=load_hdf5(directory+datafilename,verbose=True)


# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\07192024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_0kAwire_density_10shots.h5'
data1=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_2kAwire_density_10shots.h5'#without disk
data2=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_4kAwire_density_10shots.h5'#without disk
data3=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_6kAwire_density_10shots.h5'#without disk
data4=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_8kAwire_density_10shots.h5'#without disk
data5=load_hdf5(directory+datafilename,verbose=True)


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

#Bdot Triplets
numshots=10
directions=3
probes = 3


#fft of squid probes
probes = 1
loops = 6

ave_isat_1p5kV_0kA = np.zeros([25004])
ave_isat_1p5kV_2kA = np.zeros([25004])
ave_isat_1p5kV_4kA = np.zeros([25004])
ave_isat_1p5kV_6kA = np.zeros([25004])
ave_isat_1p5kV_8kA = np.zeros([25004])
#spec_squid_frombdot_wodisk = np.zeros([probes,loops,numshots,fsize])


for shot in np.arange(numshots):
    dataisat1=data1['isat probe']['isat'][shot,:]
    ave_isat_1p5kV_0kA=ave_isat_1p5kV_0kA+dataisat1
    
    dataisat2=data2['isat probe']['isat'][shot,:]
    ave_isat_1p5kV_2kA=ave_isat_1p5kV_2kA+dataisat2
    
    dataisat3=data3['isat probe']['isat'][shot,:]
    ave_isat_1p5kV_4kA=ave_isat_1p5kV_4kA+dataisat3
    
    dataisat4=data4['isat probe']['isat'][shot,:]
    ave_isat_1p5kV_6kA=ave_isat_1p5kV_6kA+dataisat4
    
    dataisat5=data5['isat probe']['isat'][shot,:]
    ave_isat_1p5kV_8kA=ave_isat_1p5kV_8kA+dataisat5

ave_isat_1p5kV_0kA=ave_isat_1p5kV_0kA/10.0
ave_isat_1p5kV_2kA=ave_isat_1p5kV_2kA/10.0
ave_isat_1p5kV_4kA=ave_isat_1p5kV_4kA/10.0
ave_isat_1p5kV_6kA=ave_isat_1p5kV_6kA/10.0
ave_isat_1p5kV_8kA=ave_isat_1p5kV_8kA/10.0

plt.figure(1)
plt.plot(time_us[:],ave_isat_1p5kV_0kA,label='0kA')
plt.plot(time_us[:],ave_isat_1p5kV_2kA,label='2kA')
plt.plot(time_us[:],ave_isat_1p5kV_4kA,label='4kA')
plt.plot(time_us[:],ave_isat_1p5kV_6kA,label='6kA')
plt.plot(time_us[:],ave_isat_1p5kV_8kA,label='8kA')
plt.legend()




#### Bfield Analysis

ave_b25_1p5kV_0kA = np.zeros([25003])
ave_b25_1p5kV_2kA = np.zeros([25003])
ave_b25_1p5kV_4kA = np.zeros([25003])
ave_b25_1p5kV_6kA = np.zeros([25003])
ave_b25_1p5kV_8kA = np.zeros([25003])

for shot in np.arange(numshots):
    databr=data1['mag_probe']['r']['b'][2,shot,:]
    databt=data1['mag_probe']['t']['b'][2,shot,:]
    databz=data1['mag_probe']['z']['b'][2,shot,:]
    dataB1 = np.sqrt(databr**2+databt**2+databz**2)
    ave_b25_1p5kV_0kA=ave_b25_1p5kV_0kA+dataB1
    
    databr=data2['mag_probe']['r']['b'][2,shot,:]
    databt=data2['mag_probe']['t']['b'][2,shot,:]
    databz=data2['mag_probe']['z']['b'][2,shot,:]
    dataB2 = np.sqrt(databr**2+databt**2+databz**2)
    ave_b25_1p5kV_2kA=ave_b25_1p5kV_2kA+dataB2
    
    databr=data3['mag_probe']['r']['b'][2,shot,:]
    databt=data3['mag_probe']['t']['b'][2,shot,:]
    databz=data3['mag_probe']['z']['b'][2,shot,:]
    dataB3 = np.sqrt(databr**2+databt**2+databz**2)
    ave_b25_1p5kV_4kA=ave_b25_1p5kV_4kA+dataB3
    
    databr=data4['mag_probe']['r']['b'][2,shot,:]
    databt=data4['mag_probe']['t']['b'][2,shot,:]
    databz=data4['mag_probe']['z']['b'][2,shot,:]
    dataB4 = np.sqrt(databr**2+databt**2+databz**2)
    ave_b25_1p5kV_6kA=ave_b25_1p5kV_6kA+dataB4
    
    databr=data5['mag_probe']['r']['b'][2,shot,:]
    databt=data5['mag_probe']['t']['b'][2,shot,:]
    databz=data5['mag_probe']['z']['b'][2,shot,:]
    dataB5 = np.sqrt(databr**2+databt**2+databz**2)
    ave_b25_1p5kV_8kA=ave_b25_1p5kV_8kA+dataB5

ave_b25_1p5kV_0kA=ave_b25_1p5kV_0kA/10.0
ave_b25_1p5kV_2kA=ave_b25_1p5kV_2kA/10.0
ave_b25_1p5kV_4kA=ave_b25_1p5kV_4kA/10.0
ave_b25_1p5kV_6kA=ave_b25_1p5kV_6kA/10.0
ave_b25_1p5kV_8kA=ave_b25_1p5kV_8kA/10.0

plt.figure(2)
plt.plot(timeB_us[:],ave_b25_1p5kV_0kA,label='0kA')
plt.plot(timeB_us[:],ave_b25_1p5kV_2kA,label='2kA')
plt.plot(timeB_us[:],ave_b25_1p5kV_4kA,label='4kA')
plt.plot(timeB_us[:],ave_b25_1p5kV_6kA,label='6kA')
plt.plot(timeB_us[:],ave_b25_1p5kV_8kA,label='8kA')
plt.legend()
"""
fig = plt.figure(num=2, figsize=(5, 3), dpi=300,
                 facecolor='w', edgecolor='k')

# call first subplot object
ax = plt.subplot(2, 1, 1)  # (num rows, num columns, subplot position)
plt.plot(timeB_us[:],aveb_squid_1kV_0kA[1,2,:],label='0kA')
plt.plot(timeB_us[:],aveb_squid_1kV_4kA[1,2,:],label='4kA')
plt.xlabel('Time[us]', fontsize=7)
plt.xticks(fontsize=7)
plt.ylabel('Magnetic Field [G]', fontsize=7)
plt.yticks(fontsize=7)
plt.legend(fontsize=4)

ax2 = plt.subplot(2,1,2)
plt.plot(timeB_us[:],aveb_squid_2kV_0kA[1,2,:],label='0kA')
plt.plot(timeB_us[:],aveb_squid_2kV_4kA[1,2,:],label='4kA')
plt.xlabel('Time[us]', fontsize=7)
plt.xticks(fontsize=7)
plt.ylabel('Magnetic Field [G]', fontsize=7)
plt.yticks(fontsize=7)
plt.legend(fontsize=4)
"""