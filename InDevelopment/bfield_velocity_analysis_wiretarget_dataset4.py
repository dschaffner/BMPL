
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

#Bdot Triplets
numshots=40
directions=3
probes = 3

#determine FFT size and generate an output array
fsize=int((data1['mag_probe']['r']['b'][0,0,:].shape[0]))

aveb_1kV_0kA = np.zeros([probes,directions,fsize])
aveb_1kV_2kA = np.zeros([probes,directions,fsize])
aveb_1kV_4kA = np.zeros([probes,directions,fsize])
aveb_1kV_8kA = np.zeros([probes,directions,fsize])
#spec_frombdot_wdisk = np.zeros([probes,directions,numshots,fsize])
#spec_frombdot_wodisk = np.zeros([probes,directions,numshots,fsize])

direction_list = ['r','t','z']
"""
#loop over shots to read in data and compute FFT
savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\Dataset 4 Summer 2024 Wire Target\\Plots\\Raw_Bfield\\1p5kV_8kA\\'
for probe in np.arange(probes):
    for direction_index, direction in enumerate(direction_list):
        for shot in np.arange(numshots):    
            #data1ts=data6['mag_probe'][direction]['b'][probe,shot,:]
            #plt.plot(timeB_us,data1ts)
            #plt.title('Probe '+str(probe+1)+' Direction '+direction+' Shot '+str(shot))
            #savefile = 'Probe'+str(probe+1)+'_Direction_'+direction+'_Shot'+str(shot)+'.png'
            #plt.savefig(savedirectory+savefile, dpi=300, facecolor='w', edgecolor='k')
            #plt.clf()
            

            aveb_2kV_0kA[probe,direction,:]=aveb_2kV_0kA[probe,direction,:]+
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_frombdot_2k_0kA[probe,direction_index,shot,:]=pwr/(f*f)
            avebspec_frombdot_2kV_0kA[probe,direction_index,:]=avebspec_frombdot_2kV_0kA[probe,direction_index,:]+pwr/(f*f)
            
            data2ts=data2['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data2ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_frombdot_wodisk[probe,direction_index,shot,:]=pwr/(f*f)       
            avebspec_frombdot_2kV_4kA[probe,direction_index,:]=avebspec_frombdot_2kV_4kA[probe,direction_index,:]+pwr/(f*f)
            
            data3ts=data3['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data3ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_frombdot_wodisk[probe,direction_index,shot,:]=pwr/(f*f)       
            avebspec_frombdot_1kV_0kA[probe,direction_index,:]=avebspec_frombdot_1kV_0kA[probe,direction_index,:]+pwr/(f*f)
            
            data4ts=data4['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data4ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_frombdot_wodisk[probe,direction_index,shot,:]=pwr/(f*f)       
            avebspec_frombdot_1kV_4kA[probe,direction_index,:]=avebspec_frombdot_1kV_4kA[probe,direction_index,:]+pwr/(f*f)
"""

#fft of squid probes
probes = 2
loops = 6

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