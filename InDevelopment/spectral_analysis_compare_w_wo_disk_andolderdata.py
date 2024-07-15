
import scipy.io as spio
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff
import numpy as np
import spectrum_wwind as spec

#load dataset 1
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06192024\\'
datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
data1=load_hdf5(directory+datafilename,verbose=True)

#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06202024\\'
#datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
#data2=load_hdf5(directory+datafilename,verbose=True)


# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06262024\\'
# ***** CHANGE Pathway!! *****************************************************
#pick a file from four possible (by commenting out the other files)
datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_nocurrent_30shots.h5'#without disk
#datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_4kA_20shots.h5'#without disk
#datafilename= 'Dataset_06262024_1kV_1stuff_wire_nocurrent_20shots.h5'#without disk
#datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_4kA_20shots.h5'
#datafilename= 'Dataset_06262024_1kV_1stuff_wire_nocurrent_20shots.h5'
#datafilename= 'Dataset_06262024_1kV_1stuff_wire_4kA_20shots.h5'
data2=load_hdf5(directory+datafilename,verbose=True)

#magnetic probe time
time_s = data2['time']['time_s']
timeB_s = time_s[1:]
time_us = data2['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m

#Bdot Triplets
numshots=20
directions=3
probes = 3

#determine FFT size and generate an output array
fsize=int((data1['mag_probe']['r']['bdot'][0,0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot_wdisk = np.zeros([probes,directions,fsize])
avebspec_frombdot_wodisk = np.zeros([probes,directions,fsize])
spec_frombdot_wdisk = np.zeros([probes,directions,numshots,fsize])
spec_frombdot_wodisk = np.zeros([probes,directions,numshots,fsize])

direction_list = ['r','t','z']

#loop over shots to read in data and compute FFT

for probe in np.arange(probes):
    for direction_index, direction in enumerate(direction_list):
        for shot in np.arange(numshots):    
            data1ts=data1['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            spec_frombdot_wdisk[probe,direction_index,shot,:]=pwr/(f*f)
            avebspec_frombdot_wdisk[probe,direction_index,:]=avebspec_frombdot_wdisk[probe,direction_index,:]+pwr/(f*f)
            
            data2ts=data2['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data2ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            spec_frombdot_wodisk[probe,direction_index,shot,:]=pwr/(f*f)       
            avebspec_frombdot_wodisk[probe,direction_index,:]=avebspec_frombdot_wodisk[probe,direction_index,:]+pwr/(f*f)

#fft of squid probes
probes = 2
loops = 6

avebspec_squid_frombdot_wodisk = np.zeros([probes,loops,fsize])
spec_squid_frombdot_wodisk = np.zeros([probes,loops,numshots,fsize])

for probe in np.arange(probes):
    for loop in np.arange(loops):
        for shot in np.arange(numshots):
            datasquid=data2['mag_probe']['squid']['bdot'][probe,loop,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datasquid[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            spec_squid_frombdot_wodisk[probe,loop,shot,:]=pwr/(f*f)
            avebspec_squid_frombdot_wodisk[probe,loop,:]=avebspec_squid_frombdot_wodisk[probe,loop,:]+pwr/(f*f)

"""
#Data from 2022

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
datafilename='Dataset_01122022.h5'
data3=load_hdf5(directory+datafilename,verbose=True)

time_s = data3['time']['time_s']
timeB_s = time_s[1:]
time_us = data3['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 40
analysis_end_time = 60
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m

numshots=94
direction_list=['r','t','z']
probelist=['probe5','probe7','probe19','probe21','probe33','probe35']
directions = len(direction_list)
numprobes = len(probelist)

#determine FFT size and generate an output array
fsize=int((data3['mag_probe']['positions']['probe5']['r']['bdot'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot_old = np.zeros([numprobes,directions,fsize])
#avebspec_direct_old = np.zeros([numprobes,directions,fsize])
#avebmagspec = np.zeros([numprobes,fsize])
#spec_frombdot_old = np.zeros([numshots,numprobes,directions,fsize])
#spec_frombdot_sm = np.zeros([numshots,numprobes,directions,fsize])

#loop over shots to read in data and compute FFT
for shot in np.arange(numshots):
    for probe_index, probe in enumerate(probelist):
        for direction_index, direction in enumerate(direction_list):
            data3ts=data3['mag_probe']['positions'][probe][direction]['bdot'][shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data3ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hanning')
            #spec_frombdot[shot,probe_index,direction_index,:]=pwr/(f*f)
            avebspec_frombdot_old[probe_index,direction_index,:]=avebspec_frombdot_old[probe_index,direction_index]+(pwr/(f*f))

spec1=avebspec_frombdot_wodisk[0,1,1:]+avebspec_frombdot_wodisk[0,1,1:]#port1
spec2=avebspec_frombdot_wdisk[0,1,1:]+avebspec_frombdot_wdisk[0,1,1:]
spec3=avebspec_frombdot_old[0,1,1:]+avebspec_frombdot_old[0,1,1:]
spec4=avebspec_frombdot_wodisk[2,1,1:]+avebspec_frombdot_wodisk[2,1,1:]#port5

max1=np.max(avebspec_frombdot_wodisk[0,1,1:]+avebspec_frombdot_wodisk[0,1,1:])
max2=np.max(avebspec_frombdot_wdisk[0,1,1:]+avebspec_frombdot_wdisk[0,1,1:])
max3=np.max(avebspec_frombdot_old[0,1,1:]+avebspec_frombdot_old[0,1,1:])
max4=np.max(avebspec_frombdot_wodisk[2,1,1:]+avebspec_frombdot_wodisk[2,1,1:])#port5

plt.figure(1)
plt.loglog(f[1:],spec1/max1,label='wo-disk port1')
plt.loglog(f[1:],spec2/max2,label='w-disk port1')
plt.loglog(f[1:],spec3/max3,label='2022 data port5' )
plt.loglog(f[1:],spec4/max4,label='wo-disk port5')

plt.xlim(3e4,5e7)
plt.ylim(1e-10,2e0)
plt.legend()
"""