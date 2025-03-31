
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
analysis_start_time = 75
analysis_end_time = 175
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
fsize=int((data1['isat probe']['isat'][0,start_time_index:end_time_index].shape[0])/2)+1
ave_isat_1p5kV_0kA_spectra = np.zeros([fsize])
ave_isat_1p5kV_2kA_spectra = np.zeros([fsize])
ave_isat_1p5kV_4kA_spectra = np.zeros([fsize])
ave_isat_1p5kV_6kA_spectra = np.zeros([fsize])
ave_isat_1p5kV_8kA_spectra = np.zeros([fsize])


for shot in np.arange(numshots):
    dataisat1=data1['isat probe']['isat'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(dataisat1[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    ave_isat_1p5kV_0kA_spectra=ave_isat_1p5kV_0kA_spectra+pwr
    
    dataisat2=data2['isat probe']['isat'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(dataisat2[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    ave_isat_1p5kV_2kA_spectra=ave_isat_1p5kV_2kA_spectra+pwr
    
    dataisat3=data3['isat probe']['isat'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(dataisat3[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    ave_isat_1p5kV_4kA_spectra=ave_isat_1p5kV_4kA_spectra+pwr
    
    dataisat4=data4['isat probe']['isat'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(dataisat4[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    ave_isat_1p5kV_6kA_spectra=ave_isat_1p5kV_6kA_spectra+pwr
    
    dataisat5=data5['isat probe']['isat'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(dataisat5[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    ave_isat_1p5kV_8kA_spectra=ave_isat_1p5kV_8kA_spectra+pwr
    

ave_isat_1p5kV_0kA_spectra=ave_isat_1p5kV_0kA_spectra/10.0
ave_isat_1p5kV_2kA_spectra=ave_isat_1p5kV_2kA_spectra/10.0
ave_isat_1p5kV_4kA_spectra=ave_isat_1p5kV_4kA_spectra/10.0
ave_isat_1p5kV_6kA_spectra=ave_isat_1p5kV_6kA_spectra/10.0
ave_isat_1p5kV_8kA_spectra=ave_isat_1p5kV_8kA_spectra/10.0

plt.loglog(f,ave_isat_1p5kV_0kA_spectra,label='0kA')
plt.loglog(f,ave_isat_1p5kV_2kA_spectra,label='2kA')
plt.loglog(f,ave_isat_1p5kV_4kA_spectra,label='4kA')
plt.loglog(f,ave_isat_1p5kV_6kA_spectra,label='6kA')
plt.loglog(f,ave_isat_1p5kV_8kA_spectra,label='8kA')
plt.legend()


#Bfield spectra
directions=3
probes = 3

aveb_1p5kV_0kA_spectra = np.zeros([probes,directions,fsize])
aveb_1p5kV_2kA_spectra = np.zeros([probes,directions,fsize])
aveb_1p5kV_4kA_spectra = np.zeros([probes,directions,fsize])
aveb_1p5kV_6kA_spectra = np.zeros([probes,directions,fsize])
aveb_1p5kV_8kA_spectra = np.zeros([probes,directions,fsize])

direction_list = ['r','t','z']

#loop over shots to read in data and compute FFT

for probe in np.arange(probes):
    for direction_index, direction in enumerate(direction_list):
        for shot in np.arange(numshots):    
            datab1=data1['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datab1[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            aveb_1p5kV_0kA_spectra[probe,direction_index,:]=aveb_1p5kV_0kA_spectra[probe,direction_index,:]+pwr/(f*f)
            
            datab2=data2['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datab2[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            aveb_1p5kV_2kA_spectra[probe,direction_index,:]=aveb_1p5kV_2kA_spectra[probe,direction_index,:]+pwr/(f*f)

            datab3=data3['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datab3[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            aveb_1p5kV_4kA_spectra[probe,direction_index,:]=aveb_1p5kV_4kA_spectra[probe,direction_index,:]+pwr/(f*f)
            
            datab4=data4['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datab4[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            aveb_1p5kV_6kA_spectra[probe,direction_index,:]=aveb_1p5kV_6kA_spectra[probe,direction_index,:]+pwr/(f*f)
            
            datab5=data5['mag_probe'][direction]['bdot'][probe,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datab5[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            aveb_1p5kV_8kA_spectra[probe,direction_index,:]=aveb_1p5kV_8kA_spectra[probe,direction_index,:]+pwr/(f*f)

aveb_1p5kV_0kA_spectra=aveb_1p5kV_0kA_spectra/10.0
aveb_1p5kV_2kA_spectra=aveb_1p5kV_2kA_spectra/10.0
aveb_1p5kV_4kA_spectra=aveb_1p5kV_4kA_spectra/10.0
aveb_1p5kV_6kA_spectra=aveb_1p5kV_6kA_spectra/10.0
aveb_1p5kV_8kA_spectra=aveb_1p5kV_8kA_spectra/10.0

plt.figure(2)
plt.loglog(f,aveb_1p5kV_0kA_spectra[2,0,:]+aveb_1p5kV_0kA_spectra[2,1,:]+aveb_1p5kV_0kA_spectra[2,2,:],label='0kA')
plt.loglog(f,aveb_1p5kV_2kA_spectra[2,0,:]+aveb_1p5kV_2kA_spectra[2,1,:]+aveb_1p5kV_2kA_spectra[2,2,:],label='2kA')
plt.loglog(f,aveb_1p5kV_4kA_spectra[2,0,:]+aveb_1p5kV_4kA_spectra[2,1,:]+aveb_1p5kV_4kA_spectra[2,2,:],label='4kA')
plt.loglog(f,aveb_1p5kV_6kA_spectra[2,0,:]+aveb_1p5kV_6kA_spectra[2,1,:]+aveb_1p5kV_6kA_spectra[2,2,:],label='6kA')
plt.loglog(f,aveb_1p5kV_8kA_spectra[2,0,:]+aveb_1p5kV_8kA_spectra[2,1,:]+aveb_1p5kV_8kA_spectra[2,2,:],label='8kA')
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