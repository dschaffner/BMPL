
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
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06262024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_nocurrent_20shots.h5'
data1=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_4kA_20shots.h5'#without disk
data2=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_1kV_1stuff_wire_nocurrent_20shots.h5'#without disk
data3=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_1kV_1stuff_wire_4kA_20shots.h5'#without disk
data4=load_hdf5(directory+datafilename,verbose=True)

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
numshots=20
directions=3
probes = 3

wavetest=data1['mag_probe']['t']['bdot'][0,0,:]

import pywt
wavelet = "cmor1.5-1.0"
# logarithmic scale for scales, as suggested by Torrence & Compo:
widths = np.geomspace(1, 1024*1, num=100)
sampling_period = np.diff(time_s[:]).mean()
cwtmatr, freqs = pywt.cwt(wavetest, widths, wavelet, sampling_period=sampling_period)
pwr_cwtmatr=np.abs(cwtmatr)**2
pwr_sum=np.sum(pwr_cwtmatr,axis=1)



f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(wavetest[:],time_s[:],window='hamming')

plt.loglog(f,1e27*pwr/f**2)

#plt.loglog(freqs,pwr_sum/freqs**2)

wv_n=100
Bpwr_60t160_para = np.zeros([wv_n])
Bpwr_60t160_perp = np.zeros([wv_n])
Bpwr_60t160_fft = np.zeros([12503])

for shot in np.arange(numshots):
    if shot==18: continue    
    print('Calculating Shot ',shot)
    datar=data3['mag_probe']['r']['bdot'][1,shot,:]
    Bwv_r, freqs = pywt.cwt(datar, widths, wavelet, sampling_period=sampling_period)
    Bpwr_r=np.abs(Bwv_r)**2
    for n in np.arange(len(time_s[:])):
        Bpwr_r[:,n] = Bpwr_r[:,n]/(freqs**2)
    datat=data3['mag_probe']['t']['bdot'][1,shot,:]    
    Bwv_t, freqs = pywt.cwt(datat, widths, wavelet, sampling_period=sampling_period)
    Bpwr_t=np.abs(Bwv_t)**2
    for n in np.arange(len(time_s[:])):
        Bpwr_t[:,n] = Bpwr_t[:,n]/(freqs**2)
    dataz=data3['mag_probe']['z']['bdot'][1,shot,:]    
    Bwv_z, freqs = pywt.cwt(dataz, widths, wavelet, sampling_period=sampling_period)
    Bpwr_z=np.abs(Bwv_z)**2
    for n in np.arange(len(time_s[:])):
        Bpwr_z[:,n] = Bpwr_z[:,n]/(freqs**2)
    
    Bmod=np.sqrt(datar**2+datat**2+dataz**2)
    
    #B vector percentages
    Br_para = (datar/Bmod)**2
    Bt_para = (datat/Bmod)**2
    Bz_para = (dataz/Bmod)**2
        
    Br_perp = 1.0-Br_para
    Bt_perp = 1.0-Bt_para
    Bz_perp = 1.0-Bz_para
        
    Bpara = np.zeros([wv_n,len(time_s[:])])
    Bperp = np.zeros([wv_n,len(time_s[:])])
    for w in np.arange(wv_n):
        Bpara[w,:] = (Bpwr_r[w,:]*Br_para)+(Bpwr_t[w,:]*Bt_para)+(Bpwr_z[w,:]*Bz_para)
        Bperp[w,:] = (Bpwr_r[w,:]*Br_perp)+(Bpwr_t[w,:]*Bt_perp)+(Bpwr_z[w,:]*Bz_perp)
        
    #para
    Bpwrslice = Bpara[:,start_time_index:end_time_index]
    Bpwrslice_sum = np.sum(Bpwrslice,axis=1)
    Bpwr_60t160_para = Bpwr_60t160_para+Bpwrslice_sum
    #perp
    Bpwrslice = Bperp[:,start_time_index:end_time_index]
    Bpwrslice_sum = np.sum(Bpwrslice,axis=1)
    Bpwr_60t160_perp = Bpwr_60t160_perp+Bpwrslice_sum
    
plt.loglog(freqs,Bpwr_60t160_para)
plt.loglog(freqs,Bpwr_60t160_perp)
#plt.contourf(time_s, np.log10(freqs), np.log10(pwr_cwtmatr)  )
"""
import compute_wavelet as cw

plt.figure(1)
plt.plot(wavetest)
Bpwr_r,Bscalespec,wvfreq,Bfft,fftfreq = cw.compute_wavelet(wavetest,time_us[:])
plt.figure(2)
plt.loglog(wvfreq,Bscalespec)
plt.loglog(fftfreq,Bfft)
"""
"""
#determine FFT size and generate an output array
fsize=int((data1['mag_probe']['r']['bdot'][0,0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot_2kV_0kA = np.zeros([probes,directions,fsize])
avebspec_frombdot_2kV_4kA = np.zeros([probes,directions,fsize])
avebspec_frombdot_1kV_0kA = np.zeros([probes,directions,fsize])
avebspec_frombdot_1kV_4kA = np.zeros([probes,directions,fsize])
#spec_frombdot_wdisk = np.zeros([probes,directions,numshots,fsize])
#spec_frombdot_wodisk = np.zeros([probes,directions,numshots,fsize])

direction_list = ['r','t','z']

#loop over shots to read in data and compute FFT

for probe in np.arange(probes):
    for direction_index, direction in enumerate(direction_list):
        for shot in np.arange(numshots):    
            data1ts=data1['mag_probe'][direction]['bdot'][probe,shot,:]
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

#fft of squid probes
probes = 2
loops = 6

avebspec_squid_frombdot_2kV_0kA = np.zeros([probes,loops,fsize])
avebspec_squid_frombdot_2kV_4kA = np.zeros([probes,loops,fsize])
avebspec_squid_frombdot_1kV_0kA = np.zeros([probes,loops,fsize])
avebspec_squid_frombdot_1kV_4kA = np.zeros([probes,loops,fsize])
#spec_squid_frombdot_wodisk = np.zeros([probes,loops,numshots,fsize])

for probe in np.arange(probes):
    for loop in np.arange(loops):
        for shot in np.arange(numshots):
            datasquid1=data1['mag_probe']['squid']['bdot'][probe,loop,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datasquid1[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_squid_frombdot_wodisk[probe,loop,shot,:]=pwr/(f*f)
            avebspec_squid_frombdot_2kV_0kA[probe,loop,:]=avebspec_squid_frombdot_2kV_0kA[probe,loop,:]+pwr/(f*f)

            datasquid2=data2['mag_probe']['squid']['bdot'][probe,loop,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datasquid2[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_squid_frombdot_wodisk[probe,loop,shot,:]=pwr/(f*f)
            avebspec_squid_frombdot_2kV_4kA[probe,loop,:]=avebspec_squid_frombdot_2kV_4kA[probe,loop,:]+pwr/(f*f)
            
            datasquid3=data3['mag_probe']['squid']['bdot'][probe,loop,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datasquid3[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_squid_frombdot_wodisk[probe,loop,shot,:]=pwr/(f*f)
            avebspec_squid_frombdot_1kV_0kA[probe,loop,:]=avebspec_squid_frombdot_1kV_0kA[probe,loop,:]+pwr/(f*f)

            datasquid4=data4['mag_probe']['squid']['bdot'][probe,loop,shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(datasquid4[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
            #spec_squid_frombdot_wodisk[probe,loop,shot,:]=pwr/(f*f)
            avebspec_squid_frombdot_1kV_4kA[probe,loop,:]=avebspec_squid_frombdot_1kV_4kA[probe,loop,:]+pwr/(f*f)
"""