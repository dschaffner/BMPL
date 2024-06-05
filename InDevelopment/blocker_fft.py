
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

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y

import smooth as sm

def iter_smooth(array,loops=6,window_len=3):
    for l in np.arange(loops):
        array = sm.smooth(array,window_len=window_len)
    return array


# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2023\\07182023\\'
#datafilename='Dataset_06282023_course_radscan.h5'
#datafilename='Dataset_07132023_blockertest.h5'
datafilename= 'Dataset_07182023_full_radscan.h5'
data=load_hdf5(directory+datafilename,verbose=True)


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

numshots=10




directions = ['r','t','z']
spectra=np.zeros([7,5,20,3,6251])
for r in np.arange(7):
    for probe in np.arange(5):
        for n in np.arange(20):
            for d in np.arange(len(directions)): 
                arr=data['mag_probe'][directions[d]]['b'][r,probe,n,:]
                f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(arr[start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
                spectra[r,probe,n,d,:]=pwr
                
spectra_dirsum = np.sum(spectra,axis=3)
spectra_shotsum_dirsum = np.sum(spectra_dirsum,axis=2)
spectra_shotsum = np.sum(spectra,axis=2)         
                

        
     






"""
cs=plt.contourf(x,y,z,levels=100,cmap='plasma')
cbar=plt.colorbar()
cbar.set_label(r'$\Delta B / B$', rotation=270, labelpad=15)

plt.title(r"Average $\Delta B/B$ over 20 shots, 60-160us",fontsize=9)
plt.ylabel("Radial Position from Central Axis [in]",fontsize=8)
plt.xlabel("Axial Distance from Block [in]", fontsize=8)
plt.xlim(-1.0,3)
plt.xticks(np.array([-1,0,1,2,3]), [
           0,1,2,3,4], fontsize=9)

plt.vlines(-1, 1.5, 3, color='grey', linestyle='solid', linewidth=50.0)

savefile = 'blocker_scan_contour_deltaBoverB.png'
plt.savefig(savefile, dpi=300, facecolor='w', edgecolor='k')


#plt.yticks(np.array([0, 0.5, 1, 1.5,2, 2.5, 3]), fontsize=9)
"""