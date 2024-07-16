# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:49:14 2022

@author: Josh0
"""

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
from scipy.signal import savgol_filter as sv
from scipy.stats import norm

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

fs = 125e6
lowcut = 50e3
highcut = 2e6

#plt.style.use('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
plt.style.use('C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots_new.py')

directory='C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
datafilename1='1122022_111shots_somebad_011822.h5'
datafilename2='7152022_3p5kV_51shots.h5'
datafilename='7202022_3p5kV_203shots_Mixof225V_250V_or5ms_Nozzle.h5'

datafilename923='9232022_3p5kV_NozzleSweep_144shots.h5'

data=load_hdf5(directory+datafilename,verbose=True)
data1=load_hdf5(directory+datafilename1,verbose=True)
data923=load_hdf5(directory+datafilename923,verbose=True)

time_s = data['time']['time_s'][:]
timeB_s = data['time']['timeB_s'][:]
time_us = data['time']['time_us'][:]
timeB_us = data['time']['timeB_us'][:]
analysis_start_time = 50
analysis_end_time = 150
start_time_index = iff.tindex_min(analysis_start_time,timeB_us[0,:])
end_time_index = iff.tindex_min(analysis_end_time,timeB_us[0,:])
timerange_limit = 3e-6#s
port_sep = 0.0254#m

numshots=30
direction_list=['r','theta','z']
num_pairs = 2
directions = 3
tde_pairs = [['pos19','pos21']]

Bmag_means5 = np.zeros(96)
Bmag_means7 = np.zeros(96)
Bmag_means19 = np.zeros(96)
Bmag_means21 = np.zeros(96)
Bmag_means33 = np.zeros(96)
Bmag_means35 = np.zeros(96)

v_alfven5 = np.zeros(96)
v_alfven7 = np.zeros(96)
v_alfven19 = np.zeros(96)
v_alfven21 = np.zeros(96)
v_alfven33 = np.zeros(96)
v_alfven35 = np.zeros(96)

Bmag1_means5 = np.zeros(111)
Bmag1_means7 = np.zeros(111)
Bmag1_means19 = np.zeros(111)
Bmag1_means21 = np.zeros(111)
Bmag1_means33 = np.zeros(111)
Bmag1_means35 = np.zeros(111)

v1_alfven5 = np.zeros(111)
v1_alfven7 = np.zeros(111)
v1_alfven19 = np.zeros(111)
v1_alfven21 = np.zeros(111)
v1_alfven33 = np.zeros(111)
v1_alfven35 = np.zeros(111)

Bmag2_means5 = np.zeros(51)
Bmag2_means7 = np.zeros(51)
Bmag2_means19 = np.zeros(51)
Bmag2_means21 = np.zeros(51)
Bmag2_means33 = np.zeros(51)
Bmag2_means35 = np.zeros(51)

v2_alfven5 = np.zeros(51)
v2_alfven7 = np.zeros(51)
v2_alfven19 = np.zeros(51)
v2_alfven21 = np.zeros(51)
v2_alfven33 = np.zeros(51)
v2_alfven35 = np.zeros(51)

Bmag923_means5 = np.zeros(145)
Bmag923_means6 = np.zeros(145)
Bmag923_means7 = np.zeros(145)
Bmag923_means8 = np.zeros(145)

v923_alfven5 = np.zeros(145)
v923_alfven6 = np.zeros(145)
v923_alfven7 = np.zeros(145)
v923_alfven8 = np.zeros(145)

plotshots=1
"""
plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2.5,markeredgewidth=0.0)
"""
for shot in np.arange(1,145):
    
    datar5_9=data923['pos5']['B']['r'][shot,start_time_index:end_time_index]
    datat5_9=data923['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    dataz5_9=data923['pos5']['B']['z'][shot,start_time_index:end_time_index]
    datamod5_9=np.sqrt(datar5_9**2+datat5_9**2+dataz5_9**2)

    datar6_9=data923['pos5']['B']['r'][shot,start_time_index:end_time_index]
    datat6_9=data923['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    dataz6_9=data923['pos5']['B']['z'][shot,start_time_index:end_time_index]
    datamod6_9=np.sqrt(datar6_9**2+datat6_9**2+dataz6_9**2)

    datar7_9=data923['pos7']['B']['r'][shot,start_time_index:end_time_index]
    datat7_9=data923['pos7']['B']['theta'][shot,start_time_index:end_time_index]
    dataz7_9=data923['pos7']['B']['z'][shot,start_time_index:end_time_index]
    datamod7_9=np.sqrt(datar7_9**2+datat7_9**2+dataz7_9**2)

    datar8_9=data923['pos5']['B']['r'][shot,start_time_index:end_time_index]
    datat8_9=data923['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    dataz8_9=data923['pos5']['B']['z'][shot,start_time_index:end_time_index]
    datamod8_9=np.sqrt(datar8_9**2+datat8_9**2+dataz8_9**2)


    Bmag923_means5[shot] = np.mean(datamod5_9[shot])
    Bmag923_means6[shot] = np.mean(datamod6_9[shot])
    Bmag923_means7[shot] = np.mean(datamod7_9[shot])
    Bmag923_means8[shot] = np.mean(datamod8_9[shot])

    alfven_const = 0.004587
    v923_alfven5[shot] = Bmag923_means5[shot]/alfven_const*1e-3
    v923_alfven6[shot] = Bmag923_means6[shot]/alfven_const*1e-3
    v923_alfven7[shot] = Bmag923_means7[shot]/alfven_const*1e-3
    v923_alfven8[shot] = Bmag923_means8[shot]/alfven_const*1e-3


for shot in np.arange(1,96):

    datar5=data['pos5']['B']['r'][shot,start_time_index:end_time_index]
    
    #datar5=data['pos5']['B']['r'][40,1000:250000]

    
    datat5=data['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    dataz5=data['pos5']['B']['z'][shot,start_time_index:end_time_index]
    datamod5=np.sqrt(datar5**2+datat5**2+dataz5**2)
    
    datar7=data['pos7']['B']['r'][shot,start_time_index:end_time_index]
    datat7=data['pos7']['B']['theta'][shot,start_time_index:end_time_index]
    dataz7=data['pos7']['B']['z'][shot,start_time_index:end_time_index]
    datamod7=np.sqrt(datar7**2+datat7**2+dataz7**2)

    
    datar19=data['pos19']['B']['r'][shot,start_time_index:end_time_index]
    datat19=data['pos19']['B']['theta'][shot,start_time_index:end_time_index]
    dataz19=data['pos19']['B']['z'][shot,start_time_index:end_time_index]
    datamod19=np.sqrt(datar19**2+datat19**2+dataz19**2)
    
    datar21=data['pos21']['B']['r'][shot,start_time_index:end_time_index]
    datat21=data['pos21']['B']['theta'][shot,start_time_index:end_time_index]
    dataz21=data['pos21']['B']['z'][shot,start_time_index:end_time_index]
    datamod21=np.sqrt(datar21**2+datat21**2+dataz21**2)

    datar33=data['pos33']['B']['r'][shot,start_time_index:end_time_index]
    datat33=data['pos33']['B']['theta'][shot,start_time_index:end_time_index]
    dataz33=data['pos33']['B']['z'][shot,start_time_index:end_time_index]
    datamod33=np.sqrt(datar33**2+datat33**2+dataz33**2)
    
    datar35=data['pos35']['B']['r'][shot,start_time_index:end_time_index]
    datat35=data['pos35']['B']['theta'][shot,start_time_index:end_time_index]
    dataz35=data['pos35']['B']['z'][shot,start_time_index:end_time_index]
    datamod35=np.sqrt(datar35**2+datat35**2+dataz35**2)

    Bmag_means5[shot] = np.mean(datamod5[shot])
    Bmag_means7[shot] = np.mean(datamod7[shot])
    Bmag_means19[shot] = np.mean(datamod19[shot])
    Bmag_means21[shot] = np.mean(datamod21[shot])
    Bmag_means33[shot] = np.mean(datamod33[shot])
    Bmag_means35[shot] = np.mean(datamod35[shot])

    alfven_const = 0.004587
    v_alfven5[shot] = Bmag_means5[shot]/alfven_const*1e-3
    v_alfven7[shot] = Bmag_means7[shot]/alfven_const*1e-3
    v_alfven19[shot] = Bmag_means19[shot]/alfven_const*1e-3
    v_alfven21[shot] = Bmag_means21[shot]/alfven_const*1e-3
    v_alfven33[shot] = Bmag_means33[shot]/alfven_const*1e-3
    v_alfven35[shot] = Bmag_means35[shot]/alfven_const*1e-3

for shot in np.arange(1,111):

    data1r5=data1['pos5']['B']['r'][shot,start_time_index:end_time_index]
    
    #datar5=data['pos5']['B']['r'][40,1000:250000]

    
    data1t5=data1['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    data1z5=data1['pos5']['B']['z'][shot,start_time_index:end_time_index]
    data1mod5=np.sqrt(data1r5**2+data1t5**2+data1z5**2)
    
    data1r7=data1['pos7']['B']['r'][shot,start_time_index:end_time_index]
    data1t7=data1['pos7']['B']['theta'][shot,start_time_index:end_time_index]
    data1z7=data1['pos7']['B']['z'][shot,start_time_index:end_time_index]
    data1mod7=np.sqrt(data1r7**2+data1t7**2+data1z7**2)

    
    data1r19=data['pos19']['B']['r'][shot,start_time_index:end_time_index]
    data1t19=data['pos19']['B']['theta'][shot,start_time_index:end_time_index]
    data1z19=data['pos19']['B']['z'][shot,start_time_index:end_time_index]
    data1mod19=np.sqrt(data1r19**2+data1t19**2+data1z19**2)
    
    data1r21=data['pos21']['B']['r'][shot,start_time_index:end_time_index]
    data1t21=data['pos21']['B']['theta'][shot,start_time_index:end_time_index]
    data1z21=data['pos21']['B']['z'][shot,start_time_index:end_time_index]
    data1mod21=np.sqrt(data1r21**2+data1t21**2+data1z21**2)

    data1r33=data['pos33']['B']['r'][shot,start_time_index:end_time_index]
    data1t33=data['pos33']['B']['theta'][shot,start_time_index:end_time_index]
    data1z33=data['pos33']['B']['z'][shot,start_time_index:end_time_index]
    data1mod33=np.sqrt(data1r33**2+data1t33**2+data1z33**2)
    
    data1r35=data['pos35']['B']['r'][shot,start_time_index:end_time_index]
    data1t35=data['pos35']['B']['theta'][shot,start_time_index:end_time_index]
    data1z35=data['pos35']['B']['z'][shot,start_time_index:end_time_index]
    data1mod35=np.sqrt(data1r35**2+data1t35**2+data1z35**2)

    Bmag1_means5[shot] = np.mean(data1mod5[shot])
    Bmag1_means7[shot] = np.mean(data1mod7[shot])
    Bmag1_means19[shot] = np.mean(data1mod19[shot])
    Bmag1_means21[shot] = np.mean(data1mod21[shot])
    Bmag1_means33[shot] = np.mean(data1mod33[shot])
    Bmag1_means35[shot] = np.mean(data1mod35[shot])

    alfven_const = 0.004587
    v1_alfven5[shot] = Bmag1_means5[shot]/alfven_const*1e-3
    v1_alfven7[shot] = Bmag1_means7[shot]/alfven_const*1e-3
    v1_alfven19[shot] = Bmag1_means19[shot]/alfven_const*1e-3
    v1_alfven21[shot] = Bmag1_means21[shot]/alfven_const*1e-3
    v1_alfven33[shot] = Bmag1_means33[shot]/alfven_const*1e-3
    v1_alfven35[shot] = Bmag1_means35[shot]/alfven_const*1e-3
    
for shot in np.arange(1,51):

    data2r5=data1['pos5']['B']['r'][shot,start_time_index:end_time_index]
    
    #datar5=data['pos5']['B']['r'][40,1000:250000]

    
    data2t5=data1['pos5']['B']['theta'][shot,start_time_index:end_time_index]
    data2z5=data1['pos5']['B']['z'][shot,start_time_index:end_time_index]
    data2mod5=np.sqrt(data2r5**2+data2t5**2+data2z5**2)
    
    data2r7=data1['pos7']['B']['r'][shot,start_time_index:end_time_index]
    data2t7=data1['pos7']['B']['theta'][shot,start_time_index:end_time_index]
    data2z7=data1['pos7']['B']['z'][shot,start_time_index:end_time_index]
    data2mod7=np.sqrt(data2r7**2+data2t7**2+data2z7**2)

    
    data2r19=data['pos19']['B']['r'][shot,start_time_index:end_time_index]
    data2t19=data['pos19']['B']['theta'][shot,start_time_index:end_time_index]
    data2z19=data['pos19']['B']['z'][shot,start_time_index:end_time_index]
    data2mod19=np.sqrt(data2r19**2+data2t19**2+data2z19**2)
    
    data2r21=data['pos21']['B']['r'][shot,start_time_index:end_time_index]
    data2t21=data['pos21']['B']['theta'][shot,start_time_index:end_time_index]
    data2z21=data['pos21']['B']['z'][shot,start_time_index:end_time_index]
    data2mod21=np.sqrt(data2r21**2+data2t21**2+data2z21**2)

    data2r33=data['pos33']['B']['r'][shot,start_time_index:end_time_index]
    data2t33=data['pos33']['B']['theta'][shot,start_time_index:end_time_index]
    data2z33=data['pos33']['B']['z'][shot,start_time_index:end_time_index]
    data2mod33=np.sqrt(data2r33**2+data2t33**2+data2z33**2)
    
    data2r35=data['pos35']['B']['r'][shot,start_time_index:end_time_index]
    data2t35=data['pos35']['B']['theta'][shot,start_time_index:end_time_index]
    data2z35=data['pos35']['B']['z'][shot,start_time_index:end_time_index]
    data2mod35=np.sqrt(data2r35**2+data2t35**2+data2z35**2)

    Bmag2_means5[shot] = np.mean(data2mod5[shot])
    Bmag2_means7[shot] = np.mean(data2mod7[shot])
    Bmag2_means19[shot] = np.mean(data2mod19[shot])
    Bmag2_means21[shot] = np.mean(data2mod21[shot])
    Bmag2_means33[shot] = np.mean(data2mod33[shot])
    Bmag2_means35[shot] = np.mean(data2mod35[shot])

    alfven_const = 0.004587
    v2_alfven5[shot] = Bmag2_means5[shot]/alfven_const*1e-3
    v2_alfven7[shot] = Bmag2_means7[shot]/alfven_const*1e-3
    v2_alfven19[shot] = Bmag2_means19[shot]/alfven_const*1e-3
    v2_alfven21[shot] = Bmag2_means21[shot]/alfven_const*1e-3
    v2_alfven33[shot] = Bmag2_means33[shot]/alfven_const*1e-3
    v2_alfven35[shot] = Bmag2_means35[shot]/alfven_const*1e-3


#valf_3p5kV_120mT = np.array([v_alfven5[0:61],v_alfven7[0:61],v_alfven19[0:61],v_alfven21[0:61],v_alfven33[0:61],v_alfven35[0:61]])
valf_2kV_0mT = np.concatenate((v1_alfven5[0:111],v1_alfven7[0:111],v1_alfven19[0:111],v1_alfven21[0:111],v1_alfven33[0:111],v1_alfven35[0:111]))
valf_2kV_57 = np.concatenate((v1_alfven5[0:111], v1_alfven7[0:111]))
valf_2kV_3335 = np.concatenate((v1_alfven33[0:111], v1_alfven35[0:111]))

valf_3p5kV_110mT = np.concatenate((v2_alfven5[0:51],v2_alfven7[0:51],v2_alfven19[0:51],v2_alfven21[0:51],v2_alfven33[0:51],v2_alfven35[0:51]))
valf_3p5_110_57 = np.concatenate((v2_alfven5[0:51], v2_alfven7[0:51]))                                
valf_3p5_110_3335 = np.concatenate((v2_alfven33[0:51], v2_alfven35[0:51]))  

## ^ This keeps the dimensionality of each array, so its like 6 arrays each of 61 values.  Put them all into an array with concatenate below.
valf_3p5kV_120mT = np.concatenate((v_alfven5[0:61],v_alfven7[0:61],v_alfven19[0:61],v_alfven21[0:61],v_alfven33[0:61],v_alfven35[0:61]))
valf_3p5_120_57 = np.concatenate((v_alfven5[0:61], v_alfven7[0:61]))
valf_3p5_120_3335 = np.concatenate((v_alfven33[0:61], v_alfven35[0:61]))

valf_3p5kV_133mT = np.concatenate((v_alfven5[61:96],v_alfven7[61:96],v_alfven19[61:96],v_alfven21[61:96],v_alfven33[61:96],v_alfven35[61:96]))
valf_3p5_133_57 = np.concatenate((v_alfven5[61:96], v_alfven7[61:96]))
valf_3p5_133_3335 = np.concatenate((v_alfven33[61:96], v_alfven35[61:96]))

# Collection on 9/23 had positions 5,6,7,8 only, no 33/35
valf_3p5kV_0mT = np.concatenate((v923_alfven5[99:102],v923_alfven6[99:102],v923_alfven7[99:102],v923_alfven8[99:102],v923_alfven5[113:119],v923_alfven6[113:119],v923_alfven7[113:119],v923_alfven8[113:119],v923_alfven5[134:144],v923_alfven6[134:144],v923_alfven7[134:144],v923_alfven8[134:144]))
valf_3p5_0_57 = np.concatenate((v923_alfven5[99:102], v923_alfven7[99:102]))
#valf_3p5_0_3335 = np.concatenate((v923_alfven33[99:102], v923_alfven35[99:102]))

valf_0mix_57 = np.concatenate((valf_2kV_57, valf_3p5_0_57))
#valf_0mix_3335 = np.concatenate((valf_2kv_3335, valf_3p5_0_3335))

valf_globalmix_57 = np.concatenate((valf_0mix_57, valf_3p5_110_57, valf_3p5_120_57, valf_3p5_133_57))
valf_globalmix_3335 = np.concatenate((valf_2kV_3335, valf_3p5_110_3335, valf_3p5_120_3335, valf_3p5_133_3335))
#valf_3p5kV_0mT = np.concatenate((v923_alfven5[134:144],v923_alfven6[134:144],v923_alfven7[134:144],v923_alfven8[134:144]))
# ***** Save this in case the data gets wonky!! ^^^^^^ ********

(mu923, sigma923) = norm.fit(valf_3p5kV_0mT)
print(round(mu923), round(sigma923))
bins = np.array(range(100, 1100, 75))
bins2 = np.array(range(100, 1100, 75))
#cut_valfven5 = v_alfven5[(v_alfven5 >= 100) & (v_alfven5 <= 1000)]
#cut_valfven7 = v_alfven7[(v_alfven7 >= 100) & (v_alfven7 <= 1000)]
#cut_valfven19 = v_alfven5[(v_alfven19 >= 100) & (v_alfven19 <= 1000)]
#cut_valfven21 = v_alfven5[(v_alfven21 >= 100) & (v_alfven21 <= 1000)]
#cut_valfven33 = v_alfven5[(v_alfven33 >= 100) & (v_alfven33 <= 1000)]
#cut_valfven35 = v_alfven5[(v_alfven35 >= 100) & (v_alfven35 <= 1000)]

cut_valf_2kV_0mT = valf_2kV_0mT[(valf_2kV_0mT >= 100) & (valf_2kV_0mT <=1100)]
cut_valf_3p5kV_0mT = valf_3p5kV_0mT[(valf_3p5kV_0mT >= 100) & (valf_3p5kV_0mT <=1100)]
cut_valf_3p5kV_110mT = valf_3p5kV_110mT[(valf_3p5kV_110mT >= 100) & (valf_3p5kV_110mT <= 1100)]
cut_valf_3p5kV_120mT = valf_3p5kV_120mT[(valf_3p5kV_120mT >= 100) & (valf_3p5kV_120mT <= 1100)]
cut_valf_3p5kV_133mT = valf_3p5kV_133mT[(valf_3p5kV_133mT >= 100) & (valf_3p5kV_133mT <= 1100)]

cut_valf_23p5kV_0mT = np.concatenate((cut_valf_2kV_0mT, cut_valf_3p5kV_0mT))

cut_valf_global57 = valf_globalmix_57[(valf_globalmix_57 >= 100) & (valf_globalmix_57 <= 1100)]
cut_valf_global3335 = valf_globalmix_3335[(valf_globalmix_3335 >= 100) & (valf_globalmix_3335 <= 1100)]

"""
(mu1, sigma1) = norm.fit(valf_3p5kV_120mT)# !=0)
(mu2, sigma2) = norm.fit(valf_3p5kV_133mT)
(mu3, sigma3) = norm.fit(valf_2kV_0mT)
(mu4, sigma4) = norm.fit(valf_3p5kV_110mT)
"""
(mu1, sigma1) = norm.fit(cut_valf_3p5kV_120mT)# !=0)
(mu2, sigma2) = norm.fit(cut_valf_3p5kV_133mT)
(mu4, sigma4) = norm.fit(cut_valf_3p5kV_110mT)
(mu923, sigma923) = norm.fit(cut_valf_3p5kV_0mT)
(mu2kv, sigma2kv) = norm.fit(cut_valf_2kV_0mT)
(mu0mT, sigma0mT) = norm.fit(cut_valf_23p5kV_0mT)

(mu57, sigma57) = norm.fit(cut_valf_global57)
(mu3335, sigma3335) = norm.fit(cut_valf_global3335)

#print(mu3, mu4, mu1, mu2)
#print(sigma3, sigma4, sigma1, sigma2)

#(mu33, sigma33) = norm.fit(v_alfven33)
#(mu35, sigma35) = norm.fit(v_alfven35)

#mu5_round = round(mu5)
#mu7_round = round(mu7)
##mu19_round = round(mu19)
#mu21_round = round(mu21)
#mu33_round = round(mu33)
#mu35_round = round(mu35)

stdev_alf = np.std([v_alfven5,v_alfven7, v_alfven19, v_alfven21,v_alfven33, v_alfven35])
mean_alf=np.mean([v_alfven5,v_alfven7, v_alfven19, v_alfven21,v_alfven33, v_alfven35])


#fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6)

"""
fig, (ax1, ax2) = plt.subplots(2)

ax1.hist(cut_valf_global57, bins=bins, color='blue', edgecolor='black', label='Pos5/7, All Nozzles', align='right', alpha=0.7)
ax1.axvline(x=mu57, color='black', label='Mean =' + str(round(mu57)) + ' km s$^{-1}$', linestyle='--', markersize=13)
ax1.set_ylabel('Counts')
ax1.set_xlim(-100, 1400)
ax1.set_ylim(0, 60)
ax1.legend(loc='best', frameon=False)

ax2.hist(cut_valf_global3335, bins=bins, color='blue', edgecolor='black', label='Pos33/35, All Nozzles', align='right', alpha=0.7)
ax2.axvline(x=mu3335, color='black', label='Mean =' + str(round(mu3335)) + ' km s$^{-1}$', linestyle='--', markersize=13)
ax2.set_ylabel('Counts')
ax2.set_xlim(-100, 1400)
ax2.set_ylim(0, 70)
ax2.legend(loc='best', frameon=False)
ax2.set_xlabel(r'$Alfv\acute{e}n~Speed~V_a~(km~s^{-1})$')

#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\AlfvenSpeedStudies\\AlfvenSpeed_PositionComparison_57_3335.png', dpi=300)
"""
"""
ax0.hist(cut_valf_2kV_0mT, bins=bins, color='green', edgecolor='black', label='2kv, 0mT Nozz', align='right', alpha=0.7)
ax0.axvline(x=mu2kv, color='black', label='Mean =' + str(round(mu2kv)) + ' km s$^{-1}$', linestyle='--', markersize=13)
"""

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)


ax1.hist(cut_valf_23p5kV_0mT, bins=bins2, color='purple', edgecolor='black', label='0mT Nozzle', align='right')# alpha=0.7)
ax1.axvline(x=mu0mT, color='black', label= 'Mean =' + str(round(mu0mT)) + ' km s$^{-1}$', linestyle='--', markersize=13)

ax2.hist(cut_valf_3p5kV_110mT, bins=bins, color='purple', edgecolor='black', label=r'100mT Nozzle', align='left')# alpha=0.7)
ax2.axvline(x=mu4, color='black', label= 'Mean =' + str(round(mu4)) + ' km s$^{-1}$', linestyle='--', markersize=13)

ax3.hist(cut_valf_3p5kV_120mT, bins=bins, color='purple', edgecolor='black', label=r'124mT Nozzle', align='left')# alpha=0.7)
ax3.axvline(x=mu1, color='black', label='Mean =' + str(round(mu1)) + ' km s$^{-1}$', linestyle='--', markersize=13)

ax4.hist(cut_valf_3p5kV_133mT, bins=bins, color='purple', edgecolor='black', label=r'133mT Nozzle', align='left')# alpha=0.7)
ax4.axvline(x=mu2, color='black', label='Mean = ' + str(round(mu2)) + ' km s$^{-1}$', linestyle='--', markersize=13)

ax4.set_xlabel(r'$Alfv\acute{e}n~Speed~V_a~(km~s^{-1})$')
ax1.set_xlim(-100, 1400)
ax1.set_ylim(0, 110)
ax4.set_ylim(0, 30)
ax1.set_ylabel('Counts')
ax2.set_ylabel('Counts')
ax2.set_ylim(0, 50)
ax3.set_ylim(0, 40)
ax3.set_ylabel('Counts')
ax4.set_ylabel('Counts')
ax2.set_xlim(-100, 1400)
ax3.set_xlim(-100, 1400)
ax4.set_xlim(-100, 1400)

#ax0.legend(loc='best', frameon=False)
ax1.legend(loc='best', frameon=False)
ax2.legend(loc='best', frameon=False)
ax3.legend(loc='best', frameon=False)
ax4.legend(loc='best', frameon=False)


#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\AlfvenSpeedStudies\\AlfvenSpeed_NozzleComparison_Allposavgd_3p5kV.png', dpi=300)
#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\EPS_filesforpaper\\AlfvenSpeed_NozzleComparison_Allposavgd_3p5kV.eps', dpi=300)

#avg_alf_allprobes = np.mean([mu5,mu7,mu19,mu21,mu33,mu35])
#avg_stdAlf_allprobes = np.mean([sigma5,sigma7,sigma19,sigma21,sigma33,sigma35])
#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\2232022_3172022\\MagStructures_BulkVelocity\\AlfvenSpeedStudies\\Pos57_alfvenHistogram.png', dpi=300)


#for shot in np.arange(1,253):
#    data1r=data['pos19']['B']['r'][shot,:]
#    data1t=data['pos19']['B']['theta'][shot,:]
#    data1z=data['pos19']['B']['z'][shot,:]
#    data1mod=np.sqrt(data1r**2+data1t**2+data1z**2)
#    data1mod_max=np.max(data1mod)
#    data1mod_norm=data1mod/data1mod_max
#   data1mod_sm = sv(data1mod_norm, 11, 1)

#    Bmag_range = data1mod[start_time_index:end_time_index]


