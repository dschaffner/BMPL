# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
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


# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
datafilename='Dataset_01122022.h5'
data=load_hdf5(directory+datafilename,verbose=True)


time_s = data['time']['time_s']
timeB_s = time_s[1:]
time_us = data['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 2e-6#s
port_sep = 0.0254#m

numshots=94
direction_list=['r','t','z']
num_pairs = 3
directions = 3
tde_pairs=[['probe5','probe7']]#,['probe19','probe21'],['probe33','probe35']]

delaytimes = np.zeros([num_pairs, directions, numshots])
delayindex = np.zeros([num_pairs, directions, numshots])

for shot in np.arange(numshots):
    tde_pair_index=0
    for tde_pair in tde_pairs:
        direction_index=0
        for direction in direction_list:
            data1=data['mag_probe']['positions'][tde_pair[0]][direction]['b'][shot,:]
            data1_max=np.max(np.abs(data1))
            data1_norm=data1/data1_max
            
            data2=data['mag_probe']['positions'][tde_pair[1]][direction]['b'][shot,:]
            data2_max=np.max(np.abs(data2))
            data2_norm=data2/data2_max
            
            data1_filt=butter_bandpass_filter(data1_norm,lowcut,highcut,fs,order=9)
            data2_filt=butter_bandpass_filter(data2_norm,lowcut,highcut,fs,order=9)
            
            d1=data1_filt[start_time_index:end_time_index]
            d2=data2_filt[start_time_index:end_time_index]

            t=timeB_s[start_time_index:end_time_index]
            dt=timeB_s[1]-timeB_s[0]
            tau,corr=gc.get_corr(t,d2,d1,normalized=False)
            
            index_at_zero_time=iff.tindex_min(0,tau)
            indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
            #find time (in seconds) of max in correlation function within timerange limit
            delayindex[tde_pair_index,direction_index,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
            delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
            delaytimes[tde_pair_index,direction_index,shot]=delay
            direction_index+=1
        tde_pair_index+=1
            
velocities = (port_sep/delaytimes)/1000.0#km/s
mean_velocities = np.mean(velocities,axis=1)

plt.plot(velocities[0,0,:],'o')
plt.plot(velocities[0,1,:],'o')
plt.plot(velocities[0,2,:],'o')
plt.plot(mean_velocities[0,:],'x',color='red')


#########plot velocity distribution ###############
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=2,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[0,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
#plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)



#b5r=data['mag_probe']['positions']['probe5']['r']['b'][shot,:]
#b5r_max=np.max(np.abs(b5r))
#b5r_norm=b5r/b5r_max
#b7r=data['mag_probe']['positions']['probe7']['r']['b'][shot,:]
#b7r_max=np.max(np.abs(b7r))
#b7r_norm=b7r/b7r_max
#b5t=data['mag_probe']['positions']['probe5']['t']['b'][shot,:]
#b5t_max=np.max(np.abs(b5t))
#b5t_norm=b5t/b5t_max
#b7t=data['mag_probe']['positions']['probe7']['t']['b'][shot,:]
#b7t_max=np.max(np.abs(b7t))
#b7t_norm=b7t/b7t_max
#b5z=data['mag_probe']['positions']['probe5']['z']['b'][shot,:]
#b5z_max=np.max(np.abs(b5z))
#b5z_norm=b5z/b5z_max
#b7z=data['mag_probe']['positions']['probe7']['z']['b'][shot,:]
#b7z_max=np.max(np.abs(b7z))
#b7z_norm=b7z/b7z_max




#b5rfilt=butter_bandpass_filter(b5r_norm,lowcut,highcut,fs,order=9)
#b7rfilt=butter_bandpass_filter(b7r_norm,lowcut,highcut,fs,order=9)
#b5tfilt=butter_bandpass_filter(b5t_norm,lowcut,highcut,fs,order=9)
#b7tfilt=butter_bandpass_filter(b7t_norm,lowcut,highcut,fs,order=9)
#b5zfilt=butter_bandpass_filter(b5z_norm,lowcut,highcut,fs,order=9)
#b7zfilt=butter_bandpass_filter(b7z_norm,lowcut,highcut,fs,order=9)

#plt.figure(1)
#plt.clf()
#plt.plot(timeB_us,b5r)
#plt.plot(timeB_us,b5rfilt)
#plt.plot(timeB_us,b7rfilt)
# Plot the frequency response for a few different orders.
#plt.figure(1)
#plt.clf()
#for order in [3, 6, 9]:
#    sos = butter_bandpass(lowcut, highcut, fs, order=order)
#    w, h = sosfreqz(sos, worN=2000)
#    plt.semilogx((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

#d5r=b5rfilt[start_time_index:end_time_index]
#d7r=b7rfilt[start_time_index:end_time_index]
#d5t=b5tfilt[start_time_index:end_time_index]
#d7t=b7tfilt[start_time_index:end_time_index]
#d5z=b5zfilt[start_time_index:end_time_index]
#d7z=b7zfilt[start_time_index:end_time_index]

#t=timeB_us[start_time_index:end_time_index]
#tau57r,corr57r=gc.get_corr(t,d7r,d5r,normalized=False)
#tau57t,corr57t=gc.get_corr(t,d7t,d5t,normalized=False)
#tau57z,corr57z=gc.get_corr(t,d7z,d5z,normalized=False)

#plt.plot(tau57r,corr57r)
#plt.plot(tau57t,corr57t)
#plt.plot(tau57z,corr57z)


"""
plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2.5,markeredgewidth=0.0)

fig=plt.figure(num=1,figsize=(6,3),dpi=300,facecolor='w',edgecolor='k')

left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.05   # the bottom of the subplots of the figure
top = 0.97      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.05   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax=plt.axes([left,bottom,right-left,top-bottom])

for shot in np.arange(111,112):
    pico1=spio.loadmat(directory+'Pico1\\20220112-0001 ('+str(shot)+').mat')
    pico2=spio.loadmat(directory+'Pico2\\20220112-0001 ('+str(shot)+').mat')
    pico3=spio.loadmat(directory+'Pico3\\20220112-0001 ('+str(shot)+').mat')
    pico4=spio.loadmat(directory+'Pico4\\20220112-0001 ('+str(shot)+').mat')
    pico5=spio.loadmat(directory+'Pico5\\20220112-0001 ('+str(shot)+').mat')
    
    ax=plt.subplot(3,2,1)
    plt.plot(pico1['A'])
    plt.plot(pico1['B'])
    plt.plot(pico1['C'])
    plt.plot(pico1['D'])
    
    ax=plt.subplot(3,2,2)
    plt.plot(pico2['A'])
    plt.plot(pico2['B'])
    plt.plot(pico2['C'])
    plt.plot(pico2['D'])
    
    ax=plt.subplot(3,2,3)
    plt.plot(pico3['A'])
    plt.plot(pico3['B'])
    plt.plot(pico3['C'])
    plt.plot(pico3['D'])
    
    ax=plt.subplot(3,2,4)
    plt.plot(pico4['A'])
    plt.plot(pico4['B'])
    plt.plot(pico4['C'])
    plt.plot(pico4['D'])
    
    ax=plt.subplot(3,2,5)
    plt.plot(pico5['A'])
    plt.plot(pico5['B'])
    plt.plot(pico5['C'])
    plt.plot(pico5['D'])
    
    save_dir = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\'
    filename = 'Shot'+str(shot)+'.png'
    savefile = os.path.normpath(save_dir+filename)
    plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
    plt.clf()
   """ 