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


def loadnpzfile(filename, supress=False):

    savefile = filename
    savefile = os.path.normpath(savefile)

    file = np.load(savefile)
    if not supress:
        print('Arrays loaded: ')
        for arr in file.files:
            print(arr)
    return file



# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
datafilename='Dataset_01122022.h5'
data=load_hdf5(directory+datafilename,verbose=True)
vels=loadnpzfile(directory+'mean_vels_01122022.npz')
velocity=np.abs(vels['mean_velocities'])*1000.0#in m/s


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

numshots=94
direction_list=['r','t','z']
probelist=['probe5','probe7','probe19','probe21','probe33','probe35']
tde_pairs=[['probe5','probe7'],['probe19','probe21'],['probe33','probe35']]
directions = len(direction_list)
numprobes = len(probelist)


#determine FFT size and generate an output array
fsize=int((data['mag_probe']['positions']['probe5']['r']['bdot'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot = np.zeros([numprobes,directions,fsize])
avebspec_direct = np.zeros([numprobes,directions,fsize])
#avebmagspec = np.zeros([numprobes,fsize])
power_frombdot = np.zeros([numshots,numprobes,directions,fsize])
wavenum = np.zeros([numshots,numprobes,directions,fsize])

#loop over shots to read in data and compute FFT
for shot in np.arange(numshots):
    tde_pair_index=0
    velocity_index=0
    for probepair in tde_pairs:
        print ('on probe pair',probepair, 'with index,', tde_pair_index)
        for direction_index, direction in enumerate(direction_list):
            #probepair first
            data1=data['mag_probe']['positions'][probepair[0]][direction]['bdot'][shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hanning')
            power_frombdot[shot,tde_pair_index,direction_index,:]=pwr/(f*f)
            wavenum[shot,tde_pair_index,direction_index,:]=f/velocity[velocity_index,shot]
            #probepair second
            data1=data['mag_probe']['positions'][probepair[1]][direction]['bdot'][shot,:]
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hanning')
            power_frombdot[shot,tde_pair_index+1,direction_index,:]=pwr/(f*f)
            wavenum[shot,tde_pair_index+1,direction_index,:]=f/velocity[velocity_index,shot]
        velocity_index+=1
        tde_pair_index+=2

plotshot = 0
if plotshot:
    fig=plt.figure(num=1,figsize=(5,3.5),dpi=300,facecolor='w',edgecolor='k')
    left  = 0.15  # the left side of the subplots of the figure
    right = 0.94    # the right side of the subplots of the figure
    bottom = 0.2  # the bottom of the subplots of the figure
    top = 0.96      # the top of the subplots of the figure
    wspace = 0.2   # the amount of width reserved for blank space between subplots
    hspace = 0.1   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    ax1=plt.subplot(1,1,1)
#interpolate wavenumber traces to sum wavenumber spectra together
probelist=['probe5','probe7','probe19','probe21','probe33','probe35']
newwavenum = np.arange(0.3,400,0.03)
power_interp = np.zeros([numshots,numprobes,directions,len(newwavenum)])
for shot in np.arange(numshots):
    velocity_index=0
    if shot == 43: continue
    if shot == 22: continue
    for probe_index, probe in enumerate(probelist):
        for direction_index, direction in enumerate(direction_list):
            x=wavenum[shot,probe_index,direction_index,:]
            y=power_frombdot[shot,probe_index,direction_index,:]
            f=interp1d(x,y)
            power_interp[shot,probe_index,direction_index,:]=f(newwavenum)
            if plotshot:
                plt.loglog(newwavenum,power_interp[shot,probe_index,direction_index,:])
                plt.title(probe+' Shot '+str(shot)+'_'+direction)
                plt.xlim(0.1,1000)
                plt.ylim(1e-35,1e-19)
                if probe_index == 0 or probe_index == 1:
                    velocitylabel = str(np.round(velocity[0,shot]/1000.0,decimals=0))
                if probe_index == 2 or probe_index == 3:
                    velocitylabel = str(np.round(velocity[1,shot]/1000.0,decimals=0))
                if probe_index == 4 or probe_index == 5:
                    velocitylabel = str(np.round(velocity[2,shot]/1000.0,decimals=0))
                plt.text(0.95,0.95,velocitylabel,fontsize=16,horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
                savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\Spectra\\wavenum_by_shot\\'
                savefilename=probe+'_shot_'+str(shot)+'_'+direction+'.png'
                savefile = savedirectory+savefilename
                plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
                plt.clf()
        velocity_index+=1
meanpower_interp = np.mean(power_interp,axis=0)
"""

plt.figure(57)
plt.plot(velocities[0,0,:],'o')
plt.plot(velocities[0,1,:],'o')
plt.plot(velocities[0,2,:],'o')
plt.plot(mean_velocities[0,:],'x',color='red')

plt.figure(1921)
plt.plot(velocities[1,0,:],'o')
plt.plot(velocities[1,1,:],'o')
plt.plot(velocities[1,2,:],'o')
plt.plot(mean_velocities[1,:],'x',color='red')

plt.figure(3335)
plt.plot(velocities[2,0,:],'o')
plt.plot(velocities[2,1,:],'o')
plt.plot(velocities[2,2,:],'o')
plt.plot(mean_velocities[2,:],'x',color='red')



#########plot velocity distribution ###############
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=571,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
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

fig=plt.figure(num=19211,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[1,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
#plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)


fig=plt.figure(num=33351,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[2,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
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

