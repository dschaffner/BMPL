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


import smooth as sm
def iter_smooth(array,loops=6,window_len=3):
    for l in np.arange(loops):
        array = sm.smooth(array,window_len=window_len)
    return array

# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\H5Files\\'
datafilename='2022_FullCatalogueofShots.h5'
data=load_hdf5(directory+datafilename,verbose=True)


time_s = np.array(data['07202022 Dataset']['Time']['time_s'])[0]
timeB_s = np.array(data['07202022 Dataset']['Time']['timeB_s'])[0]
time_us = np.array(data['07202022 Dataset']['Time']['time_us'])[0]
timeB_us = np.array(data['07202022 Dataset']['Time']['timeB_us'])[0]

Br5=np.array(data['07202022 Dataset']['pos5']['Bfilt']['r'])
Bt5=np.array(data['07202022 Dataset']['pos5']['Bfilt']['theta'])
Bz5=np.array(data['07202022 Dataset']['pos5']['Bfilt']['z'])
Bmod5 = np.sqrt(Br5**2+Bt5**2+Bz5**2)

Br7=np.array(data['07202022 Dataset']['pos7']['Bfilt']['r'])
Bt7=np.array(data['07202022 Dataset']['pos7']['Bfilt']['theta'])
Bz7=np.array(data['07202022 Dataset']['pos7']['Bfilt']['z'])
Bmod7 = np.sqrt(Br7**2+Bt7**2+Bz7**2)
#time_us = np.array(time_us)
#timeB_us = time_us[1:]

#select time range for FFT (in us)
start_time = 60
end_time = 160
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,timeB_us)
end_time_index = iff.tindex_min(end_time,timeB_us)

tau,corr=gc.get_corr(timeB_us,Bmod7[73,start_time_index:end_time_index],Bmod5[73,start_time_index:end_time_index],normalized=True)
#plt.plot(tau,corr)

plotrange_start = 122
plotrange_end = 136

#Plot Details###############################
plt.rc('axes',linewidth=0.5)
plt.rc('xtick.major',width=0.5)
plt.rc('ytick.major',width=0.5)
plt.rc('xtick.minor',width=0.5)
plt.rc('ytick.minor',width=0.5)
plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=0.5)
fig=plt.figure(num=1,figsize=(8,3),dpi=300,facecolor='w',edgecolor='k')
left  = 0.1  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.14  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,2,1)   

labelfont=6
legendfont=5
#####################################  

plt.plot(timeB_us,Bmod5[73,:],label='z=15.6cm')
plt.plot(timeB_us,Bmod7[73,:],label='z=18.2cm')
plt.ylabel('B[G]',fontsize=labelfont)
plt.ylim(0,2200)
plt.xlim(122,136)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
plt.xlabel('Time [us]',fontsize=labelfont+2)
plt.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.1, 0.9, '(a)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes)

ax2=plt.subplot(1,2,2)
plt.plot(tau,corr)
plt.ylabel('Norm. Correlation',fontsize=labelfont)
plt.ylim(0,0.15)
plt.xlim(-5,5)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
plt.xlabel('Delay Time [us]',fontsize=labelfont+2)
plt.text(0.1, 0.9, '(b)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)
#plt.legend(loc='best', frameon=False,fontsize=legendfont)

plt.savefig('C:\\Users\\dschaffner\\Documents\\GitHub\\BMPL\\InDevelopment\\velocitypaper_plots\\timeseries_and_corr.png', dpi=600)




"""
analysis_start_time = 20
analysis_end_time = 40
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m

numshots=94
direction_list=['r','t','z']
num_pairs = 3
directions = 3
tde_pairs=[['probe5','probe7'],['probe19','probe21'],['probe33','probe35']]

delaytimes = np.zeros([num_pairs, numshots])
delayindex = np.zeros([num_pairs, numshots])

plotshots=1
plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2.5,markeredgewidth=0.0)

fig=plt.figure(num=1,figsize=(6,3),dpi=300,facecolor='w',edgecolor='k')

left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.90      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.05   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
#ax=plt.axes([left,bottom,right-left,top-bottom])
ax=plt.subplot(1,1,1)

# Function to print mouse click event coordinates
def onclick(event):
    if event.key=='shift':
        print('Time is', event.xdata)

for shot in np.arange(0,numshots):
    tde_pair_index=0
    for tde_pair in tde_pairs:

        analysis_start_time = 60
        analysis_end_time = 160
        start_time_index = iff.tindex_min(analysis_start_time,time_us)
        end_time_index = iff.tindex_min(analysis_end_time,time_us)

        
        data1r=data['mag_probe']['positions'][tde_pair[0]]['r']['b'][shot,:]
        data1t=data['mag_probe']['positions'][tde_pair[0]]['t']['b'][shot,:]
        data1z=data['mag_probe']['positions'][tde_pair[0]]['z']['b'][shot,:]
        data1mod=np.sqrt(data1r**2+data1t**2+data1z**2)
        data1mod_max=np.max(data1mod[start_time_index:end_time_index])
        data1mod_norm=data1mod/data1mod_max
        #data1mod_norm=iter_smooth(data1mod_norm,loops=3,window_len=11)
        
        data2r=data['mag_probe']['positions'][tde_pair[1]]['r']['b'][shot,:]
        data2t=data['mag_probe']['positions'][tde_pair[1]]['t']['b'][shot,:]
        data2z=data['mag_probe']['positions'][tde_pair[1]]['z']['b'][shot,:]
        data2mod=np.sqrt(data2r**2+data2t**2+data2z**2)
        data2mod_max=np.max(data2mod[start_time_index:end_time_index])
        data2mod_norm=data2mod/data2mod_max
        #data2mod_norm=iter_smooth(data2mod_norm,loops=3,window_len=11)

        #data1mod_filt=butter_bandpass_filter(data1mod_norm,lowcut,highcut,fs,order=9)
        #data2mod_filt=butter_bandpass_filter(data2mod_norm,lowcut,highcut,fs,order=9)
        
        #correlation
        t=timeB_us[start_time_index:end_time_index]
        d1=data1mod_norm[start_time_index:end_time_index]
        d2=data2mod_norm[start_time_index:end_time_index]
        tau,corr=gc.get_corr(t,d2,d1,normalized=False)
        
        
        #use filtered
        #d1=data1_filt[start_time_index:end_time_index]
        #d2=data2_filt[start_time_index:end_time_index]
        if plotshots:
            
            plt.plot(np.array(timeB_us),data1mod_norm)
            plt.plot(np.array(timeB_us),data2mod_norm)
            plt.title(tde_pair[0]+tde_pair[1]+' Shot '+str(shot))
            plt.xlim(0.0,40.0)
            plt.ylim(0,0.3)
            
            # Bind the button_press_event with the onclick() method
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()    
            plt.ginput(10,timeout=0)
            plt.clf()
            
            plt.plot(tau,corr)
            plt.title(tde_pair[0]+tde_pair[1]+' Shot '+str(shot))
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()    
            plt.ginput(20,timeout=0)
            plt.clf()
            #savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\timeseries\\probepairs\\'+tde_pair[0]+tde_pair[1]+'\\mod\\'
            #savefilename=tde_pair[0]+tde_pair[1]+'_shot_'+str(shot+1).zfill(2)+'mod.png'
            #savefile = savedirectory+savefilename
            #plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
            #plt.clf()
        #use unfiltered
        #d1=data1_norm[start_time_index:end_time_index]
        #d2=data2_norm[start_time_index:end_time_index]

        #t=timeB_s[start_time_index:end_time_index]
        #dt=timeB_s[1]-timeB_s[0]
        #tau,corr=gc.get_corr(t,data2mod_norm,data1mod_norm,normalized=False)
        
        #index_at_zero_time=iff.tindex_min(0,tau)
        #indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))

        #find time (in seconds) of max in correlation function within timerange limit
        #delayindex[tde_pair_index,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        #delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
        #delaytimes[tde_pair_index,shot]=delay
    tde_pair_index+=1
            
#velocities = (port_sep/delaytimes)/1000.0#km/s
#mean_velocities = np.mean(velocities,axis=1)

#np.savez(directory+'mean_vels_01122022_60t160_50to500kHzfilt',mean_velocities=mean_velocities)
"""
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