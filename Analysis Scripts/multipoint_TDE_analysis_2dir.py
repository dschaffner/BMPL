#sample_hdf5_data_loadin.py

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 23:31:43 2019

@author: dschaffner
"""

import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import numpy as np
import indexfinderfuncs as iff
import get_corr as gc

#######################################################################
# Directory style depends on Mac vs PC. For PC, use a double backslash.
# for a Mac, use a single forward slash.
### PC Style ###
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\10142019\\processed\\'
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\07022019\\processed\\'
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\07052019\\processed\\'

### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
#place the following file in the directory indicated above
datafilename='sample_2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'
datafilename='2kV_0p5ms_1msGas_20shots_bdot_10142019.h5'
#datafilename='2kV_0p5ms_10msGas_20shots_bdot_10142019.h5'
#datafilename='2kV_oddpos9to21zr_2ms_stuffdelay_40_07022019.h5'
datafilename='2kV_oddpos9to21zr_0p5ms_stuffdelay_30shots_07052019.h5'
#load hdf5 file
data = load_hdf5(data_directory_location+datafilename,verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6

#select time range for FFT (in us)
start_time = 90
end_time = 140
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,timeB_us)
end_time_index = iff.tindex_min(end_time,timeB_us)
#select shots to analyze
first_shot = 1
last_shot = 16
numshots = 30
shot_range = [first_shot,last_shot]

delaytimes = np.zeros([6,2,numshots])#storage for number of delay times from cross correlation function
delayshifts = np.zeros([6,2,numshots])
timerange_limit = 5e-6 #limits the range for finding the peak to +/- 2 microseconds
dt=timeB_s[1]-timeB_s[0]
port_sep = 0.0254#m
direction_list=['r','z']
#loop over shots to read in data and compute cross correlation function
for shot in np.arange(numshots):
    for direction in np.arange(2):
        d1=data['pos11']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos9']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[0,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[0,direction,shot]=delay        

        d1=data['pos13']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos11']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[1,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[1,direction,shot]=delay  

        d1=data['pos15']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos13']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[2,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[2,direction,shot]=delay 
        
        d1=data['pos17']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos15']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[3,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[3,direction,shot]=delay 
        
        
        d1=data['pos19']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos17']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[4,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[4,direction,shot]=delay 


        d1=data['pos21']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        d2=data['pos19']['b'][direction_list[direction]][shot,start_time_index:end_time_index]
        t=timeB_s[start_time_index:end_time_index]
        tau,corr = gc.get_corr(t,d1,d2,normalized=False)
        index_at_zero_time=iff.tindex_min(0,tau)
        indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))
    
    #find time (in seconds) of max in correlation function within timerange limit
        delayshifts[5,direction,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
    #delay=tau[index_at_zero_time+np.argmax(np.abs(corr))]
    
        delaytimes[5,direction,shot]=delay 
        
velocities = (port_sep/delaytimes)/1000.0#km/s

velocities = velocities.reshape(6*2*numshots)
velocities = np.abs(velocities)
velocities_gt10=velocities[np.where(velocities>10.0)]
velocities_gt10=velocities_gt10[np.where(velocities_gt10<100.0)]
mean_vel=np.mean(velocities_gt10)
std_vel=np.std(velocities_gt10)
mean_vel_str=str(np.round(mean_vel,1))
std_vel_str=str(np.round(std_vel,1))
#########plot velocity distribution ###############
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=1,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(velocities_gt10, bins=10,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)

savefilename='11msGas_0p5msStuffing_9to21rz'+str(start_time)+'to'+str(end_time)+'us.png'
plt.savefig(savefilename,dpi=300,facecolor='w',edgecolor='k')
plt.clf()
plt.close()