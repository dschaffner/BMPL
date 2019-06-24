#sample_hdf5_data_loadin.py

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 23:31:43 2019

@author: dschaffner
"""

import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import numpy as np
import spectrum_wwind as spec
import indexfinderfuncs as iff

#######################################################################
# Directory style depends on Mac vs PC. For PC, use a double backslash.
# for a Mac, use a single forward slash.
### PC Style ###
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\04232019\\processed\\'
### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
#place the following file in the directory indicated above
datafilename='sample_2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'
#load hdf5 file
data = load_hdf5(data_directory_location+datafilename,verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6

#select time range for FFT (in us)
start_time = 50
end_time = 150
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,timeB_us)
end_time_index = iff.tindex_min(end_time,timeB_us)
#select shot to analyze
shot = 5

bt5 = data['pos5']['b']['theta'][5,:]
bz5 = data['pos5']['b']['z'][5,:]
mag5 = np.sqrt(bt5**2+bz5**2)


#Plot Details###############################
plt.rc('axes',linewidth=1.0)
plt.rc('xtick.major',width=1.0)
plt.rc('ytick.major',width=1.0)
plt.rc('xtick.minor',width=1.0)
plt.rc('ytick.minor',width=1.0)
plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=1,figsize=(4,4),dpi=200,facecolor='w',edgecolor='k')
left  = 0.21  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.15  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)   
#####################################  

### Plot arrow ####
arrowtime = 80.00#in us
timestep = iff.tindex_min(arrowtime,timeB_us)
plt.arrow(0, 0, bz5[timestep], bt5[timestep], head_width=100, head_length=100, fc='k', ec='k')
plt.arrow(0, 0, bz5[timestep], 0, head_width=10, head_length=10, fc='blue', ec='blue',alpha=0.3)
plt.arrow(bz5[timestep], 0, 0, bt5[timestep], head_width=10, head_length=10, fc='red', ec='red',alpha=0.3)

plt.xlabel(r'$B_{z}$ [G]',color='blue')
plt.ylabel(r'$B_{\theta}$ [G]',color='red')
plt.xlim(-3000,3000)
plt.ylim(-3000,3000)
timelabel=str(arrowtime)+r' $\mu$ s' 
plt.text(-2500,-2500,timelabel)
poslabel='Position 5'
plt.text(-2500,2500,poslabel)
    
        