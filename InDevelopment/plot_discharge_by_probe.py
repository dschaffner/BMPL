# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
"""

import matplotlib.pylab as plt
import numpy as np
from load_hdf5 import load_hdf5
import spectrum_wwind as spec
import indexfinderfuncs as iff

#Put in your own directory below
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
timerange_limit = 3e-6#s
port_sep = 0.0254#m


#Plot voltage discharge
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
ax=plt.subplot(1,1,1)



shot = 1 #input shot number to plot (94 shots total)
data1=data['discharge']['dis_V'][shot]
plt.plot(np.array(time_us),-data1)
plt.ylim(0,2000)
plt.title('Voltage Shot '+str(shot+1))
#directory information for saving the file Change to your own directory!!!
savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\timeseries\\voltage\\20t35us\\'
savefilename='shot_'+str(shot+1).zfill(2)+'_voltage.png'
savefile = savedirectory+savefilename
plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
#plt.clf()
#plt.close()