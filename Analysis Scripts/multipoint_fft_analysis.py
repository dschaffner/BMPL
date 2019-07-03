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
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\07022019\\processed\\'
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\04232019\\processed\\'
### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
#place the following file in the directory indicated above
datafilename='2kV_oddpos9to21zr_2ms_stuffdelay_40_07022019.h5'
datafilename='2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'
#load hdf5 file
data = load_hdf5(data_directory_location+datafilename,verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6

x=data['pos1']['b']['r'][:]
y=data['pos1']['b']['z'][:]
Bmag = np.sqrt(x**2+y**2)

#select time range for FFT (in us)
start_time = 50
end_time = 150
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,timeB_us)
end_time_index = iff.tindex_min(end_time,timeB_us)
#select shots to analyze
first_shot = 1
last_shot = 17
numshots = (last_shot-first_shot)+1
shot_range = [first_shot,last_shot]

#determine FFT size and generate an output array
fsize=int((data['pos9']['b']['z'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec = np.zeros([7,2,fsize])
avebmagspec = np.zeros([fsize])



#loop over shots to read in data and compute FFT
poslist = ['pos9','pos11','pos13','pos15','pos17','pos19','pos21']
poslist = ['pos1','pos3','pos5','pos7','pos9','pos11','pos13']
dirlist = ['r','z']
#for shot in np.arange(first_shot-1,last_shot):
#    for pos in np.arange(len(poslist)):
#        for direction in np.arange(len(dirlist)):
#            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['b'][dirlist[direction]][shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
#            avebspec[pos,direction,:]=avebspec[pos,direction]+pwr
            
for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        for direction in np.arange(len(dirlist)):
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['bdot'][dirlist[direction]][shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
            avebspec[pos,direction,:]=avebspec[pos,direction]+(pwr/(f*f))
            
for shot in np.arange(first_shot-1,last_shot):
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(Bmag[shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
    avebmagspec=avebmagspec+pwr
        
#########plot loglog FFT###############
        
#Plot Details###############################
plt.rc('axes',linewidth=1.0)
plt.rc('xtick.major',width=1.0)
plt.rc('ytick.major',width=1.0)
plt.rc('xtick.minor',width=1.0)
plt.rc('ytick.minor',width=1.0)
plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=1,figsize=(4,3),dpi=300,facecolor='w',edgecolor='k')
left  = 0.2  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)   
#####################################  
    


#compute indices from freq range
start_f = 1e5
end_f = 1e6
start_f_index = iff.tindex_min(start_f,f)
end_f_index = iff.tindex_min(end_f,f)

dlog=np.log10(avebmagspec[start_f_index:end_f_index])
flog=np.log10(f[start_f_index:end_f_index])
A1=np.array([flog,np.ones(len(flog))])
w1=np.linalg.lstsq(A1.T,dlog)[0]
slope=np.round(w1[0],3)
func = f**(slope)
scale_data=avebmagspec[100]
scale_fit=func[100]
ratio=scale_data/scale_fit

#plot FFT#
#for pos in np.arange(len(poslist)):
plt.loglog(f[1:],avebmagspec[1:])
plt.loglog(f[start_f_index:end_f_index],avebmagspec[start_f_index:end_f_index])
plt.loglog(f,ratio*func)
#plt.loglog(f[1:],avebspec[0,1,1:])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power (arb)')

print ('B slope is',slope)