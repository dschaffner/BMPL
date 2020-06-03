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
#data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\04232019\\processed\\'

### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
#place the following file in the directory indicated above
datafilename='2kV_oddpos9to21zr_2ms_stuffdelay_40_07022019.h5'
#datafilename='2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'

#load hdf5 file
data = load_hdf5(data_directory_location+datafilename,verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6
time_s = data['time']['time_s'][:]

poslist = ['pos9','pos11','pos13','pos15','pos17','pos19','pos21']

#select time range for FFT (in us)
start_time = 50
end_time = 150
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,timeB_us)
end_time_index = iff.tindex_min(end_time,timeB_us)
#select shots to analyze
first_shot = 1
last_shot = 40
numshots = (last_shot-first_shot)+1
shot_range = [first_shot,last_shot]

Bmag = np.zeros([7,numshots,len(data['pos9']['b']['z'][0,:])])
for pos in np.arange(len(poslist)):
        x=data[poslist[pos]]['b']['r'][:]
        y=data[poslist[pos]]['b']['z'][:]
        Bmag[pos,:] = np.sqrt(x**2+y**2)

#determine FFT size and generate an output array
fsize=int((data['pos9']['b']['z'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot = np.zeros([7,2,fsize])
avebspec_direct = np.zeros([7,2,fsize])
avebmagspec = np.zeros([7,fsize])



#loop over shots to read in data and compute FFT

#poslist = ['pos1','pos3','pos5','pos7','pos9','pos11','pos13']
dirlist = ['r','z']
#for shot in np.arange(first_shot-1,last_shot):
#    for pos in np.arange(len(poslist)):
#        for direction in np.arange(len(dirlist)):
#            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['b'][dirlist[direction]][shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
#            avebspec[pos,direction,:]=avebspec[pos,direction]+pwr
            
for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        for direction in np.arange(len(dirlist)):
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['bdot'][dirlist[direction]][shot,start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hanning')
            avebspec_frombdot[pos,direction,:]=avebspec_frombdot[pos,direction]+(pwr/(f*f))
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['b'][dirlist[direction]][shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
            avebspec_direct[pos,direction,:]=avebspec_direct[pos,direction]+(pwr)
            
for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(Bmag[pos,shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
        avebmagspec[pos,:]=avebmagspec[pos,:]+pwr

allpos_avebspec_frombdot_2=np.sum(avebspec_frombdot,axis=0)
allpos_avebspec_frombdot_2=np.sum(allpos_avebspec_frombdot_2,axis=0)

allpos_avebmagspec_2=np.sum(avebmagspec,axis=0)
        

data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\07052019\\processed\\'
datafilename='2kV_oddpos9to21zr_0p5ms_stuffdelay_30shots_07052019.h5'
#load hdf5 file
data = load_hdf5(data_directory_location+datafilename,verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6
time_s = data['time']['time_s'][:]
#select shots to analyze
first_shot = 1
last_shot = 30
numshots = (last_shot-first_shot)+1
shot_range = [first_shot,last_shot]

Bmag = np.zeros([7,numshots,len(data['pos9']['b']['z'][0,:])])
for pos in np.arange(len(poslist)):
        x=data[poslist[pos]]['b']['r'][:]
        y=data[poslist[pos]]['b']['z'][:]
        Bmag[pos,:] = np.sqrt(x**2+y**2)
#determine FFT size and generate an output array
fsize=int((data['pos9']['b']['z'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot = np.zeros([7,2,fsize])
avebspec_direct = np.zeros([7,2,fsize])
avebmagspec = np.zeros([7,fsize])            
for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        for direction in np.arange(len(dirlist)):
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['bdot'][dirlist[direction]][shot,start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hanning')
            avebspec_frombdot[pos,direction,:]=avebspec_frombdot[pos,direction]+(pwr/(f*f))
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['b'][dirlist[direction]][shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
            avebspec_direct[pos,direction,:]=avebspec_direct[pos,direction]+(pwr)
for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(Bmag[pos,shot,start_time_index:end_time_index],timeB_s[start_time_index:end_time_index],window='hanning')
        avebmagspec[pos,:]=avebmagspec[pos,:]+pwr

allpos_avebspec_frombdot_p5=np.sum(avebspec_frombdot,axis=0)
allpos_avebspec_frombdot_p5=np.sum(allpos_avebspec_frombdot_p5,axis=0)

allpos_avebmagspec_p5=np.sum(avebmagspec,axis=0)

#########plot loglog FFT###############
        
#Plot Details###############################
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=1,figsize=(4,3),dpi=300,facecolor='w',edgecolor='k')
left  = 0.2  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.90      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)   
#####################################  
    


#compute indices from freq range
start_f = 1e5
end_f = 7e5
start_f_index = iff.tindex_min(start_f,f)
end_f_index = iff.tindex_min(end_f,f)

dlog=np.log10(allpos_avebmagspec_p5[start_f_index:end_f_index])
flog=np.log10(f[start_f_index:end_f_index])
A1=np.array([flog,np.ones(len(flog))])
w1=np.linalg.lstsq(A1.T,dlog)[0]
slope=np.round(w1[0],3)
func = f**(slope)
scale_data=allpos_avebmagspec_p5[100]
scale_fit=func[100]
ratio=scale_data/scale_fit

#plot FFT#
#for pos in np.arange(len(poslist)):
plt.loglog(f[1:],allpos_avebmagspec_p5[1:],label='0.5ms')
plt.loglog(f[start_f_index:end_f_index],allpos_avebmagspec_p5[start_f_index:end_f_index])
#plt.loglog(f,ratio*func)
#plt.loglog(f[1:],avebspec[0,1,1:])
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlabel('Frequency [Hz]',fontsize=9)
plt.ylabel('B Magnitude Power (arb)',fontsize=9)

print ('0.5ms slope is',slope)


dlog=np.log10(allpos_avebmagspec_2[start_f_index:end_f_index])
flog=np.log10(f[start_f_index:end_f_index])
A1=np.array([flog,np.ones(len(flog))])
w1=np.linalg.lstsq(A1.T,dlog)[0]
slope=np.round(w1[0],3)
func = f**(slope)
scale_data=allpos_avebmagspec_2[100]
scale_fit=func[100]
ratio=scale_data/scale_fit

#plot FFT#
#for pos in np.arange(len(poslist)):
plt.loglog(f[1:],allpos_avebmagspec_2[1:],label='2.0ms')
plt.loglog(f[start_f_index:end_f_index],allpos_avebmagspec_2[start_f_index:end_f_index])
#plt.loglog(f,ratio*func)
#plt.loglog(f[1:],avebspec[0,1,1:])
plt.title('2.0ms vs 0.5ms - 7/2/19 vs 7/5/19 - Allshots,positions',fontsize=9)

print ('2.0ms slope is',slope)

plt.legend(loc='lower left',fontsize=8)

plt.savefig('2vsp5.png',dpi=600,facecolor='w',edgecolor='k')