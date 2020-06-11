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
import get_corr as gc
from loadnpzfile import loadnpzfile
import os


#######################################################################
# Directory style depends on Mac vs PC. For PC, use a double backslash.
# for a Mac, use a single forward slash.
### PC Style ###
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2020\\03102020\\'
### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
#place the following file in the directory indicated above
#datafilename='Dataset_03102020_Bfilt_1e4.npz'
#datafilename='Dataset_03102020_Bfilt_5e4.npz'
#datafilename='Dataset_03102020_Bfilt_7e4.npz'
#datafilename='Dataset_03102020_Bfilt_8e4.npz'
datafilename='Dataset_03102020_Bfilt_1e5.npz'
#datafilename='Dataset_03102020_Bfilt_1e6.npz'

plottitle='3-10-20 Correlation B(t) 75-175ms 100kHz HPF'

#load npz file
datafile=loadnpzfile(data_directory_location+datafilename)
time_ms = datafile['time_ms']
time_s = datafile['time_s']
bfilt1 = datafile['bfilt1']

#select time range for FFT (in us)
start_time = 75
end_time = 175
time_range = [start_time,end_time]
#compute indices from time range
start_time_index = iff.tindex_min(start_time,time_ms)
end_time_index = iff.tindex_min(end_time,time_ms)
first_shot = 1
last_shot = 25
numshots = (last_shot-first_shot)+1
shot_range = [first_shot,last_shot]

"""
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
"""

plt.rc('lines',markersize=3.0,markeredgewidth=0.0,linewidth=1.0)
#Correlation Using get_corr.py
norm_corr_ave = np.zeros([8])
direction = 0
for shot in np.arange(25):
    tau,corr0,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,0],normalized=False,optimize=False,mode='valid')
    tau,corr1,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,1],normalized=False,optimize=False,mode='valid')
    tau,corr2,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,2],normalized=False,optimize=False,mode='valid')
    tau,corr3,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,3],normalized=False,optimize=False,mode='valid')
    tau,corr4,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,4],normalized=False,optimize=False,mode='valid')
    tau,corr5,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,5],normalized=False,optimize=False,mode='valid')
    tau,corr6,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,6],normalized=False,optimize=False,mode='valid')
    tau,corr7,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,7],normalized=False,optimize=False,mode='valid')
    norm_corr = np.absolute(np.array([1.0,corr1[0]/corr0[0],corr2[0]/corr0[0],corr3[0]/corr0[0],corr4[0]/corr0[0],corr5[0]/corr0[0],corr6[0]/corr0[0],corr7[0]/corr0[0]]))
    space = (np.arange(8)*2.6)
    norm_corr_ave = norm_corr_ave+norm_corr
    #plt.plot(space[1:],norm_corr[1:],'o')
norm_corr_ave_r=norm_corr_ave/numshots

norm_corr_ave = np.zeros([8])
direction = 1
for shot in np.arange(25):
    tau,corr0,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,0],normalized=False,optimize=False,mode='valid')
    tau,corr1,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,1],normalized=False,optimize=False,mode='valid')
    tau,corr2,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,2],normalized=False,optimize=False,mode='valid')
    tau,corr3,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,3],normalized=False,optimize=False,mode='valid')
    tau,corr4,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,4],normalized=False,optimize=False,mode='valid')
    tau,corr5,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,5],normalized=False,optimize=False,mode='valid')
    tau,corr6,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,6],normalized=False,optimize=False,mode='valid')
    tau,corr7,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,7],normalized=False,optimize=False,mode='valid')
    norm_corr = np.absolute(np.array([1.0,corr1[0]/corr0[0],corr2[0]/corr0[0],corr3[0]/corr0[0],corr4[0]/corr0[0],corr5[0]/corr0[0],corr6[0]/corr0[0],corr7[0]/corr0[0]]))
    space = (np.arange(8)*2.6)
    norm_corr_ave = norm_corr_ave+norm_corr
    #plt.plot(space[1:],norm_corr[1:],'o')
norm_corr_ave_t=norm_corr_ave/numshots

norm_corr_ave = np.zeros([8])
direction = 2
for shot in np.arange(25):
    tau,corr0,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,0],normalized=False,optimize=False,mode='valid')
    tau,corr1,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,1],normalized=False,optimize=False,mode='valid')
    tau,corr2,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,2],normalized=False,optimize=False,mode='valid')
    tau,corr3,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,3],normalized=False,optimize=False,mode='valid')
    tau,corr4,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,4],normalized=False,optimize=False,mode='valid')
    tau,corr5,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,5],normalized=False,optimize=False,mode='valid')
    tau,corr6,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,6],normalized=False,optimize=False,mode='valid')
    tau,corr7,mean1,mean2 = gc.get_corr_wmean(time_s[start_time_index:end_time_index],bfilt1[direction,start_time_index:end_time_index,shot,0],bfilt1[direction,start_time_index:end_time_index,shot,7],normalized=False,optimize=False,mode='valid')
    norm_corr = np.absolute(np.array([1.0,corr1[0]/corr0[0],corr2[0]/corr0[0],corr3[0]/corr0[0],corr4[0]/corr0[0],corr5[0]/corr0[0],corr6[0]/corr0[0],corr7[0]/corr0[0]]))
    space = (np.arange(8)*2.6)
    norm_corr_ave = norm_corr_ave+norm_corr
    #plt.plot(space[1:],norm_corr[1:],'o')
norm_corr_ave_z=norm_corr_ave/numshots

norm_corr_ave_all = (norm_corr_ave_r+norm_corr_ave_t+norm_corr_ave_z)/3.0

plt.rc('lines',markersize=10.0,markeredgewidth=0.0,linewidth=1.0)
plt.plot(space,norm_corr_ave_r,'D',color='blue',label=r'$r$')
plt.plot(space,norm_corr_ave_t,'D',color='green',label=r'$\theta$')
plt.plot(space,norm_corr_ave_z,'D',color='red',label=r'$z$')
plt.plot(space,norm_corr_ave_all,'o',color='purple',label='Ave all')
plt.xlim(0,20)
plt.ylim(0,1)
plt.xlabel('Distance from Probe 1 [cm]',fontsize=9)
plt.ylabel('Shot Ave. Correlation',fontsize=9)
plt.title(plottitle,fontsize=12)
leg=plt.legend(loc='upper right',fontsize=12,frameon=False,handlelength=5)


savefilename=plottitle
#savefilename='CR6_1t_Rg_Full_DeltaCount_tdyn1totdyn2_wranges.png'
savefile = os.path.normpath(data_directory_location+savefilename)
#plt.savefig(savefile,dpi=200,facecolor='w',edgecolor='k')
"""

from scipy.optimize import curve_fit
def func(x,a):
    return 1-((x*2)/(2*a*2))

#norm_corr = np.absolute(np.array([1.0,corr1[0]/corr0[0],corr2[0]/corr0[0],corr3[0]/corr0[0],corr4[0]/corr0[0],corr5[0]/corr0[0],corr6[0]/corr0[0]]))
#space = (np.arange(7)*1.3)
rich_extrap = np.zeros([6])
points = np.arange(6)+2

popt,pcov=curve_fit(func,space,norm_corr_ave)
r_array = np.arange(80)*0.1
#bigR=1-(r_array*2/(2*2.6*2))
plt.plot(r_array,func(r_array,popt[0]))
print('6 point fit:',popt[0])
rich_extrap[-1]=popt[0]

popt,pcov=curve_fit(func,space[:-1],norm_corr_ave[:-1])
plt.plot(r_array,func(r_array,popt[0]))
print('5 point fit:',popt[0])
rich_extrap[-2]=popt[0]

popt,pcov=curve_fit(func,space[:-2],norm_corr_ave[:-2])
plt.plot(r_array,func(r_array,popt[0]))
print('4 point fit:',popt[0])
rich_extrap[-3]=popt[0]

popt,pcov=curve_fit(func,space[:-3],norm_corr_ave[:-3])
plt.plot(r_array,func(r_array,popt[0]))
print('3 point fit:',popt[0])
rich_extrap[-4]=popt[0]

popt,pcov=curve_fit(func,space[:-4],norm_corr_ave[:-4])
plt.plot(r_array,func(r_array,popt[0]))
print('2 point fit:',popt[0])
rich_extrap[-5]=popt[0]

popt,pcov=curve_fit(func,space[:-5],norm_corr_ave[:-5])
plt.plot(r_array,func(r_array,popt[0]))
print('1 point fit:',popt[0])
rich_extrap[-6]=popt[0]

plt.figure(2)
plt.plot(points,rich_extrap)


def lin_func(x,a,b):
    return (x*a)+b
popt,pcov=curve_fit(lin_func,points,rich_extrap)
lin_array=np.arange(60)*0.1
plt.plot(lin_array,lin_func(lin_array,popt[0],popt[1]))
print('$\lambda_{0}=$',popt[1])

"""
#Correlation Using manual calculation
#data1_mean = np.mean(data['pos9']['b']['z'][0,start_time_index:end_time_index])
#data2_mean = np.mean(data['pos11']['b']['z'][0,start_time_index:end_time_index])
#data1_sub = data['pos9']['b']['z'][0,start_time_index:end_time_index]-data1_mean
#data2_sub = data['pos11']['b']['z'][0,start_time_index:end_time_index]-data2_mean
#corr1_man = np.sum(data1_sub*data2_sub)
#corr1_man = corr1_man/(timeB_s[end_time_index]-timeB_s[start_time_index])

#data1_mean = np.mean(data['pos9']['b']['z'][0,start_time_index:end_time_index])
#data2_mean = np.mean(data['pos13']['b']['z'][0,start_time_index:end_time_index])
#data1_sub = data['pos9']['b']['z'][0,start_time_index:end_time_index]-data1_mean
#data2_sub = data['pos13']['b']['z'][0,start_time_index:end_time_index]-data2_mean
#corr2_man = np.sum(data1_sub*data2_sub)

"""
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

allpos_avebspec_frombdot=np.sum(avebspec_frombdot,axis=0)
allpos_avebspec_frombdot=np.sum(allpos_avebspec_frombdot,axis=0)

allpos_avebmagspec=np.sum(avebmagspec,axis=0)
        
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

dlog=np.log10(allpos_avebmagspec[start_f_index:end_f_index])
flog=np.log10(f[start_f_index:end_f_index])
A1=np.array([flog,np.ones(len(flog))])
w1=np.linalg.lstsq(A1.T,dlog)[0]
slope=np.round(w1[0],3)
func = f**(slope)
scale_data=allpos_avebmagspec[100]
scale_fit=func[100]
ratio=scale_data/scale_fit

#plot FFT#
#for pos in np.arange(len(poslist)):
plt.loglog(f[1:],allpos_avebmagspec[1:])
plt.loglog(f[start_f_index:end_f_index],allpos_avebmagspec[start_f_index:end_f_index])
plt.loglog(f,ratio*func)
#plt.loglog(f[1:],avebspec[0,1,1:])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power (arb)')

print ('B slope is',slope)
"""