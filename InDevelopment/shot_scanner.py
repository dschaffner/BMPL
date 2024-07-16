# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
"""
import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import os

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\03172022\\'
#shots that should be skipped
#43, 71, 72, 73, 76, 80, 82, 87
#89, 90, 91, 92, 93, 94, 100
#102, 106

#test=spio.loadmat(directory+'Pico1\\20220112-0001 (2).mat')

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

for shot in np.arange(65,66):
    pico1=spio.loadmat(directory+'Pico1\\20220317-0001 ('+str(shot)+').mat')
    pico2=spio.loadmat(directory+'Pico2\\20220317-0001 ('+str(shot)+').mat')
    pico3=spio.loadmat(directory+'Pico3\\20220317-0001 ('+str(shot)+').mat')
    pico4=spio.loadmat(directory+'Pico4\\20220317-0001 ('+str(shot)+').mat')
    pico5=spio.loadmat(directory+'Pico5\\20220317-0001 ('+str(shot)+').mat')
    
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
    #plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
    #plt.clf()
    
    