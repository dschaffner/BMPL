# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 14:37:12 2022

@author: dschaffner
"""

import numpy as np
import pandas as pd
import matplotlib.pylab as plt
sheetdirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
spreadsheetfile = 'Separation Times by Shot.xlsx'
sheet = pd.read_excel(sheetdirectory+spreadsheetfile,header=1)
vels57 = np.array(sheet['57v'])
vels1921 = np.array(sheet['1921v'])
vels3335 = np.array(sheet['3335v'])
velsflags = np.array(sheet['flag'])
mean_vels57 = np.mean(vels57)
mean_vels1921 = np.mean(vels1921)
mean_vels3335 = np.mean(vels3335)

numshots=94
vdeltas_57to1921 = np.zeros([numshots])
vdeltas_1921to3335 = np.zeros([numshots])
for shot in np.arange(numshots):
    if velsflags[shot]==1: continue
    vdeltas_57to1921[shot]=vels1921[shot]-vels57[shot]
    vdeltas_1921to3335[shot]=vels3335[shot]-vels1921[shot]

mean_vdeltas_57to1921=np.mean(vdeltas_57to1921)
mean_vdeltas_1921to3335=np.mean(vdeltas_1921to3335)

#########plot velocity distribution ###############
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=571,figsize=(5,12),dpi=600,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.1  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(3,1,1)
plt.hist(vels57, bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.title('Probe 5 to 7')
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
plt.xlim(0,100)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
plt.ylim(0,25)
plt.vlines(mean_vels57,0,25,linestyle='dashed',color='black')
plt.text(mean_vels57+5,20,'Mean = '+str(round(mean_vels57,1)),fontsize=12)

ax1=plt.subplot(3,1,2)
plt.hist(vels1921, bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.title('Probe 19 to 21')
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
plt.xlim(0,100)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
plt.ylim(0,25)
plt.vlines(mean_vels1921,0,25,linestyle='dashed',color='black')
plt.text(mean_vels1921+5,20,'Mean = '+str(round(mean_vels1921,1)),fontsize=12)


ax1=plt.subplot(3,1,3)
plt.hist(vels3335, bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.title('Probe 33 to 35')
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
plt.xlim(0,100)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
plt.ylim(0,25)
plt.vlines(mean_vels3335,0,25,linestyle='dashed',color='black')
plt.text(mean_vels3335+5,20,'Mean = '+str(round(mean_vels3335,1)),fontsize=12)


save_dir = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\velocity\\'
savefilename = 'Velocity_Distribution_fromfluxarrival.png'
plt.savefig(save_dir+savefilename,dpi=600,facecolor='white',edgecolor='black')
plt.clf()
plt.close()

#########plot velocity deltas per shot ###############
fig=plt.figure(num=2,figsize=(5,8),dpi=600,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.1  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(2,1,1)
plt.hist(vdeltas_57to1921,bins=40,range=(-100,100))
plt.title('Delta V from Probes 5-7 to Probes 19-21')
plt.xlabel(r'$\Delta$ v [km/s]',fontsize=16)
plt.xlim(-100,100)
plt.ylabel('Count',fontsize=16)
plt.ylim(0,20)
plt.vlines(mean_vdeltas_57to1921,0,20,linestyle='dashed',color='black')
plt.text(mean_vdeltas_57to1921+10,15,'Mean = '+str(round(mean_vdeltas_57to1921,1)),fontsize=12)


ax1=plt.subplot(2,1,2)
plt.hist(vdeltas_1921to3335,bins=40,range=(-100,100))
plt.title('Delta V from Probes 19-21 to Probes 33-35')
plt.xlabel(r'$\Delta$ v [km/s]',fontsize=16)
plt.xlim(-100,100)
plt.ylabel('Count',fontsize=16)
plt.ylim(0,20)
plt.vlines(mean_vdeltas_1921to3335,0,20,linestyle='dashed',color='black')
plt.text(mean_vdeltas_1921to3335+10,15,'Mean = '+str(round(mean_vdeltas_1921to3335,1)),fontsize=12)


savefilename = 'Velocity_deltas_pershot_fromfluxarrival.png'
plt.savefig(save_dir+savefilename,dpi=600,facecolor='white',edgecolor='black')
plt.clf()
plt.close()
