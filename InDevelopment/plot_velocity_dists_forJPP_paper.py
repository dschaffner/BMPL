# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:44:37 2023

@author: dschaffner
"""

import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import scipy.integrate as sp
import os
import pandas as pd

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\Processed Data\\'
data1=np.loadtxt(datadirectory+'Full3172022_ShotsandVelocities_P5P7.txt',skiprows = 1, unpack = True)
data2=np.loadtxt(datadirectory+'Full6202022_ShotsandVelocities_P5P7.txt',skiprows = 1, unpack = True)
data3=np.loadtxt(datadirectory+'Full6202022_ShotsandVelocities_P19P21.txt',skiprows = 1, unpack = True)
data4=np.loadtxt(datadirectory+'Full7152022_ShotsandVelocities_P5P7.txt',skiprows = 1, unpack = True)
data5=np.loadtxt(datadirectory+'Full7152022_ShotsandVelocities_P19P21.txt',skiprows = 1, unpack = True)
data6=np.loadtxt(datadirectory+'Full7202022_ShotsandVelocities_P5P7.txt',skiprows = 1, unpack = True)
data7=np.loadtxt(datadirectory+'Full7202022_ShotsandVelocities_P19P21.txt',skiprows = 1, unpack = True)
data8=np.loadtxt(datadirectory+'Full9232022_ShotsandVelocities_P5P7.txt',skiprows = 1, unpack = True)
data9=np.loadtxt(datadirectory+'Full9232022_ShotsandVelocities_P6P8.txt',skiprows = 1, unpack = True)


#resort the data into arrays by voltage, position, nozzle
#from data1
vels_2p5kV_pos57_110mT = data1[2][0:72]
vels_2p5kV_pos57_55mT = data1[2][72:]

#from data2 and data4
vels_3p0kV_pos57_110mT = data2[2][0:18]
vels_3p5kV_pos57_110mT = np.concatenate((data2[2][18:],data4[2][:]))

#from data3
vels_3p0kV_pos1921_110mT = data3[2][:]

#from data5
vels_3p5kV_pos1921_110mT = data5[2][:]

#from data6
vels_3p5kV_pos57_124mT = data6[2][0:53]
vels_3p5kV_pos57_133mT = data6[2][53:]

#from data7
vels_3p5kV_pos1921_124mT = data7[2][0:53]
vels_3p5kV_pos1921_133mT = data7[2][53:]

#from data8
vels_3p5kV_pos57_0mT = data8[2][90:]

#from data9
vels_3p5kV_pos68_0mT = data9[2][90:]

cutoff_vel = 251

#compute means
mean_vels_2p5kV_pos57_110mT = np.mean(vels_2p5kV_pos57_110mT,where=vels_2p5kV_pos57_110mT<=cutoff_vel)
mean_vels_2p5kV_pos57_55mT = np.mean(vels_2p5kV_pos57_55mT, where=vels_2p5kV_pos57_55mT<=cutoff_vel)
mean_vels_3p0kV_pos57_110mT = np.mean(vels_3p0kV_pos57_110mT, where=vels_3p0kV_pos57_110mT<=cutoff_vel)
mean_vels_3p5kV_pos57_110mT = np.mean(vels_3p5kV_pos57_110mT, where=vels_3p5kV_pos57_110mT<=cutoff_vel)
mean_vels_3p0kV_pos1921_110mT = np.mean(vels_3p0kV_pos1921_110mT, where=vels_3p0kV_pos1921_110mT<=cutoff_vel)
mean_vels_3p5kV_pos1921_110mT = np.mean(vels_3p5kV_pos1921_110mT, where=vels_3p5kV_pos1921_110mT<=cutoff_vel)
mean_vels_3p5kV_pos57_124mT = np.mean(vels_3p5kV_pos57_124mT, where=vels_3p5kV_pos57_124mT<=cutoff_vel)
mean_vels_3p5kV_pos57_133mT= np.mean(vels_3p5kV_pos57_133mT, where=vels_3p5kV_pos57_133mT<=cutoff_vel)
mean_vels_3p5kV_pos1921_124mT = np.mean(vels_3p5kV_pos1921_124mT, where=vels_3p5kV_pos1921_124mT<=cutoff_vel)
mean_vels_3p5kV_pos1921_133mT = np.mean(vels_3p5kV_pos1921_133mT, where=vels_3p5kV_pos1921_133mT<=cutoff_vel)
mean_vels_3p5kV_pos57_0mT = np.mean(vels_3p5kV_pos57_0mT, where=vels_3p5kV_pos57_0mT<=cutoff_vel)
mean_vels_3p5kV_pos68_0mT = np.mean(vels_3p5kV_pos68_0mT, where=vels_3p5kV_pos68_0mT<=cutoff_vel)

#compute std
std_vels_2p5kV_pos57_110mT = np.std(vels_2p5kV_pos57_110mT,where=vels_2p5kV_pos57_110mT<=cutoff_vel)
std_vels_2p5kV_pos57_55mT = np.std(vels_2p5kV_pos57_55mT, where=vels_2p5kV_pos57_55mT<=cutoff_vel)
std_vels_3p0kV_pos57_110mT = np.std(vels_3p0kV_pos57_110mT, where=vels_3p0kV_pos57_110mT<=cutoff_vel)
std_vels_3p5kV_pos57_110mT = np.std(vels_3p5kV_pos57_110mT, where=vels_3p5kV_pos57_110mT<=cutoff_vel)
std_vels_3p0kV_pos1921_110mT = np.std(vels_3p0kV_pos1921_110mT, where=vels_3p0kV_pos1921_110mT<=cutoff_vel)
std_vels_3p5kV_pos1921_110mT = np.std(vels_3p5kV_pos1921_110mT, where=vels_3p5kV_pos1921_110mT<=cutoff_vel)
std_vels_3p5kV_pos57_124mT = np.std(vels_3p5kV_pos57_124mT, where=vels_3p5kV_pos57_124mT<=cutoff_vel)
std_vels_3p5kV_pos57_133mT= np.std(vels_3p5kV_pos57_133mT, where=vels_3p5kV_pos57_133mT<=cutoff_vel)
std_vels_3p5kV_pos1921_124mT = np.std(vels_3p5kV_pos1921_124mT, where=vels_3p5kV_pos1921_124mT<=cutoff_vel)
std_vels_3p5kV_pos1921_133mT = np.std(vels_3p5kV_pos1921_133mT, where=vels_3p5kV_pos1921_133mT<=cutoff_vel)
std_vels_3p5kV_pos57_0mT = np.std(vels_3p5kV_pos57_0mT, where=vels_3p5kV_pos57_0mT<=cutoff_vel)
std_vels_3p5kV_pos68_0mT = np.std(vels_3p5kV_pos68_0mT, where=vels_3p5kV_pos68_0mT<=cutoff_vel)
#print means

print('2p5kV_pos57_110mT mean =', mean_vels_2p5kV_pos57_110mT)
print('2p5kV_pos57_55mT mean = ',mean_vels_2p5kV_pos57_55mT)
print('3p0kV_pos57_110mT mean =',mean_vels_3p0kV_pos57_110mT)
print('3p5kV_pos57_110mT mean = ',mean_vels_3p5kV_pos57_110mT)
print('3p0kV_pos1921_110mT mean =',mean_vels_3p0kV_pos1921_110mT) 
print('3p5kV_pos1921_110mT mean =',mean_vels_3p5kV_pos1921_110mT) 
print('3p5kV_pos57_124mT mean =',mean_vels_3p5kV_pos57_124mT) 
print('3p5kV_pos57_133mT mean =',mean_vels_3p5kV_pos57_133mT)
print('3p5kV_pos1921_124mT mean =',mean_vels_3p5kV_pos1921_124mT)
print('3p5kV_pos1921_133mT mean =',mean_vels_3p5kV_pos1921_133mT)
print('3p5kV_pos57_0mT mean = ',mean_vels_3p5kV_pos57_0mT)
print('3p5kV_pos68_0mT mean = ',mean_vels_3p5kV_pos68_0mT)




bins1 = np.array(range(20, 250, 15))

#Plot Details###############################
plt.rc('axes',linewidth=0.5)
plt.rc('xtick.major',width=0.5)
plt.rc('ytick.major',width=0.5)
plt.rc('xtick.minor',width=0.5)
plt.rc('ytick.minor',width=0.5)
plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=0.5)
fig=plt.figure(num=1,figsize=(4,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.1  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.1  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(4,1,1)   
#####################################  


labelfont=6
legendfont=5
histheight=19

ax1.hist(vels_3p5kV_pos57_0mT, bins=bins1, label='0G Nozzle', align='left', edgecolor='black', color='teal',linewidth=0.5)
ax1.axvline(x=mean_vels_3p5kV_pos57_0mT, color='black', label='Mean = '+str(np.round(mean_vels_3p5kV_pos57_0mT,1)) + ' km/s', linestyle='--', markersize='13')
ax1.set_ylabel('Shots',fontsize=labelfont)
ax1.set_ylim(0, histheight)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
ax1.set_xticklabels([])
ax1.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.05, 0.9, '(a)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes,fontsize=6)

ax2=plt.subplot(4,1,2)
ax2.hist(vels_3p5kV_pos57_110mT, bins=bins1, label='1110G Nozzle', align='left', edgecolor='black', color='teal',linewidth=0.5)
ax2.axvline(x=mean_vels_3p5kV_pos57_110mT, color='black', label='Mean = '+str(np.round(mean_vels_3p5kV_pos57_110mT,1)) + ' km/s', linestyle='--', markersize='13')
ax2.set_ylabel('Shots',fontsize=labelfont)
ax2.set_ylim(0, histheight)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
ax2.set_xticklabels([])
ax2.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.05, 0.9, '(b)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes,fontsize=6)

ax3=plt.subplot(4,1,3)
ax3.hist(vels_3p5kV_pos57_124mT, bins=bins1, label='1240G Nozzle', align='left', edgecolor='black', color='teal',linewidth=0.5)
ax3.axvline(x=mean_vels_3p5kV_pos57_124mT, color='black', label='Mean = '+str(np.round(mean_vels_3p5kV_pos57_124mT,1)) + ' km/s', linestyle='--', markersize='13')
ax3.set_ylabel('Shots',fontsize=labelfont)
ax3.set_ylim(0, histheight)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
ax3.set_xticklabels([])
ax3.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.05, 0.9, '(c)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax3.transAxes,fontsize=6)

ax4=plt.subplot(4,1,4)
ax4.hist(vels_3p5kV_pos57_133mT, bins=bins1, label='1330G Nozzle', align='left', edgecolor='black', color='teal',linewidth=0.5)
ax4.axvline(x=mean_vels_3p5kV_pos57_133mT, color='black', label='Mean = '+str(np.round(mean_vels_3p5kV_pos57_133mT,1)) + ' km/s', linestyle='--', markersize='13')
ax4.set_ylabel('Shots',fontsize=labelfont)
ax4.set_ylim(0, histheight)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
ax4.set_xlabel(r'Plume Velocity $[km/s]$',fontsize=labelfont+2)
ax4.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.05, 0.9, '(d)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax4.transAxes,fontsize=6)

plt.savefig('C:\\Users\\dschaffner\\Documents\\GitHub\\BMPL\\InDevelopment\\velocitypaper_plots\\nozzlescan.png', dpi=600)





scan_3p5kV_pos57_G = np.array([0,1110,1240,1330])
scan_3p5kV_pos57_vels = np.array([mean_vels_3p5kV_pos57_0mT,mean_vels_3p5kV_pos57_110mT,
                                      mean_vels_3p5kV_pos57_124mT,mean_vels_3p5kV_pos57_133mT])
scan_3p5kV_pos57_vels_err = np.array([std_vels_3p5kV_pos57_0mT,std_vels_3p5kV_pos57_110mT,
                                      std_vels_3p5kV_pos57_124mT,std_vels_3p5kV_pos57_133mT])

scan_pos57_1100G_kV = np.array([2.5,3.0,3.5])
scan_pos57_1100G_vels = np.array([mean_vels_2p5kV_pos57_110mT,mean_vels_3p0kV_pos57_110mT,
                                      mean_vels_3p5kV_pos57_110mT])
scan_pos57_1100G_vels_err = np.array([std_vels_2p5kV_pos57_110mT,std_vels_3p0kV_pos57_110mT,
                                      std_vels_3p5kV_pos57_110mT])


#Plot Details###############################
plt.rc('axes',linewidth=0.5)
plt.rc('xtick.major',width=0.5)
plt.rc('ytick.major',width=0.5)
plt.rc('xtick.minor',width=0.5)
plt.rc('ytick.minor',width=0.5)
plt.rc('lines',markersize=4.0,markeredgewidth=0.0,linewidth=0.5)
fig=plt.figure(num=2,figsize=(6,3),dpi=300,facecolor='w',edgecolor='k')
left  = 0.1  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.15  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.1  # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,2,1)   
#####################################  

plt.plot(scan_3p5kV_pos57_G,scan_3p5kV_pos57_vels,'--o',linewidth=0.25,label='Dis. Voltage = 3.5kV')
plt.errorbar(scan_3p5kV_pos57_G,scan_3p5kV_pos57_vels,yerr=scan_3p5kV_pos57_vels_err,capsize=3, capthick=1.5, color='black', linewidth=1, ls='none')
plt.ylabel('Ave. Plume Velocity [km/s]',fontsize=labelfont)
plt.ylim(20,200)
plt.yticks(fontsize=labelfont)
plt.xticks(fontsize=labelfont)
plt.xlabel('Nozzle Field [G]',fontsize=labelfont+2)
plt.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.1, 0.9, '(a)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes)

ax2=plt.subplot(1,2,2)
plt.plot(scan_pos57_1100G_kV,scan_pos57_1100G_vels,'--o',linewidth=0.25,label='Nozzle Field = 1110G')
plt.errorbar(scan_pos57_1100G_kV,scan_pos57_1100G_vels,yerr=scan_pos57_1100G_vels_err,capsize=3, capthick=1.5, color='black', linewidth=1, ls='none')
ax2.set_yticklabels([])
plt.ylim(20, 200)
plt.xticks(fontsize=labelfont)
plt.xlabel('Discharge Voltage [kV]',fontsize=labelfont+2)
plt.legend(loc='best', frameon=False,fontsize=legendfont)
plt.text(0.1, 0.9, '(b)',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)

plt.savefig('C:\\Users\\dschaffner\\Documents\\GitHub\\BMPL\\InDevelopment\\velocitypaper_plots\\scatterplots.png', dpi=600)
