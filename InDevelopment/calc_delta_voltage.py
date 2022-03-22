# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:09:50 2022

@author: dschaffner
"""
import scipy.io as spio
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import os
from load_hdf5 import load_hdf5
import spectrum_wwind as spec
import indexfinderfuncs as iff
import get_corr as gc
import scipy.integrate as sp
from scipy.interpolate import interp1d
import plasmapy as plasmapy
from astropy import units as u
from astropy import constants as const

from scipy.signal import butter, sosfiltfilt, sosfreqz




def loadnpzfile(filename, supress=False):

    savefile = filename
    savefile = os.path.normpath(savefile)

    file = np.load(savefile)
    if not supress:
        print('Arrays loaded: ')
        for arr in file.files:
            print(arr)
    return file



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



#load velocities
sheetdirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
spreadsheetfile = 'Separation Times by Shot.xlsx'
sheet = pd.read_excel(sheetdirectory+spreadsheetfile,header=1)
vels57 = np.array(sheet['57v'])*1000.0#m/s
vels1921 = np.array(sheet['1921v'])*1000.0#m/s
vels3335 = np.array(sheet['3335v'])*1000.0#m/s
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


numshots=94
direction_list=['r','t','z']
probelist=['probe5','probe7','probe19','probe21','probe33','probe35']
tde_pairs=[['probe5','probe7'],['probe19','probe21'],['probe33','probe35']]
directions = len(direction_list)
numprobes = len(probelist)


stuffing = 2.0
delta_V = np.zeros([94])
time_s = np.array(time_s)

for shot in np.arange(numshots):
    if velsflags[shot] == 1: continue
    data1=-data['discharge']['dis_V'][shot]
    #v_offset=np.mean(data1[iff.tindex_min(10.0,timeB_us):iff.tindex_min(15.0,timeB_us)])
    #data1=data1-v_offset
    #integrate 10 to 50us
    analysis_start_time = 10
    analysis_end_time = 30
    start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
    end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
    delta_V[shot] = np.std(data1[start_time_index:end_time_index])
    
    
        
        #plt.plot(timeB_us,b_energy_density)
        #plt.ylim(-700,7000)
        #plt.title(probe+' Shot '+str(shot+1))
        #savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\timeseries\\'+probe+'\\bmod\\'
        #savefilename=probe+'_shot_'+str(shot+1).zfill(2)+'_EnDens.png'
        #savefile = savedirectory+savefilename
        #plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
        #plt.clf()

plt.plot(delta_V)

pickprobe=0
pickdeltas = vdeltas_57to1921#to3335
#pickdeltas = vdeltas_1921to3335
#pickdeltas = vels57
fig=plt.figure(num=2,figsize=(5,3.5),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)

#y=np.where(pickdeltas!=0)
#pickdeltas=pickdeltas[y]
#x=np.where(inj_hel>0)
#inj_hel=inj_hel[y]

plt.scatter(delta_V,pickdeltas)
plt.title(r'$\Delta V$ 10 to 50us')
fig1corr = np.corrcoef(delta_V,pickdeltas)
print(fig1corr)
#plt.plot(inj_hel,inj_hel*fig1corr[0][1],'o',color='red')

plt.figure(333)
plt.plot(pickdeltas,'o-')