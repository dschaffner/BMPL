#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""
import scipy.io as spio
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import spectrum_wwind as spec

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\bdotcalibration\\'

#load frequencies
#load velocities
sheetdirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\bdotcalibration\\'
spreadsheetfile = 'Calibration Test.xlsx'
sheet = pd.read_excel(sheetdirectory+spreadsheetfile,header=0)
freqs = np.array(sheet['Frequency (kHz)'])*1000.0#in Hz

#storage arrays
current = np.zeros([33,12500])
bdot = np.zeros([33,12500])

#make time array
bcal=spio.loadmat(datadirectory+'bdotcal_01.mat')
time=np.arange(12500)*bcal['Tinterval'][0]

#load fluctuation data
for n in np.arange(1,34):
    bcal=spio.loadmat(datadirectory+'bdotcal_'+str(n).zfill(2)+'.mat')
    current[n-1,:]=np.array(bcal['A'][:,0])
    bdot[n-1,:]=np.array(bcal['average_B_'][:,0])

neginfs = np.isneginf(bdot)
bdot[np.where(neginfs)]=0.0
posinfs = np.isposinf(bdot)
bdot[np.where(posinfs)]=0.0

inputpower=np.zeros([33])
bdotpower=np.zeros([33])
    
#plt.plot(time,bdot[0,:])
for n in np.arange(0,33):
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(current[n,:],time,window='hanning')
    inputpower[n]=np.max(pwr)
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(bdot[n,:],time,window='hanning')
    print('For ',n, ' max power is ', np.max(pwr[1:])/freqs[n]**2)
    bdotpower[n]=np.max(pwr[1:])/((freqs[n])**2)
    #plt.loglog(f,pwr/(f**2),'o')


plt.figure(3)    
#plt.loglog(freqs,inputpower,'o-',color='blue')
plt.loglog(freqs,bdotpower,'o-',color='red')

plt.figure(333)
#plt.semilogx(freqs,inputpower/inputpower[0],'o-',color='blue')
#plt.semilogx(freqs,bdotpower/bdotpower[0],'o-',color='red')
plt.semilogx(freqs,(bdotpower/bdotpower[0])/(inputpower/inputpower[0]),'o-',color='black')
plt.xlabel('Frequency')
plt.ylabel('Fraction of Max Power')

save_dir = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\bdotcalibration\\'
filename = 'bdot_calibration_04042022.png'
plt.savefig(save_dir+filename,dpi=300,facecolor='w',edgecolor='k')
plt.clf()
plt.close() 

bdotcalibrationarray=(bdotpower/bdotpower[0])/(inputpower/inputpower[0])
np.savez(datadirectory+'bdot_calibration_array_04042022.npz',
         bdotcal=bdotcalibrationarray,freqs=freqs)