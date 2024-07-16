#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""
import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import scipy.integrate as sp
import os
import pandas as pd

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\03262024\\'
maxshot=1
startshot=77
skipshots=[1000]
numskipped=len(skipshots)
numshots=maxshot-(startshot-1)-numskipped

runfilename = 'Dataset_03262024'


"""
high voltage - Pico1A
current - Pico1B
R2 - Pico 1C
T2 - Pico 1D

Z2 - Pico 2A
R1 - Pico 2B
T1 - Pico 2C
Z1 - Pico 2D



RTZ2 is closest to source, RTZ1 is furthest

06272023 data was probe at center
06282023 data was probe at 1in, 2in, and 3in away from center radially
"""

from scipy.signal import butter, sosfiltfilt, sosfreqz

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y
    
# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 50e6

import h5py
h5f = h5py.File(datadirectory+runfilename+'.h5','a')

#double stalked bdot probe
probe_dia = 0.00158755#m (1/16'' probe)
hole_sep = 0.001016#m (1/16''probe)
probe_info=h5f.create_group("probe_info")
probe_info.create_dataset('probe_stalk_diameter',data=probe_dia)
probe_info.create_dataset('probe_hole_separation',data=hole_sep)
probe_info.create_dataset('rloop_area',data=np.pi*(probe_dia/2)**2)
probe_info.create_dataset('tloop_area',data=probe_dia*hole_sep)
probe_info.create_dataset('zloop_area',data=probe_dia*hole_sep)
#probe_info.create_dataset('radial_location',data='center')
probe_info.create_dataset('voltage_probe_factor',data=1000.0)#in Volts
probe_info.create_dataset('current_probe_factor',data=1000.0*2)# in Amps
rloop_area = np.pi*((probe_dia/2)**2)
tloop_area = probe_dia*hole_sep
zloop_area = probe_dia*hole_sep

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data=2.2e3)
run_info.create_dataset('Collection Dates',data='06282023')
run_info.create_dataset('Stuffing Delay (ms)',data=-1.5)
run_info.create_dataset('Gas Open (ms)',data=-2.0)
run_info.create_dataset('Gas Close (ms)',data=0.001)
startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=5.0


#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
#testpicoshot=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(1)+').mat')
#testpicoshot=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(10)+').txt',skiprows = 2, unpack = True)
testpicoshot=pd.read_csv(datadirectory+'\pico1\\20240326-0001 ('+str(6)+').csv',header=[0],skiprows=[1])


time=h5f.create_group("time")
tstart = testpicoshot['Time'][0]*1e-6#pico reads out in us, convert to seconds
time.create_dataset('tstart',data=tstart)
tinterval=np.abs((testpicoshot['Time'][1]-testpicoshot['Time'][0])*1e-6)
time.create_dataset('tinterval',data=tinterval)
numsamples=int(np.shape(testpicoshot['Time'])[0])
time.create_dataset('numsamples',data=numsamples)
time_s = tstart+(np.arange(numsamples)*tinterval)
time.create_dataset('time_s',data=time_s)
time.create_dataset('time_us',data=time_s*1e6)

#wirecurrent time
testpicoshot=pd.read_csv(datadirectory+'\pico6\\20240326-0001 ('+str(6)+').csv',header=[0],skiprows=[1])


time_wc=h5f.create_group("time_wc")
tstart_wc = testpicoshot['Time'][0]*1e-3#pico reads out in us, convert to seconds
time_wc.create_dataset('tstart',data=tstart_wc)
tinterval_wc=np.abs(testpicoshot['Time'][0]*1e-3)
time_wc.create_dataset('tinterval',data=tinterval_wc)
numsamples_wc=int(np.shape(testpicoshot['Time'])[0])
time_wc.create_dataset('numsamples',data=numsamples_wc)
time_s_wc = tstart_wc+(np.arange(numsamples_wc)*tinterval_wc)
time_wc.create_dataset('time_s',data=time_s_wc)
time_wc.create_dataset('time_ms',data=time_s_wc*1e3)




###BEGIN Radial Scan READIN#####
num_rad_pos = 1 #by 0.5in from center
numshots=72 #per radial position
numprobes = 5
#radposarray_in = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])#radial positions from center in inches
#radposarray_cm = radposarray_in*2.54

Discharge_raw=np.zeros([num_rad_pos,numshots,numsamples])
HV_raw=np.zeros([num_rad_pos,numshots,numsamples])
Discharge=np.zeros([num_rad_pos,numshots,numsamples])
HV=np.zeros([num_rad_pos,numshots,numsamples])
wirecurrent=np.zeros([num_rad_pos,numshots,numsamples_wc])
Bdotr = np.zeros([num_rad_pos,numprobes,numshots,numsamples])
Bdott = np.zeros([num_rad_pos,numprobes,numshots,numsamples])
Bdotz = np.zeros([num_rad_pos,numprobes,numshots,numsamples])
Br = np.zeros([num_rad_pos,numprobes,numshots,numsamples-1])
Bt = np.zeros([num_rad_pos,numprobes,numshots,numsamples-1])
Bz = np.zeros([num_rad_pos,numprobes,numshots,numsamples-1])
term50=2.0


### Radial Position Center ###
radpos=0
startshot=6
maxshot=77

### Extract PICO 1 ###
savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico1\\20240326-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    Discharge_raw[radpos,savenumber,:]=np.array(pico['Channel A'],dtype=float)
    HV_raw[radpos,savenumber,:]=np.array(pico['Channel B'],dtype=float)
    Discharge[radpos,savenumber,:]=Discharge_raw[radpos,savenumber,:]*1000.0*2.0
    HV[radpos,savenumber,:]=HV_raw[radpos,savenumber,:]*1000.0
    
    Bdotr[radpos,0,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50
    Bdott[radpos,0,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50
    
    #filter Bdot
    Bdotr[radpos,0,savenumber,:]=butter_bandpass_filter(Bdotr[radpos,0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[radpos,0,savenumber,:]=butter_bandpass_filter(Bdott[radpos,0,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br[radpos,0,savenumber,:]= sp.cumtrapz(Bdotr[radpos,0,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[radpos,0,savenumber,:]= sp.cumtrapz(Bdott[radpos,0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO1 READIN ##### 
        
### Extract PICO 2 ###
savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico2\\20240326-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotz[radpos,0,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50
    Bdotr[radpos,1,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50
    Bdott[radpos,1,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50
    Bdotz[radpos,1,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50
    
    #filter Bdot
    Bdotz[radpos,0,savenumber,:]=butter_bandpass_filter(Bdotz[radpos,0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr[radpos,1,savenumber,:]=butter_bandpass_filter(Bdotr[radpos,1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[radpos,1,savenumber,:]=butter_bandpass_filter(Bdott[radpos,1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[radpos,1,savenumber,:]=butter_bandpass_filter(Bdotz[radpos,1,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bz[radpos,0,savenumber,:]= sp.cumtrapz(Bdotz[radpos,0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Br[radpos,1,savenumber,:]= sp.cumtrapz(Bdotr[radpos,1,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[radpos,1,savenumber,:]= sp.cumtrapz(Bdott[radpos,1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[radpos,1,savenumber,:]= sp.cumtrapz(Bdotz[radpos,1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####
    
### Extract PICO 3 ###
savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico3\\20240326-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotr[radpos,2,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50
    Bdott[radpos,2,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50
    Bdotz[radpos,2,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50
    Bdotr[radpos,3,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50
    
    #filter Bdot
    Bdotr[radpos,2,savenumber,:]=butter_bandpass_filter(Bdotr[radpos,2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[radpos,2,savenumber,:]=butter_bandpass_filter(Bdott[radpos,2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[radpos,2,savenumber,:]=butter_bandpass_filter(Bdotz[radpos,2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr[radpos,3,savenumber,:]=butter_bandpass_filter(Bdotr[radpos,3,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br[radpos,2,savenumber,:]= sp.cumtrapz(Bdotr[radpos,2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bt[radpos,2,savenumber,:]= sp.cumtrapz(Bdott[radpos,2,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bz[radpos,2,savenumber,:]= sp.cumtrapz(Bdotz[radpos,2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Br[radpos,3,savenumber,:]= sp.cumtrapz(Bdotr[radpos,3,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO3 READIN #####

### Extract PICO 4 ###
savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico4\\20240326-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdott[radpos,3,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50
    Bdotz[radpos,3,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50
    Bdotr[radpos,4,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50
    Bdott[radpos,4,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50
    
    #filter Bdot
    Bdott[radpos,3,savenumber,:]=butter_bandpass_filter(Bdott[radpos,3,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[radpos,3,savenumber,:]=butter_bandpass_filter(Bdotz[radpos,3,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[radpos,4,savenumber,:]=butter_bandpass_filter(Bdotr[radpos,4,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr[radpos,4,savenumber,:]=butter_bandpass_filter(Bdott[radpos,4,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bt[radpos,3,savenumber,:]= sp.cumtrapz(Bdott[radpos,3,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[radpos,3,savenumber,:]= sp.cumtrapz(Bdotz[radpos,3,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Br[radpos,4,savenumber,:]= sp.cumtrapz(Bdotr[radpos,4,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bt[radpos,4,savenumber,:]= sp.cumtrapz(Bdott[radpos,4,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO4 READIN #####
    
### Extract PICO 6 ###
savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico6\\20240326-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    
    wirecurrent[radpos,savenumber,:]=np.array(pico['Channel A'],dtype=float)*1000.0

    #End of loop
    savenumber+=1
##END PICO5 READIN #####

















#create magnetic probe group
mag = h5f.create_group("mag_probe")
r=mag.create_group("r")
r.create_dataset("bdot",data=Bdotr)
r.create_dataset("b",data=Br)
t=mag.create_group("t")
t.create_dataset("bdot",data=Bdott)
t.create_dataset("b",data=Bt)
z=mag.create_group("z")
z.create_dataset("bdot",data=Bdotz)
z.create_dataset("b",data=Bz)

#create wirecurrent probe group
wc=h5f.create_group("wirecurrent")
wc.create_dataset("wirecurrent",data=wirecurrent)
                       
#create discharage group
dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()
