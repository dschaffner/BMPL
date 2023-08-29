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

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2023\\07132023\\'
maxshot=1
startshot=53
skipshots=[53,112,113,114,158,161,162,163]
numskipped=len(skipshots)
numshots=maxshot-(startshot-1)-numskipped

runfilename = 'Dataset_07132023_blockertest'


"""
high voltage - Pico1A
current - Pico1B
R2 - Pico 1C
T2 - Pico 1D

Z2 - Pico 2A
R1 - Pico 2B
T1 - Pico 2C
Z1 - Pico 2D

vf - Pico 3A

RTZ2 is closest to source, RTZ1 is furthest
all bdots and vf 50ohm terminated (mult signal by 2)


"""
term50=2.0


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
run_info.create_dataset('Discharge_Voltage (kV)',data=2.0e3)
run_info.create_dataset('Collection Dates',data='07132023')
run_info.create_dataset('Stuffing Delay (ms)',data=-1.5)
run_info.create_dataset('Gas Open (ms)',data=-2.0)
run_info.create_dataset('Gas Close (ms)',data=0.001)
startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=5.0


#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
#testpicoshot=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(1)+').mat')
#testpicoshot=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(10)+').txt',skiprows = 2, unpack = True)
testpicoshot=pd.read_csv(datadirectory+'\pico1\\20230713-0001 ('+str(6)+').csv',header=[0],skiprows=(1,2))


time=h5f.create_group("time")
tstart = testpicoshot['Time'][0]*1e-6#pico reads out in us, convert to seconds
time.create_dataset('tstart',data=tstart)
tinterval=np.abs(testpicoshot['Time'][0]*1e-6)
time.create_dataset('tinterval',data=tinterval)
numsamples=int(np.shape(testpicoshot['Time'])[0])
time.create_dataset('numsamples',data=numsamples)
time_s = tstart+(np.arange(numsamples)*tinterval)
time.create_dataset('time_s',data=time_s)
time.create_dataset('time_us',data=time_s*1e6)


testpicoshot=pd.read_csv(datadirectory+'\pico3\\20230713-0001 ('+str(6)+').csv',header=[0],skiprows=(1,2))
tstart_vf = testpicoshot['Time'][0]*1e-6#pico reads out in us, convert to seconds
time.create_dataset('tstart_vf',data=tstart_vf)
tinterval_vf=np.abs(testpicoshot['Time'][0]*1e-6)
time.create_dataset('tinterval_vf',data=tinterval_vf)
numsamples_vf=int(np.shape(testpicoshot['Time'])[0])
time.create_dataset('numsamples_vf',data=numsamples_vf)
time_s_vf = tstart_vf+(np.arange(numsamples_vf)*tinterval_vf)
time.create_dataset('time_s_vf',data=time_s_vf)
time.create_dataset('time_us_vf',data=time_s_vf*1e6)

###BEGIN Radial Scan READIN#####
num_rad_pos = 2
numshots=10
radposarray_in = np.array([0.0,2.0])#radial positions from center in inches
radposarray_cm = radposarray_in*2.54

Discharge_raw=np.zeros([num_rad_pos,numshots,numsamples])
HV_raw=np.zeros([num_rad_pos,numshots,numsamples])
Discharge=np.zeros([num_rad_pos,numshots,numsamples])
HV=np.zeros([num_rad_pos, numshots,numsamples])
Bdotr1 = np.zeros([num_rad_pos,numshots,numsamples])
Bdotr2 = np.zeros([num_rad_pos,numshots,numsamples])
Bdott1 = np.zeros([num_rad_pos,numshots,numsamples])
Bdott2 = np.zeros([num_rad_pos,numshots,numsamples])
Bdotz1 = np.zeros([num_rad_pos,numshots,numsamples])
Bdotz2 = np.zeros([num_rad_pos,numshots,numsamples])
Br1 = np.zeros([num_rad_pos,numshots,numsamples-1])
Br2 = np.zeros([num_rad_pos,numshots,numsamples-1])
Bt1 = np.zeros([num_rad_pos,numshots,numsamples-1])
Bt2 = np.zeros([num_rad_pos,numshots,numsamples-1])
Bz1 = np.zeros([num_rad_pos,numshots,numsamples-1])
Bz2 = np.zeros([num_rad_pos,numshots,numsamples-1])
vf = np.zeros([num_rad_pos,numshots,numsamples_vf])


### Radial Position Center ###
### Extract PICO 1 ###
radpos=0
savenumber=0
startshot=6
maxshot=15
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico1\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    Discharge_raw[radpos,savenumber,:]=np.array(pico['Channel A'])
    HV_raw[radpos,savenumber,:]=np.array(pico['Channel B'])
    Discharge[radpos,savenumber,:]=Discharge_raw[radpos,savenumber,:]*1000.0*2.0
    HV[radpos,savenumber,:]=HV_raw[radpos,savenumber,:]*1000.0
    
    Bdotr2[radpos,savenumber,:]=np.array(pico['Channel C'])*term50
    Bdott2[radpos,savenumber,:]=np.array(pico['Channel D'])*term50
    
    #filter Bdot
    Bdotr2[radpos,savenumber,:]=butter_bandpass_filter(Bdotr2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott2[radpos,savenumber,:]=butter_bandpass_filter(Bdott2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br2[radpos,savenumber,:]= sp.cumtrapz(Bdotr2[radpos,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt2[radpos,savenumber,:]= sp.cumtrapz(Bdott2[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO1 READIN #####
    
### Extract PICO 2 ###
savenumber=0
startshot=6
maxshot=15
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico2\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotz2[radpos,savenumber,:]=np.array(pico['Channel A'])*term50
    Bdotr1[radpos,savenumber,:]=np.array(pico['Channel B'])*term50
    Bdott1[radpos,savenumber,:]=np.array(pico['Channel C'])*term50
    Bdotz1[radpos,savenumber,:]=np.array(pico['Channel D'])*term50
    
    #filter Bdot
    Bdotz2[radpos,savenumber,:]=butter_bandpass_filter(Bdotz2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr1[radpos,savenumber,:]=butter_bandpass_filter(Bdotr1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott1[radpos,savenumber,:]=butter_bandpass_filter(Bdott1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz1[radpos,savenumber,:]=butter_bandpass_filter(Bdotz1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bz2[radpos,savenumber,:]= sp.cumtrapz(Bdotz2[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Br1[radpos,savenumber,:]= sp.cumtrapz(Bdotr1[radpos,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt1[radpos,savenumber,:]= sp.cumtrapz(Bdott1[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz1[radpos,savenumber,:]= sp.cumtrapz(Bdotz1[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####
    
### Extract PICO 3 ###
savenumber=0
startshot=6
maxshot=15
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico3\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    vf[radpos,savenumber,:]=np.array(pico['Channel A'])*200.0#diff vf amp setting
    #End of loop
    savenumber+=1
##END PICO3 READIN #####


### Radial Position 2in ###
### Extract PICO 1 ###
radpos=1
savenumber=0
startshot=24
maxshot=35
skipshots=[27,32]
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\pico1\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    Discharge_raw[radpos,savenumber,:]=np.array(pico['Channel A'])
    HV_raw[radpos,savenumber,:]=np.array(pico['Channel B'])
    Discharge[radpos,savenumber,:]=Discharge_raw[radpos,savenumber,:]*1000.0*2.0
    HV[radpos,savenumber,:]=HV_raw[radpos,savenumber,:]*1000.0
    
    Bdotr2[radpos,savenumber,:]=np.array(pico['Channel C'])*term50
    Bdott2[radpos,savenumber,:]=np.array(pico['Channel D'])*term50
    
    #filter Bdot
    Bdotr2[radpos,savenumber,:]=butter_bandpass_filter(Bdotr2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott2[radpos,savenumber,:]=butter_bandpass_filter(Bdott2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br2[radpos,savenumber,:]= sp.cumtrapz(Bdotr2[radpos,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt2[radpos,savenumber,:]= sp.cumtrapz(Bdott2[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO1 READIN #####
    
### Extract PICO 2 ###
savenumber=0
startshot=24
maxshot=35
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\pico2\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotz2[radpos,savenumber,:]=np.array(pico['Channel A'])*term50
    Bdotr1[radpos,savenumber,:]=np.array(pico['Channel B'])*term50
    Bdott1[radpos,savenumber,:]=np.array(pico['Channel C'])*term50
    Bdotz1[radpos,savenumber,:]=np.array(pico['Channel D'])*term50
    
    #filter Bdot
    Bdotz2[radpos,savenumber,:]=butter_bandpass_filter(Bdotz2[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr1[radpos,savenumber,:]=butter_bandpass_filter(Bdotr1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott1[radpos,savenumber,:]=butter_bandpass_filter(Bdott1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz1[radpos,savenumber,:]=butter_bandpass_filter(Bdotz1[radpos,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bz2[radpos,savenumber,:]= sp.cumtrapz(Bdotz2[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Br1[radpos,savenumber,:]= sp.cumtrapz(Bdotr1[radpos,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt1[radpos,savenumber,:]= sp.cumtrapz(Bdott1[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz1[radpos,savenumber,:]= sp.cumtrapz(Bdotz1[radpos,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####
    
    
### Extract PICO 3 ###
savenumber=0
startshot=24
maxshot=35
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    #pico=np.loadtxt(datadirectory+'\pico1\\20230627-0001 ('+str(shot)+').txt',skiprows = 2, unpack = True)
    pico=pd.read_csv(datadirectory+'\pico3\\20230713-0001 ('+str(shot)+').csv',header=[0],skiprows=(1,2))
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    vf[radpos,savenumber,:]=np.array(pico['Channel A'])*200.0#diff vf amp setting
    #End of loop
    savenumber+=1
##END PICO3 READIN #####    
    
    
    
    
    

    
 

#create magnetic probe group
mag = h5f.create_group("mag_probe")
r1=mag.create_group("r1")
r1.create_dataset("bdot",data=Bdotr1)
r1.create_dataset("b",data=Br1)
r2=mag.create_group("r2")
r2.create_dataset("bdot",data=Bdotr2)
r2.create_dataset("b",data=Br2)
t1=mag.create_group("t1")
t1.create_dataset("bdot",data=Bdott1)
t1.create_dataset("b",data=Bt1)
t2=mag.create_group("t2")
t2.create_dataset("bdot",data=Bdott2)
t2.create_dataset("b",data=Bt2)
z1=mag.create_group("z1")
z1.create_dataset("bdot",data=Bdotz1)
z1.create_dataset("b",data=Bz1)
z2=mag.create_group("z2")
z2.create_dataset("bdot",data=Bdotz2)
z2.create_dataset("b",data=Bz2)

#create floating potential probe group
vfloat=h5f.create_group("vf_probe")
vfloat.create_dataset("vf",data=vf)
                       
dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()
