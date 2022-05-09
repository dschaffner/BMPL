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

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\03172022\\'
maxshot=253
startshot=181
skipshots=[181,182,183]
numskipped=len(skipshots)
numshots=maxshot-(startshot-1)-numskipped

runfilename = 'Dataset_03172022_halfouter_inner'

import h5py
h5f = h5py.File(datadirectory+runfilename+'.h5','a')

probe_dia = 0.00158755#m (1/16'' probe)
hole_sep = 0.001016#m (1/16''probe)
probe_info=h5f.create_group("probe_info")
probe_info.create_dataset('probe_stalk_diameter',data=probe_dia)
probe_info.create_dataset('probe_hole_separation',data=hole_sep)
probe_info.create_dataset('rloop_area',data=np.pi*(probe_dia/2)**2)
probe_info.create_dataset('tloop_area',data=probe_dia*hole_sep)
probe_info.create_dataset('zloop_area',data=probe_dia*hole_sep)
probe_info.create_dataset('radial_location',data='center')
probe_info.create_dataset('voltage_probe_factor',data=1000.0)#in Volts
probe_info.create_dataset('current_probe_factor',data=1000.0*2)# in Amps
rloop_area = np.pi*((probe_dia/2)**2)
tloop_area = probe_dia*hole_sep
zloop_area = probe_dia*hole_sep

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data=2e3)
run_info.create_dataset('Collection Dates',data='04112022')
run_info.create_dataset('Stuffing Delay (ms)',data=0.5)
startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=2.0


#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
testpicoshot=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(1)+').mat')

time=h5f.create_group("time")
tstart = testpicoshot['Tstart'][0][0]
time.create_dataset('tstart',data=tstart)
time.create_dataset('tinterval',data=testpicoshot['Tinterval'][0][0])
numsamples=testpicoshot['Length'][0][0]
time.create_dataset('numsamples',data=testpicoshot['Length'][0][0])
time_s = testpicoshot['Tstart'][0][0]+(np.arange(testpicoshot['Length'][0][0])*testpicoshot['Tinterval'][0][0])
time.create_dataset('time_s',data=time_s)
time.create_dataset('time_us',data=time_s*1e6)

###BEGIN PICO1 READIN#####
Discharge_raw=np.zeros([numshots,numsamples])
HV_raw=np.zeros([numshots,numsamples])
Discharge=np.zeros([numshots,numsamples])
HV=np.zeros([numshots,numsamples])
Bdot5r = np.zeros([numshots,numsamples])
Bdot5t = np.zeros([numshots,numsamples])
B5r = np.zeros([numshots,numsamples-1])
B5t = np.zeros([numshots,numsamples-1])

savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico1=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(shot)+').mat')
    Discharge_raw[savenumber,:]=pico1['B'][:,0]
    HV_raw[savenumber,:]=pico1['A'][:,0]
    Discharge[savenumber,:]=Discharge_raw[savenumber,:]*1000.0*2.0
    HV[savenumber,:]=HV_raw[savenumber,:]*1000.0
    
    Bdot5r[savenumber,:]=pico1['C'][:,0]#-np.mean(pico1['C'][0:meancutoff,0])
    Bdot5t[savenumber,:]=pico1['D'][:,0]#-np.mean(pico1['D'][0:meancutoff,0])
    
    #Remove infs
    neginfs = np.isneginf(Bdot5r[savenumber,:])
    Bdot5r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot5r[savenumber,:])
    Bdot5r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot5t[savenumber,:])
    Bdot5t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot5t[savenumber,:])
    Bdot5t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    #mean subtract 0 to meancutoff
    meansub=np.mean(Bdot5r[savenumber,0:meancutoff])
    Bdot5r[savenumber,:]=Bdot5r[savenumber,:]-meansub
    meansub=np.mean(Bdot5t[savenumber,0:meancutoff])
    Bdot5t[savenumber,:]=Bdot5t[savenumber,:]-meansub
    
    #Compute Magnetic Field for Bdots
    B5r[savenumber,:]= sp.cumtrapz(Bdot5r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    B5t[savenumber,:]= sp.cumtrapz(Bdot5t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO1 READIN #####
    
    
    
###BEGIN PICO2 READIN#####
Bdot5z = np.zeros([numshots,numsamples])
Bdot7r = np.zeros([numshots,numsamples])
Bdot7t = np.zeros([numshots,numsamples])
Bdot7z = np.zeros([numshots,numsamples])
B5z = np.zeros([numshots,numsamples-1])
B7r = np.zeros([numshots,numsamples-1])
B7t = np.zeros([numshots,numsamples-1])
B7z = np.zeros([numshots,numsamples-1])

savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico2=spio.loadmat(datadirectory+'Pico2\\20220317-0001 ('+str(shot)+').mat')

    Bdot5z[savenumber,:]=pico2['A'][:,0]#-np.mean(pico2['A'][0:meancutoff,0])
    Bdot7r[savenumber,:]=pico2['B'][:,0]#-np.mean(pico2['B'][0:meancutoff,0])
    Bdot7t[savenumber,:]=pico2['C'][:,0]#-np.mean(pico2['C'][0:meancutoff,0])
    Bdot7z[savenumber,:]=pico2['D'][:,0]#-np.mean(pico2['D'][0:meancutoff,0])
    
    #Remove infs
    neginfs = np.isneginf(Bdot5z[savenumber,:])
    Bdot5z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot5z[savenumber,:])
    Bdot5z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot7r[savenumber,:])
    Bdot7r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot7r[savenumber,:])
    Bdot7r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot7t[savenumber,:])
    Bdot7t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot7t[savenumber,:])
    Bdot7t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot7z[savenumber,:])
    Bdot7z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot7z[savenumber,:])
    Bdot7z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    #mean subtract 0 to meancutoff
    meansub=np.mean(Bdot5z[savenumber,0:meancutoff])
    Bdot5z[savenumber,:]=Bdot5z[savenumber,:]-meansub
    meansub=np.mean(Bdot7r[savenumber,0:meancutoff])
    Bdot7r[savenumber,:]=Bdot7r[savenumber,:]-meansub
    meansub=np.mean(Bdot7t[savenumber,0:meancutoff])
    Bdot7t[savenumber,:]=Bdot7t[savenumber,:]-meansub
    meansub=np.mean(Bdot7z[savenumber,0:meancutoff])
    Bdot7z[savenumber,:]=Bdot7z[savenumber,:]-meansub
    
    #Compute Magnetic Field for Bdots
    B5z[savenumber,:]= sp.cumtrapz(Bdot5z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    B7r[savenumber,:]= sp.cumtrapz(Bdot7r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    B7t[savenumber,:]= sp.cumtrapz(Bdot7t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    B7z[savenumber,:]= sp.cumtrapz(Bdot7z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####

###BEGIN PICO3 READIN#####
Bdot19r = np.zeros([numshots,numsamples])
Bdot19t = np.zeros([numshots,numsamples])
Bdot19z = np.zeros([numshots,numsamples])
Bdot21r = np.zeros([numshots,numsamples])
B19r = np.zeros([numshots,numsamples-1])
B19t = np.zeros([numshots,numsamples-1])
B19z = np.zeros([numshots,numsamples-1])
B21r = np.zeros([numshots,numsamples-1])

savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico3=spio.loadmat(datadirectory+'Pico3\\20220317-0001 ('+str(shot)+').mat')

    Bdot19r[savenumber,:]=pico3['A'][:,0]#-np.mean(pico3['A'][0:meancutoff,0])
    Bdot19t[savenumber,:]=pico3['B'][:,0]#-np.mean(pico3['B'][0:meancutoff,0])
    Bdot19z[savenumber,:]=pico3['C'][:,0]#-np.mean(pico3['C'][0:meancutoff,0])
    Bdot21r[savenumber,:]=pico3['D'][:,0]#-np.mean(pico3['D'][0:meancutoff,0])
    
    #Remove infs
    neginfs = np.isneginf(Bdot19r[savenumber,:])
    Bdot19r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot19r[savenumber,:])
    Bdot19r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot19t[savenumber,:])
    Bdot19t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot19t[savenumber,:])
    Bdot19t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot19z[savenumber,:])
    Bdot19z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot19z[savenumber,:])
    Bdot19z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot21r[savenumber,:])
    Bdot21r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot21r[savenumber,:])
    Bdot21r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    #mean subtract 0 to meancutoff
    meansub=np.mean(Bdot19r[savenumber,0:meancutoff])
    Bdot19r[savenumber,:]=Bdot19r[savenumber,:]-meansub
    meansub=np.mean(Bdot19t[savenumber,0:meancutoff])
    Bdot19t[savenumber,:]=Bdot19t[savenumber,:]-meansub
    meansub=np.mean(Bdot19z[savenumber,0:meancutoff])
    Bdot19z[savenumber,:]=Bdot19z[savenumber,:]-meansub
    meansub=np.mean(Bdot21r[savenumber,0:meancutoff])
    Bdot21r[savenumber,:]=Bdot21r[savenumber,:]-meansub
    
    #Compute Magnetic Field for Bdots
    B19r[savenumber,:]= sp.cumtrapz(Bdot19r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    B19t[savenumber,:]= sp.cumtrapz(Bdot19t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    B19z[savenumber,:]= sp.cumtrapz(Bdot19z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    B21r[savenumber,:]= sp.cumtrapz(Bdot21r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO3 READIN #####
    
###BEGIN PICO4 READIN#####
Bdot21t = np.zeros([numshots,numsamples])
Bdot21z = np.zeros([numshots,numsamples])
Bdot33r = np.zeros([numshots,numsamples])
Bdot33t = np.zeros([numshots,numsamples])
B21t = np.zeros([numshots,numsamples-1])
B21z = np.zeros([numshots,numsamples-1])
B33r = np.zeros([numshots,numsamples-1])
B33t = np.zeros([numshots,numsamples-1])

savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico4=spio.loadmat(datadirectory+'Pico4\\20220317-0001 ('+str(shot)+').mat')

    Bdot21t[savenumber,:]=pico4['A'][:,0]#-np.mean(pico4['A'][0:meancutoff,0])
    Bdot21z[savenumber,:]=pico4['B'][:,0]#-np.mean(pico4['B'][0:meancutoff,0])
    Bdot33r[savenumber,:]=pico4['C'][:,0]#-np.mean(pico4['C'][0:meancutoff,0])
    Bdot33t[savenumber,:]=pico4['D'][:,0]#-np.mean(pico4['D'][0:meancutoff,0])
    
    #Remove infs
    neginfs = np.isneginf(Bdot21t[savenumber,:])
    Bdot21t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot21t[savenumber,:])
    Bdot21t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot21z[savenumber,:])
    Bdot21z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot21z[savenumber,:])
    Bdot21z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot33r[savenumber,:])
    Bdot33r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot33r[savenumber,:])
    Bdot33r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot33t[savenumber,:])
    Bdot33t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot33t[savenumber,:])
    Bdot33t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    #mean subtract 0 to meancutoff
    meansub=np.mean(Bdot21t[savenumber,0:meancutoff])
    Bdot21t[savenumber,:]=Bdot21t[savenumber,:]-meansub
    meansub=np.mean(Bdot21z[savenumber,0:meancutoff])
    Bdot21z[savenumber,:]=Bdot21z[savenumber,:]-meansub
    meansub=np.mean(Bdot33r[savenumber,0:meancutoff])
    Bdot33r[savenumber,:]=Bdot33r[savenumber,:]-meansub
    meansub=np.mean(Bdot33t[savenumber,0:meancutoff])
    Bdot33t[savenumber,:]=Bdot33t[savenumber,:]-meansub
    
    
    #Compute Magnetic Field for Bdots
    B21t[savenumber,:]= sp.cumtrapz(Bdot21t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    B21z[savenumber,:]= sp.cumtrapz(Bdot21z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    B33r[savenumber,:]= sp.cumtrapz(Bdot33r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    B33t[savenumber,:]= sp.cumtrapz(Bdot33t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO4 READIN #####

###BEGIN PICO5 READIN#####
Bdot33z = np.zeros([numshots,numsamples])
Bdot35r = np.zeros([numshots,numsamples])
Bdot35t = np.zeros([numshots,numsamples])
Bdot35z = np.zeros([numshots,numsamples])
B33z = np.zeros([numshots,numsamples-1])
B35r = np.zeros([numshots,numsamples-1])
B35t = np.zeros([numshots,numsamples-1])
B35z = np.zeros([numshots,numsamples-1])

savenumber=0
for shot in np.arange(startshot,maxshot+1):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico5=spio.loadmat(datadirectory+'Pico5\\20220317-0001 ('+str(shot)+').mat')

    Bdot33z[savenumber,:]=pico5['A'][:,0]#-np.mean(pico5['A'][0:meancutoff,0])
    Bdot35r[savenumber,:]=pico5['B'][:,0]#-np.mean(pico5['B'][0:meancutoff,0])
    Bdot35t[savenumber,:]=pico5['C'][:,0]#-np.mean(pico5['C'][0:meancutoff,0])
    Bdot35z[savenumber,:]=pico5['D'][:,0]#-np.mean(pico5['D'][0:meancutoff,0])
    
    #Remove infs
    neginfs = np.isneginf(Bdot33z[savenumber,:])
    Bdot33z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot33z[savenumber,:])
    Bdot33z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot35r[savenumber,:])
    Bdot35r[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot35r[savenumber,:])
    Bdot35r[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot35t[savenumber,:])
    Bdot35t[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot35t[savenumber,:])
    Bdot35t[savenumber,np.where(posinfs)] = bdotmaxrange
    
    neginfs = np.isneginf(Bdot35z[savenumber,:])
    Bdot35z[savenumber,np.where(neginfs)] = -bdotmaxrange
    posinfs = np.isinf(Bdot35z[savenumber,:])
    Bdot35z[savenumber,np.where(posinfs)] = bdotmaxrange
    
    #mean subtract 0 to meancutoff
    meansub=np.mean(Bdot33z[savenumber,0:meancutoff])
    Bdot33z[savenumber,:]=Bdot33z[savenumber,:]-meansub
    meansub=np.mean(Bdot35r[savenumber,0:meancutoff])
    Bdot35r[savenumber,:]=Bdot35r[savenumber,:]-meansub
    meansub=np.mean(Bdot35t[savenumber,0:meancutoff])
    Bdot35t[savenumber,:]=Bdot35t[savenumber,:]-meansub
    meansub=np.mean(Bdot35z[savenumber,0:meancutoff])
    Bdot35z[savenumber,:]=Bdot35z[savenumber,:]-meansub
    
    #Compute Magnetic Field for Bdots
    B33z[savenumber,:]= sp.cumtrapz(Bdot33z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    B35r[savenumber,:]= sp.cumtrapz(Bdot35r[savenumber,:]/rloop_area,time_s)*1e4#Gauss
    B35t[savenumber,:]= sp.cumtrapz(Bdot35t[savenumber,:]/tloop_area,time_s)*1e4#Gauss
    B35z[savenumber,:]= sp.cumtrapz(Bdot35z[savenumber,:]/zloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO5 READIN ####

#create magnetic probe group
mag = h5f.create_group("mag_probe")
positions = mag.create_group("positions")
loc5 = positions.create_group("probe5")
r5=loc5.create_group("r")
t5=loc5.create_group("t")
z5=loc5.create_group("z")
r5.create_dataset("bdot",data=Bdot5r)
r5.create_dataset("b",data=B5r)
t5.create_dataset("bdot",data=Bdot5t)
t5.create_dataset("b",data=B5t)
z5.create_dataset("bdot",data=Bdot5z)
z5.create_dataset("b",data=B5z)

loc7 = positions.create_group("probe7")
r7=loc7.create_group("r")
t7=loc7.create_group("t")
z7=loc7.create_group("z")
r7.create_dataset("bdot",data=Bdot7r)
r7.create_dataset("b",data=B7r)
t7.create_dataset("bdot",data=Bdot7t)
t7.create_dataset("b",data=B7t)
z7.create_dataset("bdot",data=Bdot7z)
z7.create_dataset("b",data=B7z)

loc19 = positions.create_group("probe19")
r19=loc19.create_group("r")
t19=loc19.create_group("t")
z19=loc19.create_group("z")
r19.create_dataset("bdot",data=Bdot19r)
r19.create_dataset("b",data=B19r)
t19.create_dataset("bdot",data=Bdot19t)
t19.create_dataset("b",data=B19t)
z19.create_dataset("bdot",data=Bdot19z)
z19.create_dataset("b",data=B19z)

loc21 = positions.create_group("probe21")
r21=loc21.create_group("r")
t21=loc21.create_group("t")
z21=loc21.create_group("z")
r21.create_dataset("bdot",data=Bdot21r)
r21.create_dataset("b",data=B21r)
t21.create_dataset("bdot",data=Bdot21t)
t21.create_dataset("b",data=B21t)
z21.create_dataset("bdot",data=Bdot21z)
z21.create_dataset("b",data=B21z)

loc33 = positions.create_group("probe33")
r33=loc33.create_group("r")
t33=loc33.create_group("t")
z33=loc33.create_group("z")
r33.create_dataset("bdot",data=Bdot33r)
r33.create_dataset("b",data=B33r)
t33.create_dataset("bdot",data=Bdot33t)
t33.create_dataset("b",data=B33t)
z33.create_dataset("bdot",data=Bdot33z)
z33.create_dataset("b",data=B33z)

loc35 = positions.create_group("probe35")
r35=loc35.create_group("r")
t35=loc35.create_group("t")
z35=loc35.create_group("z")
r35.create_dataset("bdot",data=Bdot35r)
r35.create_dataset("b",data=B35r)
t35.create_dataset("bdot",data=Bdot35t)
t35.create_dataset("b",data=B35t)
z35.create_dataset("bdot",data=Bdot35z)
z35.create_dataset("b",data=B35z)
                       
dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()
