#data_conversion_hd5f_06192025_twodens.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""
import numpy as np
import scipy.integrate as sp
import pandas as pd

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2025\\07242025\\'

#dataset1
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive0cm' #shots 15 to 24
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive0p5cm' #shots 25 to 34
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive1cm' #shots 35 to 44
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive1p5cm' #shots 45 to 54
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive2cm' #shots 55 to 64
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive2p5cm' #shots 65 to 74
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive3cm' #shots 75 to 85 skip 79
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive3p5cm' #shots 86 to 104 skip 90,92,93,94,95,96,97,100,103
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive4cm' #shots 105 to 115 skip 111
runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive4p5cm' #shots 116 to 126 skip 123

startshot=116
maxshot=126+1
skipshots=[123]
"""
5 probes:
    3 Triplet Bdot
        Index0 Bdot in first KF16 port between quick connect port arrays (Pulled out 1.5in from center)
        Index1 Bdot in Port 15 (in front of blocker)
        Index2 Bdot in Port 19 (behind blocker)
    2 Langmuir Probes
        Index0 isat in Port 2
        Index1 isat on Probe Drive in 3rd KF40 port on Top (about Port 25)
H-alpha detector at 3rd KF16 port on top
MaCor plate at port 17

Langmuir Probe circuit: 6 12V batteries in series (~78V) across 100ohm resistor, diff amp at div100
    
Base Pressure ~1x10^-7 Torr

20us/div 14bit Trig 2V Ext 3% Rising
Picoscope 1 Chan A HV Probe div1000 at +/- 2V
            Chan B Current Probe div1000 div2 at +/- 5V
            Chan C Isat Port 2 Double Probe, Pos Tip shadowed, Res=100ohm, +/-1V (no motion)
            Chan D Isat Probe Drive Double Probe, Pos Tip shadowed, Res=100ohm, +/-1V 

Bdots r, theta, z in Chan A, B, C (all 50ohm term)
Picoscope 2 Bdot in KF16 at +/- 1V behind plate pulled out 1.5in from center 
Picoscope 3 Bdot at port 15 at +/- 1V in front of plate
Picoscope 4 Bdot at port 19 at +/- 1V behind plate

1ms/div 16bit
Picoscope 5 H-alpha filtered light detector Chan A +/-10V with 50ohm term
"""

from scipy.signal import butter, sosfiltfilt

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
probe_info.create_dataset('radial_location',data='center')
probe_info.create_dataset('voltage_probe_factor',data=1000.0)#in Volts
probe_info.create_dataset('current_probe_factor',data=1000.0*2)# in Amps

rloop_area = np.pi*((probe_dia/2)**2)
tloop_area = probe_dia*hole_sep
zloop_area = probe_dia*hole_sep

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data=2.5e3)
run_info.create_dataset('Collection Dates',data='03272025')
run_info.create_dataset('Stuffing Delay (ms)',data=-0.5)
run_info.create_dataset('Gas Open (ms)',data=-2.0)
run_info.create_dataset('Gas Close (ms)',data=0.001)

startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=2#mV


#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
#testpicoshot=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(1)+').mat')
#testpicoshot=np.loadtxt(datadirectory+'\\pico1\\20230627-0001 ('+str(10)+').txt',skiprows = 2, unpack = True)
testpicoshot=pd.read_csv(datadirectory+'\\Pico1\\20250724-0001 ('+str(10)+').csv',header=[0],skiprows=[1])


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

testpicoshot=pd.read_csv(datadirectory+'\\Pico5\\20250724-0001 ('+str(50)+').csv',header=[0],skiprows=[1])
tstart2 = testpicoshot['Time'][0]*1e-3#pico reads out in ms, convert to seconds
#time.create_dataset('tstart',data=tstart)
tinterval2=np.abs((testpicoshot['Time'][1]-testpicoshot['Time'][0])*1e-3)
#time.create_dataset('tinterval',data=tinterval)
numsamples2=int(np.shape(testpicoshot['Time'])[0])
#time.create_dataset('numsamples',data=numsamples)
time_s2 = tstart2+(np.arange(numsamples2)*tinterval2)
time.create_dataset('timeHA_s',data=time_s2)
time.create_dataset('timeHA_ms',data=time_s2*1e3)

###BEGIN Scan READIN#####
numshots=10
numprobes_trip = 3
numprobes_isat = 2
Discharge_raw=np.zeros([numshots,numsamples])
HV_raw=np.zeros([numshots,numsamples])
Discharge=np.zeros([numshots,numsamples])
HV=np.zeros([numshots,numsamples])
#wirecurrent=np.zeros([num_rad_pos,num_currents,numshots,numsamples_wc])
Bdotr = np.zeros([numprobes_trip,numshots,numsamples])
Bdott = np.zeros([numprobes_trip,numshots,numsamples])
Bdotz = np.zeros([numprobes_trip,numshots,numsamples])
Br = np.zeros([numprobes_trip,numshots,numsamples-1])
Bt = np.zeros([numprobes_trip,numshots,numsamples-1])
Bz = np.zeros([numprobes_trip,numshots,numsamples-1])
isat = np.zeros([numprobes_isat,numshots,numsamples])
halpha = np.zeros([numshots,numsamples])
term50=2.0

#

### Extract PICO 1 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\\Pico1\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 1: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-5,inplace=True)
    pico.replace("∞",5,inplace=True)

    HV_raw[savenumber,:]=np.array(pico['Channel A'],dtype=float)
    Discharge_raw[savenumber,:]=np.array(pico['Channel B'],dtype=float)
    Discharge[savenumber,:]=Discharge_raw[savenumber,:]*1000.0*2.0
    HV[savenumber,:]=HV_raw[savenumber,:]*1000.0
    
    pico=pd.read_csv(datadirectory+'\\Pico1\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-1,inplace=True)
    pico.replace("∞",1,inplace=True)
    
    isat[0,savenumber,:]=np.array(pico['Channel C'],dtype=float)*100.0/100.0#Amps
    isat[1,savenumber,:]=np.array(pico['Channel D'],dtype=float)*100.0/100.0#Amps
    
    #End of loop
    savenumber+=1
##END PICO1 READIN ##### 
   
### Extract PICO 2 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\\Pico2\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 2: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-1,inplace=True)
    pico.replace("∞",1,inplace=True)
    
    Bdotr[0,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50#*1e-3#to convert to volts
    Bdott[0,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50#*1e-3#to convert to volts
    Bdotz[0,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50#*1e-3#to convert to volts
    
    #filter Bdot
    Bdotr[0,savenumber,:]=butter_bandpass_filter(Bdotr[0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[0,savenumber,:]=butter_bandpass_filter(Bdott[0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[0,savenumber,:]=butter_bandpass_filter(Bdotz[0,savenumber,:],lowcut,highcut,fs,order=9)

    #Compute Magnetic Field for Bdots
    Br[0,savenumber,:]= sp.cumulative_trapezoid(Bdotr[0,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[0,savenumber,:]= sp.cumulative_trapezoid(Bdott[0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[0,savenumber,:]= sp.cumulative_trapezoid(Bdotz[0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####
    
### Extract PICO 3 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\\Pico3\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 3: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    #pico.replace("-∞",-200,inplace=True)
    #pico.replace("∞",200,inplace=True)
    #dataset3 and on:
    pico.replace("-∞",-1,inplace=True)
    pico.replace("∞",1,inplace=True)
    
    Bdotr[1,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50#*1e-3#to convert to volts
    Bdott[1,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50#*1e-3#to convert to volts
    Bdotz[1,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50#*1e-3#to convert to volts
    
    #filter Bdot
    Bdotr[1,savenumber,:]=butter_bandpass_filter(Bdotr[1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[1,savenumber,:]=butter_bandpass_filter(Bdott[1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[1,savenumber,:]=butter_bandpass_filter(Bdotz[1,savenumber,:],lowcut,highcut,fs,order=9)

    #Compute Magnetic Field for Bdots
    Br[1,savenumber,:]= sp.cumulative_trapezoid(Bdotr[1,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[1,savenumber,:]= sp.cumulative_trapezoid(Bdott[1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[1,savenumber,:]= sp.cumulative_trapezoid(Bdotz[1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO3 READIN #####

### Extract PICO 4 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\\Pico4\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 4: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    #pico.replace("-∞",-200,inplace=True)
    #pico.replace("∞",200,inplace=True)
    #dataset3 and on:
    pico.replace("-∞",-1,inplace=True)
    pico.replace("∞",1,inplace=True)
    
    Bdotr[2,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50#*1e-3#to convert to volts
    Bdott[2,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50#*1e-3#to convert to volts
    Bdotz[2,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50#*1e-3#to convert to volts
    
    #filter Bdot
    Bdotr[2,savenumber,:]=butter_bandpass_filter(Bdotr[2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[2,savenumber,:]=butter_bandpass_filter(Bdott[2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[2,savenumber,:]=butter_bandpass_filter(Bdotz[2,savenumber,:],lowcut,highcut,fs,order=9)

    #Compute Magnetic Field for Bdots
    Br[2,savenumber,:]= sp.cumulative_trapezoid(Bdotr[2,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[2,savenumber,:]= sp.cumulative_trapezoid(Bdott[2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[2,savenumber,:]= sp.cumulative_trapezoid(Bdotz[2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO4 READIN #####


### Extract PICO 5 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\\Pico5\\20250724-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 5: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    #pico.replace("-∞",-200,inplace=True)
    #pico.replace("∞",200,inplace=True)
    #dataset3 and on:
    pico.replace("-∞",-10,inplace=True)
    pico.replace("∞",10,inplace=True)
    
    #halpha[savenumber,:]=np.array(pico['Channel A'],dtype=float)
    
    #End of loop
    savenumber+=1
##END PICO4 READIN #####


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

#create langmuir probe group
isat_probe=h5f.create_group("isat probe")
isat_probe.create_dataset("isat",data=isat)

#create H-alpha probe group
halpha_probe=h5f.create_group('HAlpha')
halpha_probe.create_dataset("HA",data=halpha)
                       
#create discharage group
dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()
