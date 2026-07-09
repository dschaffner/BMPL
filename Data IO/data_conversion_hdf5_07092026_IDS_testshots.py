#data_conversion_hd5f_06192025_twodens.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""
import numpy as np
import scipy.integrate as sp
import pandas as pd
import sys
import matplotlib.pylab as plot

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2026\\07082026\\'

#dataset1
runfilename = 'IDS_testshots'
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive0cm' #shots 15 to 24
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive0p5cm' #shots 25 to 34
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive1cm' #shots 35 to 44
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive1p5cm' #shots 45 to 54
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive2cm' #shots 55 to 64
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive2p5cm' #shots 65 to 74
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive3cm' #shots 75 to 85 skip 79
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive3p5cm' #shots 86 to 104 skip 90,92,93,94,95,96,97,100,103
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive4cm' #shots 105 to 115 skip 111
#runfilename = 'TwoDens_ThreeBdot_2kV_p8msstuff_2mst1usgas_halfplateblocker_10shots_probedrive4p5cm' #shots 116 to 126 skip 123

startshot=1
maxshot=8+1
skipshots=[]
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
#probe_dia = 0.00158755#m (1/16'' probe)
#hole_sep = 0.001016#m (1/16''probe)

#probe_info=h5f.create_group("probe_info")
#probe_info.create_dataset('probe_stalk_diameter',data=probe_dia)
#probe_info.create_dataset('probe_hole_separation',data=hole_sep)
#probe_info.create_dataset('rloop_area',data=np.pi*(probe_dia/2)**2)
#probe_info.create_dataset('tloop_area',data=probe_dia*hole_sep)
#probe_info.create_dataset('zloop_area',data=probe_dia*hole_sep)
#probe_info.create_dataset('radial_location',data='center')
#probe_info.create_dataset('voltage_probe_factor',data=1000.0)#in Volts
#probe_info.create_dataset('current_probe_factor',data=1000.0*2)# in Amps

#rloop_area = np.pi*((probe_dia/2)**2)
#tloop_area = probe_dia*hole_sep
#zloop_area = probe_dia*hole_sep

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data='variable')
run_info.create_dataset('Collection Dates',data='07082026')
run_info.create_dataset('Stuffing Delay (ms)',data=-1.0)
run_info.create_dataset('Gas Open (ms)',data=-2.0)
run_info.create_dataset('Gas Close (ms)',data=0.001)

startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=2#mV




#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
#testpicoshot=spio.loadmat(datadirectory+'Pico1\\20220317-0001 ('+str(1)+').mat')
#testpicoshot=np.loadtxt(datadirectory+'\\pico1\\20230627-0001 ('+str(10)+').txt',skiprows = 2, unpack' = True)
testpicoshot=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(2)+').csv',
                         header='infer',
                         names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                         skiprows=2)



#pico1 - 62503 rows
#pico8_1 - 20004 rows


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


###BEGIN Scan READIN#####
numshots=7
channels=16
voltage_range = 1.0
ids = np.zeros([numshots,channels,numsamples])

#load shot 1
shotcount = 0
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)

#load shot 2
shotcount = 1
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)


#load shot 3
shotcount = 2
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)


#load shot 4
shotcount = 3
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)


#load shot 5
shotcount = 4
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)


#load shot 6
shotcount = 5
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)


#load shot 7
shotcount = 6
shotfilenum = shotcount+2
pico1=pd.read_csv(datadirectory+'\\Pico8_1\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)
pico2=pd.read_csv(datadirectory+'\\Pico8_2\\20260708-0001 ('+str(shotfilenum)+').csv',
                 header='infer',
                 names=['Time','Channel A', 'Channel B', 'Channel C', 'Channel D', 'Channel E', 'Channel F', 'Channel G', 'Channel H'],
                 skiprows=2)

#replace infinities (from clipping) with max range values
pico1.replace("-∞",-voltage_range,inplace=True)
pico1.replace("∞",voltage_range,inplace=True)
pico2.replace("-∞",-voltage_range,inplace=True)
pico2.replace("∞",voltage_range,inplace=True)

print('Pico 8: loading shot number ', 2)
ids[shotcount,0,:]=np.array(pico1['Channel A'],dtype=float)
ids[shotcount,1,:]=np.array(pico1['Channel B'],dtype=float)
ids[shotcount,2,:]=np.array(pico1['Channel C'],dtype=float)
ids[shotcount,3,:]=np.array(pico1['Channel D'],dtype=float)
ids[shotcount,4,:]=np.array(pico1['Channel E'],dtype=float)
ids[shotcount,5,:]=np.array(pico1['Channel F'],dtype=float)
ids[shotcount,6,:]=np.array(pico1['Channel G'],dtype=float)
ids[shotcount,7,:]=np.array(pico1['Channel H'],dtype=float)
ids[shotcount,8,:]=np.array(pico2['Channel A'],dtype=float)
ids[shotcount,9,:]=np.array(pico2['Channel B'],dtype=float)
ids[shotcount,10,:]=np.array(pico2['Channel C'],dtype=float)
ids[shotcount,11,:]=np.array(pico2['Channel D'],dtype=float)
ids[shotcount,12,:]=np.array(pico2['Channel E'],dtype=float)
ids[shotcount,13,:]=np.array(pico2['Channel F'],dtype=float)
ids[shotcount,14,:]=np.array(pico2['Channel G'],dtype=float)
ids[shotcount,15,:]=np.array(pico2['Channel H'],dtype=float)

##END PICO READIN #####


#create ids group
idsgroup = h5f.create_group("ids")
ids_array=idsgroup.create_dataset("ids_array",data=ids)

h5f.close()
