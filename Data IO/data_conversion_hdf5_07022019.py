#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""

import load_picoscope_bmx_07022019 as load
import numpy as np

filepath='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\07022019\\processed\\'
filename='2kV_oddpos9to21zr_2ms_stuffdelay_40_07022019'
import h5py
h5f = h5py.File(filepath+filename+'.h5','a')

probe_dia = 0.00158755#m (1/16'' probe)
hole_sep = 0.001016#m (1/16''probe)
probe_info=h5f.create_group("probe_info")
probe_info.create_dataset('probe_stalk_diameter',data=probe_dia)
probe_info.create_dataset('probe_hole_separation',data=hole_sep)
probe_info.create_dataset('rloop_area',data=np.pi*(probe_dia/2)**2)
probe_info.create_dataset('tloop_area',data=probe_dia*hole_sep)
probe_info.create_dataset('zloop_area',data=probe_dia*hole_sep)
probe_info.create_dataset('radial_location',data='center')

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data=2e3)
run_info.create_dataset('Collection Dates',data='07022019')
run_info.create_dataset('Stuffing Delay (ms)',data=2.0)

#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
data = load.load_picoscope(1)
shotnum = 40
time_us,time_s,timeB_s,timeB_us,timeraw=data[0],data[1],data[2],data[3],data[20]
time=h5f.create_group("time")
time.create_dataset('time_us',data=time_us)
time.create_dataset('time_s',data=time_s)
time.create_dataset('timeB_us',data=timeB_us)
time.create_dataset('timeraw',data=timeraw)

datalength=data[0].shape[0]
datalengthraw=data[20].shape[0]

Bdot1raw=np.zeros([shotnum,datalengthraw])
Bdot2raw=np.zeros([shotnum,datalengthraw])
Bdot3raw=np.zeros([shotnum,datalengthraw])
Bdot4raw=np.zeros([shotnum,datalengthraw])
Bdot1=np.zeros([shotnum,datalength])
Bdot2=np.zeros([shotnum,datalength])
Bdot3=np.zeros([shotnum,datalength])
Bdot4=np.zeros([shotnum,datalength])
B1=np.zeros([shotnum,datalength-1])
B2=np.zeros([shotnum,datalength-1])
B3=np.zeros([shotnum,datalength-1])
B4=np.zeros([shotnum,datalength-1])

#Picoscope1
for shot in np.arange(1,shotnum+1):
    data=load.load_picoscope(shot,scopenum=1)
    Bdot1raw[shot-1,:]=data[16]
    Bdot2raw[shot-1,:]=data[17]
    Bdot3raw[shot-1,:]=data[18]
    Bdot4raw[shot-1,:]=data[19]
    Bdot1[shot-1,:]=data[4]
    Bdot2[shot-1,:]=data[5]
    Bdot3[shot-1,:]=data[6]
    Bdot4[shot-1,:]=data[7]
    B1[shot-1,:]=data[8]
    B2[shot-1,:]=data[9]
    B3[shot-1,:]=data[10]
    B4[shot-1,:]=data[11]

pos=h5f.create_group("pos9")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot1raw)
r=raw.create_dataset('r',data=Bdot2raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot1)
r=bdot.create_dataset('r',data=Bdot2)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B1)
r=b.create_dataset('r',data=B2)

pos=h5f.create_group("pos11")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot3raw)
r=raw.create_dataset('r',data=Bdot4raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot3)
r=bdot.create_dataset('r',data=Bdot4)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B3)
r=b.create_dataset('r',data=B4)

#Picoscope2
for shot in np.arange(1,shotnum+1):
    data=load.load_picoscope(shot,scopenum=2)
    Bdot1raw[shot-1,:]=data[16]
    Bdot2raw[shot-1,:]=data[17]
    Bdot3raw[shot-1,:]=data[18]
    Bdot4raw[shot-1,:]=data[19]
    Bdot1[shot-1,:]=data[4]
    Bdot2[shot-1,:]=data[5]
    Bdot3[shot-1,:]=data[6]
    Bdot4[shot-1,:]=data[7]
    B1[shot-1,:]=data[8]
    B2[shot-1,:]=data[9]
    B3[shot-1,:]=data[10]
    B4[shot-1,:]=data[11]

pos=h5f.create_group("pos13")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot1raw)
r=raw.create_dataset('r',data=Bdot2raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot1)
r=bdot.create_dataset('r',data=Bdot2)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B1)
r=b.create_dataset('r',data=B2)

pos=h5f.create_group("pos15")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot3raw)
r=raw.create_dataset('r',data=Bdot4raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot3)
r=bdot.create_dataset('r',data=Bdot4)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B3)
r=b.create_dataset('r',data=B4)

#Picoscope3
for shot in np.arange(1,shotnum+1):
    data=load.load_picoscope(shot,scopenum=3)
    Bdot1raw[shot-1,:]=data[16]
    Bdot2raw[shot-1,:]=data[17]
    Bdot3raw[shot-1,:]=data[18]
    Bdot4raw[shot-1,:]=data[19]
    Bdot1[shot-1,:]=data[4]
    Bdot2[shot-1,:]=data[5]
    Bdot3[shot-1,:]=data[6]
    Bdot4[shot-1,:]=data[7]
    B1[shot-1,:]=data[8]
    B2[shot-1,:]=data[9]
    B3[shot-1,:]=data[10]
    B4[shot-1,:]=data[11]

pos=h5f.create_group("pos17")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot1raw)
r=raw.create_dataset('r',data=Bdot2raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot1)
r=bdot.create_dataset('r',data=Bdot2)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B1)
r=b.create_dataset('r',data=B2)

pos=h5f.create_group("pos19")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot3raw)
r=raw.create_dataset('r',data=Bdot4raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot3)
r=bdot.create_dataset('r',data=Bdot4)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B3)
r=b.create_dataset('r',data=B4)

#Picoscope4
for shot in np.arange(1,shotnum+1):
    data=load.load_picoscope(shot,scopenum=4)
    Bdot1raw[shot-1,:]=data[16]
    Bdot2raw[shot-1,:]=data[17]
    Bdot3raw[shot-1,:]=data[18]
    Bdot4raw[shot-1,:]=data[19]
    Bdot1[shot-1,:]=data[4]
    Bdot2[shot-1,:]=data[5]
    Bdot3[shot-1,:]=data[6]
    Bdot4[shot-1,:]=data[7]
    B1[shot-1,:]=data[8]
    B2[shot-1,:]=data[9]
    B3[shot-1,:]=data[10]
    B4[shot-1,:]=data[11]

pos=h5f.create_group("pos21")
raw=pos.create_group("raw")
theta=raw.create_dataset('theta',data=np.zeros([shotnum,datalengthraw]))
z=raw.create_dataset('z',data=Bdot1raw)
r=raw.create_dataset('r',data=Bdot2raw)
bdot=pos.create_group("bdot")
theta=bdot.create_dataset('theta',data=np.zeros([shotnum,datalength]))
z=bdot.create_dataset('z',data=Bdot1)
r=bdot.create_dataset('r',data=Bdot2)
b=pos.create_group("b")
theta=b.create_dataset('theta',data=np.zeros([shotnum,datalength-1]))
z=b.create_dataset('z',data=B1)
r=b.create_dataset('r',data=B2)

import matplotlib.pylab as plt
channels=['pos9','pos11','pos13','pos15','pos17','pos19','pos21']
#for shot in np.arange(shotnum):
for channel in channels:
    plt.plot(h5f[channel]['b']['z'][0,:])

