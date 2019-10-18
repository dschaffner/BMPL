#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""

import load_picoscope_bmx_10142019_forLang as load
import numpy as np

filepath='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\10142019\\processed\\'
filename='2kV_LangPort4_0p5ms_1msGas_8shots_langmuircurrent_10142019'
import h5py
h5f = h5py.File(filepath+filename+'.h5','a')

probe_dia = 0.00158755#m (1/16'' probe)
hole_sep = 0.001016#m (1/16''probe)
probe_info=h5f.create_group("probe_info")
probe_info.create_dataset('Current Monitor Factor',data='10A/V')
#probe_info.create_dataset('probe_stalk_diameter',data=probe_dia)
#probe_info.create_dataset('probe_hole_separation',data=hole_sep)
#probe_info.create_dataset('rloop_area',data=np.pi*(probe_dia/2)**2)
#probe_info.create_dataset('tloop_area',data=probe_dia*hole_sep)
#probe_info.create_dataset('zloop_area',data=probe_dia*hole_sep)
#probe_info.create_dataset('radial_location',data='center')

run_info=h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)',data=2e3)
run_info.create_dataset('Collection Dates',data='07162019')
run_info.create_dataset('Stuffing Delay (ms)',data=0.5)

#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
data = load.load_picoscope(6,scopenum=4)
shotnum = 8
time_us,time_s,timeB_s,timeB_us=data[0],data[1],data[2],data[3]
time=h5f.create_group("time")
time.create_dataset('time_us',data=time_us)
time.create_dataset('time_s',data=time_s)
#time.create_dataset('timeB_us',data=timeB_us)
#time.create_dataset('timeraw',data=timeraw)

datalength=data[0].shape[0]
#datalengthraw=data[20].shape[0]

#Bdot1raw=np.zeros([shotnum,datalengthraw])
#Bdot2raw=np.zeros([shotnum,datalengthraw])
#Bdot3raw=np.zeros([shotnum,datalengthraw])
#Bdot4raw=np.zeros([shotnum,datalengthraw])
lang_isat=np.zeros([shotnum,datalength])
#Bdot2=np.zeros([shotnum,datalength])
#Bdot3=np.zeros([shotnum,datalength])
#Bdot4=np.zeros([shotnum,datalength])
#B1=np.zeros([shotnum,datalength-1])
#B2=np.zeros([shotnum,datalength-1])
#B3=np.zeros([shotnum,datalength-1])
#B4=np.zeros([shotnum,datalength-1])

#Picoscope1
for shot in np.arange(6,14):
    data=load.load_picoscope(shot,scopenum=4)

    lang_isat[shot-6,:]=data[4]*10.0*5.0


lang_isat=h5f.create_dataset('lang_isat',data=lang_isat)

import matplotlib.pylab as plt
plt.plot(time_us,h5f['lang_isat'][0,:])
plt.plot(time_us,h5f['lang_isat'][1,:])
plt.plot(time_us,h5f['lang_isat'][2,:])
plt.plot(time_us,h5f['lang_isat'][3,:])
plt.plot(time_us,h5f['lang_isat'][4,:])
plt.plot(time_us,h5f['lang_isat'][5,:])
plt.plot(time_us,h5f['lang_isat'][6,:])
plt.plot(time_us,h5f['lang_isat'][7,:])