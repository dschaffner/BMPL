#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""

import load_picoscope_bmx_10142019_DisHV as load
import numpy as np

filepath='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\10142019\\processed\\'
filename='2kV_LangPort4_0p5ms_1msGas_8shots_bdot_10142019'
filename='2kV_0p5ms_1msGas_20shots_bdot_10142019' #shots 7 to 26
filename='2kV_0p5ms_10msGas_20shots_bdot_10142019' #27 to 47 skip 30
filename='2kV_0p5ms_1msGas_20shots_DisHV_10142019' #shots 7 to 26
filename='2kV_0p5ms_10msGas_20shots_DisHV_10142019' #27 to 47 skip 30

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
run_info.create_dataset('Collection Dates',data='07052019')
run_info.create_dataset('Stuffing Delay (ms)',data=0.5)

#time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,Bdot4,B1,B2,B3,B4,B1filt,B2filt,B3filt,B4filt
data = load.load_picoscope(1,scopenum=4)


shotnum = 21
time_us,time_s=data[0],data[1]
time=h5f.create_group("time")
time.create_dataset('time_us',data=time_us)
time.create_dataset('time_s',data=time_s)

datalength=data[0].shape[0]

Discharge_raw=np.zeros([shotnum,datalength])
HV_raw=np.zeros([shotnum,datalength])
Discharge=np.zeros([shotnum,datalength])
HV=np.zeros([shotnum,datalength])

#Picoscope4 only
for shot in np.arange(27,48):
    data=load.load_picoscope(shot,scopenum=4)
    Discharge_raw[shot-27,:]=data[2]
    HV_raw[shot-27,:]=data[3]
    Discharge[shot-27,:]=data[4]
    HV[shot-27,:]=data[5]


dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()