# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 20:05:20 2022

@author: Josh0
"""

import load_picoscope_one_01122022 as load1
import load_picoscope_rest_01122022 as loadrest
import numpy as np
import h5py
import matplotlib.pylab as plt

#####################################################################################
"""
    The purpose of this script is to generate an HDF5 file from the 01122022 dataset.
"""
#####################################################################################

###################################################
""" The picoscope-probe structure 01122022 """
###################################################
"""          A           B           C           D
pico1       HV          DC          5R          5T
pico2       5Z          7R          7T          7Z
pico3       19R         19T         19Z         21R
pico4       21T         21Z         33R         33T
pico5       33Z         35R         35T         35Z

DC --- Discharge Current
HV --- High Voltage
Total Shots = 111 (some shots are bad, notation to follow)
"""
###################################################


new_filepath = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
new_filename = '1122022_111shots_somebad_plotting_attempt1'
h5f = h5py.File(new_filepath + new_filename + '.h5', 'a')


probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info = h5f.create_group('Probe Info')
probe_info.create_dataset('Probe Stalk Diameter', data=probe_dia)
probe_info.create_dataset('Probe_Hole_Seperation', data=hole_sep)
probe_info.create_dataset('rloop_area', data = np.pi*(probe_dia/2)**2)
probe_info.create_dataset('tloop_area', data = probe_dia*hole_sep)
probe_info.create_dataset('zloop_area', data = probe_dia*hole_sep)
probe_info.create_dataset('radial_location', data = 'center')

run_info = h5f.create_group('Run Info')
run_info.create_dataset('Discharge Voltage (kV)', data=2.5)
run_info.create_dataset('Collection Dates', data='01122022')
run_info.create_dataset('Stuffing Delay (ms)', data = "check this")


# Initial output of the load_picoscope_one():
# time_s, time_us, timeB_s, timeB_us, DisCurrent_raw, HighV_raw, DisCurrent, HighV, Bdot1raw, Bdot2raw, Bdot1, Bdot2, B1, B2, B1filt, B2filt
data_pico1=load1.load_picoscope_one_01122022(1, scopenum=1)
data_pico_rest = loadrest.load_picoscope_rest_01122022(1, scopenum=2)
shotnum = 112
time_s, time_us, timeB_s, timeB_us = data_pico1[0],data_pico1[1],data_pico1[2],data_pico1[3]

# This is the output of the load_picoscope_rest():
# time_s, time_us, timeB_s, timeB_us, Bdot3raw, Bdot4raw, Bdot5raw, Bdot6raw, Bdot3, Bdot4, Bdot5, Bdot6, B3, B4, B5, B6, B3filt, B4filt, B5filt, B6filt

# Creation of time HDF5 group
time = h5f.create_group("time")
time.create_dataset('time_s', data = time_s)
time.create_dataset('time_us', data = time_us)
time.create_dataset('timeB_s', data = timeB_s)
time.create_dataset('timeB_us', data = timeB_us)

#datalength1 = data_pico1[0]
#datalengthraw1 = data_pico1[15]

#datalength_rest = data_pico_rest[0]
#datalengthraw_rest = data_pico_rest[19]

datalength1 = int(25000)
datalengthraw1 = int(25000)
datalength_rest = int(25000)
datalengthraw_rest = int(25000)

# The following initialization makes it easier to tweak the code for different datasets
DisCurrent_raw = np.zeros([shotnum,datalengthraw_rest])
HighV_raw = np.zeros([shotnum,datalengthraw_rest])

DisCurrent = np.zeros([shotnum,datalength1])
HighV = np.zeros([shotnum,datalength1])

Bdot1raw=np.zeros([shotnum,datalengthraw1])
Bdot2raw=np.zeros([shotnum,datalengthraw1])
Bdot3raw=np.zeros([shotnum,datalengthraw_rest])
Bdot4raw=np.zeros([shotnum,datalengthraw_rest])
Bdot5raw=np.zeros([shotnum,datalengthraw_rest])
Bdot6raw=np.zeros([shotnum,datalengthraw_rest])


Bdot1=np.zeros([shotnum,datalength1])
Bdot2=np.zeros([shotnum,datalength1])
Bdot3=np.zeros([shotnum,datalength_rest])
Bdot4=np.zeros([shotnum,datalength_rest])
Bdot5=np.zeros([shotnum,datalength_rest])
Bdot6=np.zeros([shotnum,datalength_rest])

B1=np.zeros([shotnum,datalength1-1])
B2=np.zeros([shotnum,datalength1-1])
B3=np.zeros([shotnum,datalength_rest-1])
B4=np.zeros([shotnum,datalength_rest-1])
B5=np.zeros([shotnum,datalength_rest-1])
B6=np.zeros([shotnum,datalength_rest-1])

B1filt=np.zeros([shotnum,datalength1-1])
B2filt=np.zeros([shotnum,datalength1-1])
B3filt=np.zeros([shotnum,datalength_rest-1])
B4filt=np.zeros([shotnum,datalength_rest-1])
B5filt=np.zeros([shotnum,datalength_rest-1])
B6filt=np.zeros([shotnum,datalength_rest-1])


# Picoscope 1, will only use data_pico1
for shot in np.arange (1, 112):
    data_pico1 = load1.load_picoscope_one_01122022(shot, scopenum=1)
    
    DisCurrent_raw[shot,:] = data_pico1[4].T
    HighV_raw[shot,:] = data_pico1[5].T
    DisCurrent[shot,:] = data_pico1[6]
    HighV[shot,:] = data_pico1[7]
    Bdot1raw[shot,:] = data_pico1[8].T
    Bdot2raw[shot,:] = data_pico1[9].T
    Bdot1[shot,:] = data_pico1[10].T
    Bdot2[shot,:] = data_pico1[11].T
    B1[shot,:] = data_pico1[12]
    B2[shot,:] = data_pico1[13]
    B1filt[shot,:] = data_pico1[14]
    B2filt[shot,:] = data_pico1[15]
    
# Take care of the discharge data first, then B fluctuations
dis = h5f.create_group("Discharge Data")
disraw = dis.create_group("Discharge Raw")
DisCurrent_raw = disraw.create_dataset("DisCurrent_raw", data=DisCurrent_raw)
HighV_raw = disraw.create_dataset("HighV_raw", data=HighV_raw)
Dis = dis.create_group("Discharge")
DisCurrent = Dis.create_dataset("DisCurrent", data=DisCurrent)
HighV = Dis.create_dataset("HighV", data=HighV)

# Position 1 B fluctuations now; for Pico1, only R and T components of Pos1
pos = h5f.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw)
theta = raw.create_dataset('theta', data=Bdot2raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot1)
theta = bdot.create_dataset('theta', data=Bdot2)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1)
theta = b.create_dataset('theta', data=B2)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt)
theta =  Bfilt.create_dataset('theta', data=B2filt)


# Picoscope 2
for shot in np.arange(1, 112):
    data_pico_rest = loadrest.load_picoscope_rest_01122022(shot, scopenum=2)
    
    Bdot3raw[shot,:] = data_pico_rest[4].T
    Bdot4raw[shot,:] = data_pico_rest[5].T
    Bdot5raw[shot,:] = data_pico_rest[6].T
    Bdot6raw[shot,:] = data_pico_rest[7].T
    Bdot3[shot,:] = data_pico_rest[8].T
    Bdot4[shot,:] = data_pico_rest[9].T
    Bdot5[shot,:] = data_pico_rest[10].T
    Bdot6[shot,:] = data_pico_rest[11].T
    B3[shot,:] = data_pico_rest[12]
    B4[shot,:] = data_pico_rest[13]
    B5[shot,:] = data_pico_rest[14]
    B6[shot,:] = data_pico_rest[15]
    B3filt[shot,:] = data_pico_rest[16]
    B4filt[shot,:] = data_pico_rest[17]
    B5filt[shot,:] = data_pico_rest[18]
    B6filt[shot,:] = data_pico_rest[19]

# We need the z component for position 5, r and theta were created above.

z = raw.create_dataset('z', data=Bdot3raw)
z = bdot.create_dataset('z', data=Bdot3raw)
z = b.create_dataset('z', data=B3)
z = Bfilt.create_dataset('z', data=B3filt)

# Done with position 5, move on to position 7

pos = h5f.create_group("pos7")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw)
theta = raw.create_dataset('theta', data=Bdot5raw)
z = raw.create_dataset('z', data=Bdot6raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4)
theta = bdot.create_dataset('theta', data=Bdot5)
z = bdot.create_dataset('z', data=Bdot6)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4)
theta = b.create_dataset('theta', data=B5)
z = b.create_dataset('z', data=B6)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt)
theta = Bfilt.create_dataset('theta', data=B5filt)
z = Bfilt.create_dataset('z', data=B6filt)


# Picoscope 3
for shot in np.arange(1,112):
    data_pico_rest = loadrest.load_picoscope_rest_01122022(shot, scopenum=3)

    Bdot3raw[shot,:] = data_pico_rest[4].T
    Bdot4raw[shot,:] = data_pico_rest[5].T
    Bdot5raw[shot,:] = data_pico_rest[6].T
    Bdot6raw[shot,:] = data_pico_rest[7].T
    Bdot3[shot,:] = data_pico_rest[8].T
    Bdot4[shot,:] = data_pico_rest[9].T
    Bdot5[shot,:] = data_pico_rest[10].T
    Bdot6[shot,:] = data_pico_rest[11].T
    B3[shot,:] = data_pico_rest[12]
    B4[shot,:] = data_pico_rest[13]
    B5[shot,:] = data_pico_rest[14]
    B6[shot,:] = data_pico_rest[15]
    B3filt[shot,:] = data_pico_rest[16]
    B4filt[shot,:] = data_pico_rest[17]
    B5filt[shot,:] = data_pico_rest[18]
    B6filt[shot,:] = data_pico_rest[19]


pos = h5f.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw)
theta = raw.create_dataset('theta', data=Bdot4raw)
z = raw.create_dataset('z', data=Bdot5raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3)
theta = bdot.create_dataset('theta', data=Bdot4)
z = bdot.create_dataset('z', data=Bdot5)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3)
theta = b.create_dataset('theta', data=B4)
z = b.create_dataset('z', data=B5)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt)
theta = Bfilt.create_dataset('theta', data=B4filt)
z = Bfilt.create_dataset('z', data=B5filt)


#Pico3 has positions 19 and 21 in it.
pos = h5f.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt)


# Picoscope 4

for shot in np.arange(1,112):
    data_pico_rest = loadrest.load_picoscope_rest_01122022(shot, scopenum=4)
    
    Bdot3raw[shot,:] = data_pico_rest[4].T
    Bdot4raw[shot,:] = data_pico_rest[5].T
    Bdot5raw[shot,:] = data_pico_rest[6].T
    Bdot6raw[shot,:] = data_pico_rest[7].T
    Bdot3[shot,:] = data_pico_rest[8].T
    Bdot4[shot,:] = data_pico_rest[9].T
    Bdot5[shot,:] = data_pico_rest[10].T
    Bdot6[shot,:] = data_pico_rest[11].T
    B3[shot,:] = data_pico_rest[12]
    B4[shot,:] = data_pico_rest[13]
    B5[shot,:] = data_pico_rest[14]
    B6[shot,:] = data_pico_rest[15]
    B3filt[shot,:] = data_pico_rest[16]
    B4filt[shot,:] = data_pico_rest[17]
    B5filt[shot,:] = data_pico_rest[18]
    B6filt[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw)
z = raw.create_dataset('z', data=Bdot4raw)

theta = bdot.create_dataset('theta', data=Bdot3)
z = bdot.create_dataset('z', data=Bdot4)

theta = b.create_dataset('theta', data=B3)
z = b.create_dataset('z', data=B4)

theta = Bfilt.create_dataset('theta', data=B3filt)
z = Bfilt.create_dataset('z', data=B4filt)

# New position: 33

pos = h5f.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw)
theta = raw.create_dataset('theta', data=Bdot6raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5)
theta = bdot.create_dataset('theta', data=Bdot6)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5)
theta = b.create_dataset('theta', data=B6)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt)
theta = Bfilt.create_dataset('theta', data=B6filt)


# Picoscope 5, start with position 33

for shot in np.arange(1,112):
    data_pico_rest = loadrest.load_picoscope_rest_01122022(shot, scopenum=5)

    Bdot3raw[shot,:] = data_pico_rest[4].T
    Bdot4raw[shot,:] = data_pico_rest[5].T
    Bdot5raw[shot,:] = data_pico_rest[6].T
    Bdot6raw[shot,:] = data_pico_rest[7].T
    Bdot3[shot,:] = data_pico_rest[8].T
    Bdot4[shot,:] = data_pico_rest[9].T
    Bdot5[shot,:] = data_pico_rest[10].T
    Bdot6[shot,:] = data_pico_rest[11].T
    B3[shot,:] = data_pico_rest[12]
    B4[shot,:] = data_pico_rest[13]
    B5[shot,:] = data_pico_rest[14]
    B6[shot,:] = data_pico_rest[15]
    B3filt[shot,:] = data_pico_rest[16]
    B4filt[shot,:] = data_pico_rest[17]
    B5filt[shot,:] = data_pico_rest[18]
    B6filt[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw)
z = bdot.create_dataset('z', data=Bdot3)
z = b.create_dataset('z', data=B3)
z = Bfilt.create_dataset('z', data=B3filt)

#Last position, 35:

pos = h5f.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw)
theta = raw.create_dataset('theta', data=Bdot5raw)
z = raw.create_dataset('z', data=Bdot6raw)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4)
theta = bdot.create_dataset('theta', data=Bdot5)
z = bdot.create_dataset('z', data=Bdot6)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4)
theta = b.create_dataset('theta', data=B5)
z = b.create_dataset('z', data=B6)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt)
theta = Bfilt.create_dataset('theta', data=B5filt)
z = Bfilt.create_dataset('z', data=B6filt)


channels = ['pos5','pos7','pos19','pos21','pos33','pos35']
for channel in channels:
    plt.plot(h5f[channel]['Bdot']['z'])
    
h5f.close()

