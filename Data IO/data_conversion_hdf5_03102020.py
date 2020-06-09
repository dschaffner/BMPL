# -*- coding: utf-8 -*-
#data_conversion_03102020.py
"""
Created on Tue 09 2020

@author: ccartagena
"""

import load_picoscope_bmx_03102020 as load
import numpy as np
import h5py
import matplotlib.pylab as plt

#####################################################################################
"""
    The purpose of this script is to generated a HDF5 file from the 03102020 dataset.
"""
#####################################################################################

###################################################
""" The picoscope-probe structure 03102020 """
###################################################
"""          A           B           C           D
pico1       1R          1T          1Z          DC
pico2       3R          3T          3Z          --
pico3       5R          5T          5Z          --
pico4       7R          7T          7Z          --
pico5       9R          9T          9Z          --
pico6       11R         11T         11Z         --
pico7       13R         13T         13Z         --
pico8       15R         15T         15Z         HV

N --- Density
DC --- Discharge Current
HV --- High Voltage
Total Shots = 25
"""
###################################################

### The filepath and filename lines will need to be changed to the location and
### filename of the data in your system.
#filepath = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\10142019\\processed\\'
#filename='2kV_LangPort4_0p5ms_1msGas_8shots_bdot_10142019'
filepath = '/Volumes/CarFlor/Research/Data/2020/03102020/HDF5/'
#filename='2kV_0p5ms_1msGas_20shots_bdot_10142019' #shots 7 to 26
filename = '2kV_1p5ms_9msGas_25shots_bdot_03102020' #27 to 47 skip 30

##### HDF5 file creation and structure setup
h5f = h5py.File(filepath + filename + '.h5', 'a')

probe_dia = 0.00158755      #m (1/16'' probe)
hole_sep = 0.001016     #m (1/16''probe)
probe_info=h5f.create_group("probe_info")
probe_info.create_dataset('probe_stalk_diameter', data = probe_dia)
probe_info.create_dataset('probe_hole_separation', data = hole_sep)
probe_info.create_dataset('rloop_area', data = np.pi*(probe_dia/2)**2)
probe_info.create_dataset('tloop_area', data = probe_dia*hole_sep)
probe_info.create_dataset('zloop_area', data = probe_dia*hole_sep)
probe_info.create_dataset('radial_location', data = 'center')

run_info = h5f.create_group("run_info")
run_info.create_dataset('Discharge_Voltage (kV)', data = 2e3)
run_info.create_dataset('Collection Dates', data = '03102020')
run_info.create_dataset('Stuffing Delay (ms)', data = 1.5)

##########################################################################
# Output structure of load_picoscope()
#:time_ms,time_s,timeB_s,timeB_ms,Bdot1,Bdot2,Bdot3,B1,B2,B3,B1filt,B2filt,B3filt,Bdotraw1,Bdotraw2,Bdotraw3,timeraw
##########################################################################
##### Obtaining temporal information time_ms, time_s, timeB_s, and timeraw
data = load.load_picoscope(1, scopenum = 1)
shotnum = 25
time_us, time_s, timeB_s, timeB_us, timeraw = data[0], data[1], data[2], data[3], data[16]
# Creation of time HDF5 group
time = h5f.create_group("time")
time.create_dataset('time_us', data = time_us)
time.create_dataset('time_s', data = time_s)
time.create_dataset('timeB_us', data = timeB_us)
time.create_dataset('timeraw', data = timeraw)

datalength = data[0].shape[0]
datalengthraw = data[16].shape[0]
###### The following initilization makes it easier to tweek the code for different datasets.
Bdot1raw = np.zeros([shotnum, datalengthraw])
Bdot2raw = np.zeros([shotnum, datalengthraw])
Bdot3raw = np.zeros([shotnum, datalengthraw])
#Bdot4raw = np.zeros([shotnum, datalengthraw])


Bdot1 = np.zeros([shotnum, datalength])
Bdot2 = np.zeros([shotnum, datalength])
Bdot3 = np.zeros([shotnum, datalength])
#Bdot4 = np.zeros([shotnum, datalength])


B1 = np.zeros([shotnum, datalength - 1])
B2 = np.zeros([shotnum, datalength - 1])
B3 = np.zeros([shotnum, datalength - 1])
#B4 = np.zeros([shotnum, datalength - 1])


#Picoscope1
for shot in np.arange(1, shotnum + 1):
    
    data = load.load_picoscope(shot, scopenum = 1)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 1 HDF5 group
pos = h5f.create_group("pos1")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)


#Picoscope2
for shot in np.arange(1, shotnum + 1):
    
    data  = load.load_picoscope(shot, scopenum =  2)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 3 HDF5 group
pos = h5f.create_group("pos3")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

#Picoscope3
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 3)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 5 HDF5 group
pos = h5f.create_group("pos5")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)


#Picoscope4
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 4)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 7 HDF5 group
pos = h5f.create_group("pos7")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

#Picoscope5
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 5)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 9 HDF5 group
pos = h5f.create_group("pos9")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

#Picoscope6
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 6)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 11 HDF5 group
pos = h5f.create_group("pos11")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

#Picoscope7
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 7)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 13 HDF5 group
pos = h5f.create_group("pos13")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

#Picoscope8
for shot in np.arange(1, 26):
    
    data = load.load_picoscope(shot, scopenum = 8)
    
    Bdot1raw[shot - 1, :] = data[13]
    Bdot2raw[shot - 1, :] = data[14]
    Bdot3raw[shot - 1, :] = data[15]
    #Bdot4raw[shot - 1, :] = data[19]
    
    Bdot1[shot - 1, :] = data[4]
    Bdot2[shot - 1, :] = data[5]
    Bdot3[shot - 1, :] = data[6]
    #Bdot4[shot - 1, :] = data[7]
    
    B1[shot - 1, :] = data[7]
    B2[shot - 1, :] = data[8]
    B3[shot - 1, :] = data[9]
    #B4[shot - 1, :] = data[11]

#Creation of probe 15 HDF5 group
pos = h5f.create_group("pos15")
raw = pos.create_group("raw")
##?
r = raw.create_dataset('r', data = Bdot1raw)
theta = raw.create_dataset('theta', data = Bdot2raw)
z = raw.create_dataset('z', data = Bdot3raw)
##?
bdot = pos.create_group("bdot")
r = bdot.create_dataset('r', data = Bdot1)
theta = bdot.create_dataset('theta', data = Bdot2)
z = bdot.create_dataset('z', data = Bdot3)
##?
b = pos.create_group("b")
r = b.create_dataset('r', data = B1)
theta = b.create_dataset('theta', data = B2)
z = b.create_dataset('z', data = B3)

channels = ['pos1', 'pos3', 'pos5', 'pos7', 'pos9', 'pos11', 'pos13', 'pos15']
#for shot in np.arange(shotnum):
for channel in channels:
    plt.plot(h5f[channel]['b']['z'][0, :])

