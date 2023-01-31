# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:26:37 2023

@author: Josh0


"""
#####################################################################################
"""
    This code will take the data taken since the BMPL move to the basement. Start with January 2022 until September 2022, all data taken in 2022.
"""
#####################################################################################


import numpy as np
import h5py
import matplotlib.pylab as plt

import load_picoscope_one_01122022 as loadone_112
import load_picoscope_one_02232022 as loadone_223
import load_picoscope_one_03172022 as loadone_317
import load_picoscope_one_06202022 as loadone_620
import load_picoscope_one_07152022 as loadone_715
import load_picoscope_one_07202022 as loadone_720
import load_picoscope_one_09232022 as loadone_923


import load_picoscope_rest_01122022 as loadrest_112
import load_picoscope_rest_02232022 as loadrest_223
import load_picoscope_rest_03172022 as loadrest_317
import load_picoscope_rest_06202022 as loadrest_620
import load_picoscope_rest_07152022 as loadrest_715
import load_picoscope_rest_07202022 as loadrest_720
import load_picoscope_rest_09232022 as loadrest_923

#File Target of Data Conversion
new_filepath = 'C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
new_filename = '2022_FullCatalogueofShots'
h5f = h5py.File(new_filepath + new_filename + '.h5', 'a')

#For Velocity Data
timeavg_data_directory_location = 'C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr COllege\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\Time_avgd_velocities\\'

data112 = h5f.create_group('01122022 Dataset')
data223 = h5f.create_group('02232022 Dataset')
data317 = h5f.create_group('03172022 Dataset')
data620 = h5f.create_group('06202022 Dataset')
data715 = h5f.create_group('07152022 Dataset')
data720 = h5f.create_group('07202022 Dataset')
data923 = h5f.create_group('09232022 Dataset')
datalength = int(25000)


Marked_Line_for_1122022 = np.zeros([1])
##########################
""" For dataset 1122022"""
##########################

probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info112 = data112.create_group('Probe Info')
probe_info112.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info112.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info112.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info112.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info112.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info112.create_dataset('radial_location', data = 'r=0cm')
probe_info112.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

run_info112 = data112.create_group('Run Info')
run_info112.create_dataset('Discharge Voltage (kV)', data=2.5)
run_info112.create_dataset('Stuffing Delay (ms)', data=1.5)

data112_pico1 = loadone_112.load_picoscope_one_01122022(1, scopenum=1)
data112_pico_rest = loadrest_112.load_picoscope_rest_01122022(1, scopenum=1)
shotnum112 = 112
time_s, time_us, timeB_s, timeB_us = data112_pico1[0],data112_pico1[1],data112_pico1[2],data112_pico1[3]

time112 = data112.create_group("Time")
time112.create_dataset('time_s', data=time_s)
time112.create_dataset('time_us', data = time_us)
time112.create_dataset('timeB_s', data = timeB_s)
time112.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw112 = np.zeros([shotnum112,datalength])
HighV_raw112 = np.zeros([shotnum112,datalength])
DisCurrent112 = np.zeros([shotnum112,datalength])
HighV112 = np.zeros([shotnum112,datalength])
Bdot1raw112 = np.zeros([shotnum112,datalength])
Bdot2raw112 = np.zeros([shotnum112,datalength])
Bdot3raw112 = np.zeros([shotnum112,datalength])
Bdot4raw112 = np.zeros([shotnum112,datalength])
Bdot5raw112 = np.zeros([shotnum112,datalength])
Bdot6raw112 = np.zeros([shotnum112,datalength])
Bdot1112 = np.zeros([shotnum112,datalength])
Bdot2112 = np.zeros([shotnum112,datalength])
Bdot3112 = np.zeros([shotnum112,datalength])
Bdot4112 = np.zeros([shotnum112,datalength])
Bdot5112 = np.zeros([shotnum112,datalength])
Bdot6112 = np.zeros([shotnum112,datalength])
B1112 = np.zeros([shotnum112,datalength-1])
B2112 = np.zeros([shotnum112,datalength-1])
B3112 = np.zeros([shotnum112,datalength-1])
B4112 = np.zeros([shotnum112,datalength-1])
B5112 = np.zeros([shotnum112,datalength-1])
B6112 = np.zeros([shotnum112,datalength-1])
B1filt112 = np.zeros([shotnum112,datalength-1])
B2filt112 = np.zeros([shotnum112,datalength-1])
B3filt112 = np.zeros([shotnum112,datalength-1])
B4filt112 = np.zeros([shotnum112,datalength-1])
B5filt112 = np.zeros([shotnum112,datalength-1])
B6filt112 = np.zeros([shotnum112,datalength-1])

# Picoscope 1, will only use data_pico1
for shot in np.arange (1, shotnum112):
    data_pico1 = loadone_112.load_picoscope_one_01122022(shot, scopenum=1)
    
    DisCurrent_raw112[shot,:] = data_pico1[4].T
    HighV_raw112[shot,:] = data_pico1[5].T
    DisCurrent112[shot,:] = data_pico1[6]
    HighV112[shot,:] = data_pico1[7]
    Bdot1raw112[shot,:] = data_pico1[8].T
    Bdot2raw112[shot,:] = data_pico1[9].T
    Bdot1112[shot,:] = data_pico1[10].T
    Bdot2112[shot,:] = data_pico1[11].T
    B1112[shot,:] = data_pico1[12]
    B2112[shot,:] = data_pico1[13]
    B1filt112[shot,:] = data_pico1[14]
    B2filt112[shot,:] = data_pico1[15]

dis112 = data112.create_group("Discharge Data")
disraw112 = dis112.create_group("Discharge Raw")
DisCurrent_raw = disraw112.create_dataset("DisCurrent_raw", data=DisCurrent_raw112)
HighV_raw = disraw112.create_dataset("HighV_raw", data=HighV_raw112)
discharge112 = dis112.create_group("Discharge")
DisCurrent = discharge112.create_dataset("DisCurrent", data=DisCurrent112)
HighV = discharge112.create_dataset("HighV", data=HighV112)

pos112 = data112.create_group("pos5")
raw112 = pos112.create_group("raw")
r = raw112.create_dataset('r', data=Bdot1raw112)
theta = raw112.create_dataset('theta', data=Bdot2raw112)

bdot112 = pos112.create_group('Bdot')
r = bdot112.create_dataset('r', data=Bdot1112)
theta = bdot112.create_dataset('theta', data=Bdot2112)

b112 = pos112.create_group('B')
r = b112.create_dataset('r', data=B1112)
theta = b112.create_dataset('theta', data=B2112)

Bfilt112 = pos112.create_group('Bfilt')
r = Bfilt112.create_dataset('r', data=B1filt112)
theta = Bfilt112.create_dataset('theta', data=B2filt112)

# Picoscope 2
for shot in np.arange(1, shotnum112):
    #if shot==45:
    #    continue
    data_pico_rest = loadrest_112.load_picoscope_rest_01122022(shot, scopenum=2)
    
    Bdot3raw112[shot,:] = data_pico_rest[4].T
    Bdot4raw112[shot,:] = data_pico_rest[5].T
    Bdot5raw112[shot,:] = data_pico_rest[6].T
    Bdot6raw112[shot,:] = data_pico_rest[7].T
    Bdot3112[shot,:] = data_pico_rest[8].T
    Bdot4112[shot,:] = data_pico_rest[9].T
    Bdot5112[shot,:] = data_pico_rest[10].T
    Bdot6112[shot,:] = data_pico_rest[11].T
    B3112[shot,:] = data_pico_rest[12]
    B4112[shot,:] = data_pico_rest[13]
    B5112[shot,:] = data_pico_rest[14]
    B6112[shot,:] = data_pico_rest[15]
    B3filt112[shot,:] = data_pico_rest[16]
    B4filt112[shot,:] = data_pico_rest[17]
    B5filt112[shot,:] = data_pico_rest[18]
    B6filt112[shot,:] = data_pico_rest[19]

z = raw112.create_dataset('z', data=Bdot3raw112)
z = bdot112.create_dataset('z', data=Bdot3112)
z = b112.create_dataset('z', data=B3112)
z = Bfilt112.create_dataset('z', data=B3filt112)

#I don't think the 112 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data112.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw112)
theta = raw.create_dataset('theta', data=Bdot5raw112)
z = raw.create_dataset('z', data=Bdot6raw112)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4112)
theta = bdot.create_dataset('theta', data=Bdot5112)
z = bdot.create_dataset('z', data=Bdot6112)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4112)
theta = b.create_dataset('theta', data=B5112)
z = b.create_dataset('z', data=B6112)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt112)
theta = bfilt.create_dataset('theta', data=B5filt112)
z = bfilt.create_dataset('z', data=B6filt112)

# Picoscope 3
for shot in np.arange(1, shotnum112):
    data_pico_rest = loadrest_112.load_picoscope_rest_01122022(shot, scopenum=3)

    Bdot3raw112[shot,:] = data_pico_rest[4].T
    Bdot4raw112[shot,:] = data_pico_rest[5].T
    Bdot5raw112[shot,:] = data_pico_rest[6].T
    Bdot6raw112[shot,:] = data_pico_rest[7].T
    Bdot3112[shot,:] = data_pico_rest[8].T
    Bdot4112[shot,:] = data_pico_rest[9].T
    Bdot5112[shot,:] = data_pico_rest[10].T
    Bdot6112[shot,:] = data_pico_rest[11].T
    B3112[shot,:] = data_pico_rest[12]
    B4112[shot,:] = data_pico_rest[13]
    B5112[shot,:] = data_pico_rest[14]
    B6112[shot,:] = data_pico_rest[15]
    B3filt112[shot,:] = data_pico_rest[16]
    B4filt112[shot,:] = data_pico_rest[17]
    B5filt112[shot,:] = data_pico_rest[18]
    B6filt112[shot,:] = data_pico_rest[19]

pos = data112.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw112)
theta = raw.create_dataset('theta', data=Bdot4raw112)
z = raw.create_dataset('z', data=Bdot5raw112)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3112)
theta = bdot.create_dataset('theta', data=Bdot4112)
z = bdot.create_dataset('z', data=Bdot5112)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3112)
theta = b.create_dataset('theta', data=B4112)
z = b.create_dataset('z', data=B5112)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt112)
theta = Bfilt.create_dataset('theta', data=B4filt112)
z = Bfilt.create_dataset('z', data=B5filt112)

#Pico3 has positions 19 and 21 in it.
pos = data112.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw112)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6112)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6112)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt112)

# Picoscope 4
for shot in np.arange(1,shotnum112):
    data_pico_rest = loadrest_112.load_picoscope_rest_01122022(shot, scopenum=4)
    
    Bdot3raw112[shot,:] = data_pico_rest[4].T
    Bdot4raw112[shot,:] = data_pico_rest[5].T
    Bdot5raw112[shot,:] = data_pico_rest[6].T
    Bdot6raw112[shot,:] = data_pico_rest[7].T
    Bdot3112[shot,:] = data_pico_rest[8].T
    Bdot4112[shot,:] = data_pico_rest[9].T
    Bdot5112[shot,:] = data_pico_rest[10].T
    Bdot6112[shot,:] = data_pico_rest[11].T
    B3112[shot,:] = data_pico_rest[12]
    B4112[shot,:] = data_pico_rest[13]
    B5112[shot,:] = data_pico_rest[14]
    B6112[shot,:] = data_pico_rest[15]
    B3filt112[shot,:] = data_pico_rest[16]
    B4filt112[shot,:] = data_pico_rest[17]
    B5filt112[shot,:] = data_pico_rest[18]
    B6filt112[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw112)
z = raw.create_dataset('z', data=Bdot4raw112)

theta = bdot.create_dataset('theta', data=Bdot3112)
z = bdot.create_dataset('z', data=Bdot4112)

theta = b.create_dataset('theta', data=B3112)
z = b.create_dataset('z', data=B4112)

theta = Bfilt.create_dataset('theta', data=B3filt112)
z = Bfilt.create_dataset('z', data=B4filt112)

# New position: 33
pos = data112.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw112)
theta = raw.create_dataset('theta', data=Bdot6raw112)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5112)
theta = bdot.create_dataset('theta', data=Bdot6112)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5112)
theta = b.create_dataset('theta', data=B6112)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt112)
theta = Bfilt.create_dataset('theta', data=B6filt112)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum112):
    data_pico_rest = loadrest_112.load_picoscope_rest_01122022(shot, scopenum=5)

    Bdot3raw112[shot,:] = data_pico_rest[4].T
    Bdot4raw112[shot,:] = data_pico_rest[5].T
    Bdot5raw112[shot,:] = data_pico_rest[6].T
    Bdot6raw112[shot,:] = data_pico_rest[7].T
    Bdot3112[shot,:] = data_pico_rest[8].T
    Bdot4112[shot,:] = data_pico_rest[9].T
    Bdot5112[shot,:] = data_pico_rest[10].T
    Bdot6112[shot,:] = data_pico_rest[11].T
    B3112[shot,:] = data_pico_rest[12]
    B4112[shot,:] = data_pico_rest[13]
    B5112[shot,:] = data_pico_rest[14]
    B6112[shot,:] = data_pico_rest[15]
    B3filt112[shot,:] = data_pico_rest[16]
    B4filt112[shot,:] = data_pico_rest[17]
    B5filt112[shot,:] = data_pico_rest[18]
    B6filt112[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw112)
z = bdot.create_dataset('z', data=Bdot3112)
z = b.create_dataset('z', data=B3112)
z = Bfilt.create_dataset('z', data=B3filt112)

#Last position, 35:
pos = data112.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw112)
theta = raw.create_dataset('theta', data=Bdot5raw112)
z = raw.create_dataset('z', data=Bdot6raw112)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4112)
theta = bdot.create_dataset('theta', data=Bdot5112)
z = bdot.create_dataset('z', data=Bdot6112)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4112)
theta = b.create_dataset('theta', data=B5112)
z = b.create_dataset('z', data=B6112)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt112)
theta = Bfilt.create_dataset('theta', data=B5filt112)
z = Bfilt.create_dataset('z', data=B6filt112)


Marked_Line_for_2232022 = np.zeros([1])
###########################
""" For dataset 2232022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info223 = data223.create_group('Probe Info')
probe_info223.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info223.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info223.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info223.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info223.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info223.create_dataset('radial_location', data = 'r=0cm')
probe_info223.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

run_info223 = data223.create_group('Run Info')
run_info223.create_dataset('Discharge Voltage (kV)', data=2.5)
run_info223.create_dataset('Stuffing Delay (ms)', data=1.5)

data223_pico1 = loadone_223.load_picoscope_one_02232022(1, scopenum=1)
data223_pico_rest = loadrest_223.load_picoscope_rest_02232022(1, scopenum=2)
shotnum223 = 22
time_s, time_us, timeB_s, timeB_us = data223_pico1[0],data223_pico1[1],data223_pico1[2],data223_pico1[3]

time223 = data223.create_group("Time")
time223.create_dataset('time_s', data=time_s)
time223.create_dataset('time_us', data = time_us)
time223.create_dataset('timeB_s', data = timeB_s)
time223.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw223 = np.zeros([shotnum223,datalength])
HighV_raw223 = np.zeros([shotnum223,datalength])
DisCurrent223 = np.zeros([shotnum223,datalength])
HighV223 = np.zeros([shotnum223,datalength])
Bdot1raw223 = np.zeros([shotnum223,datalength])
Bdot2raw223 = np.zeros([shotnum223,datalength])
Bdot3raw223 = np.zeros([shotnum223,datalength])
Bdot4raw223 = np.zeros([shotnum223,datalength])
Bdot5raw223 = np.zeros([shotnum223,datalength])
Bdot6raw223 = np.zeros([shotnum223,datalength])
Bdot1223 = np.zeros([shotnum223,datalength])
Bdot2223 = np.zeros([shotnum223,datalength])
Bdot3223 = np.zeros([shotnum223,datalength])
Bdot4223 = np.zeros([shotnum223,datalength])
Bdot5223 = np.zeros([shotnum223,datalength])
Bdot6223 = np.zeros([shotnum223,datalength])
B1223 = np.zeros([shotnum223,datalength-1])
B2223 = np.zeros([shotnum223,datalength-1])
B3223 = np.zeros([shotnum223,datalength-1])
B4223 = np.zeros([shotnum223,datalength-1])
B5223 = np.zeros([shotnum223,datalength-1])
B6223 = np.zeros([shotnum223,datalength-1])
B1filt223 = np.zeros([shotnum223,datalength-1])
B2filt223 = np.zeros([shotnum223,datalength-1])
B3filt223 = np.zeros([shotnum223,datalength-1])
B4filt223 = np.zeros([shotnum223,datalength-1])
B5filt223 = np.zeros([shotnum223,datalength-1])
B6filt223 = np.zeros([shotnum223,datalength-1])

# Picoscope 1
for shot in np.arange (1, shotnum223):
    data223_pico1 = loadone_223.load_picoscope_one_02232022(shot, scopenum=1)
    
    DisCurrent_raw223[shot,:] = data_pico1[4].T
    HighV_raw223[shot,:] = data_pico1[5].T
    DisCurrent223[shot,:] = data_pico1[6]
    HighV223[shot,:] = data_pico1[7]
    Bdot1raw223[shot,:] = data_pico1[8].T
    Bdot2raw223[shot,:] = data_pico1[9].T
    Bdot1223[shot,:] = data_pico1[10].T
    Bdot2223[shot,:] = data_pico1[11].T
    B1223[shot,:] = data_pico1[12]
    B2223[shot,:] = data_pico1[13]
    B1filt223[shot,:] = data_pico1[14]
    B2filt223[shot,:] = data_pico1[15]

dis223 = data223.create_group("Discharge Data")
disraw223 = dis223.create_group("Discharge Raw")
DisCurrent_raw = disraw223.create_dataset("DisCurrent_raw", data=DisCurrent_raw223)
HighV_raw = disraw223.create_dataset("HighV_raw", data=HighV_raw223)
discharge223 = dis223.create_group("Discharge")
DisCurrent = discharge223.create_dataset("DisCurrent", data=DisCurrent223)
HighV = discharge223.create_dataset("HighV", data=HighV223)

pos = data223.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw223)
theta = raw.create_dataset('theta', data=Bdot2raw223)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1223)
theta = bdot.create_dataset('theta', data=Bdot2223)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1223)
theta = b.create_dataset('theta', data=B2223)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt223)
theta = Bfilt.create_dataset('theta', data=B2filt223)

# Picoscope 2
for shot in np.arange(1, shotnum223):
    #if shot==45:
    #    continue
    data_pico_rest = loadrest_223.load_picoscope_rest_02232022(shot, scopenum=2)
    
    Bdot3raw223[shot,:] = data_pico_rest[4].T
    Bdot4raw223[shot,:] = data_pico_rest[5].T
    Bdot5raw223[shot,:] = data_pico_rest[6].T
    Bdot6raw223[shot,:] = data_pico_rest[7].T
    Bdot3223[shot,:] = data_pico_rest[8].T
    Bdot4223[shot,:] = data_pico_rest[9].T
    Bdot5223[shot,:] = data_pico_rest[10].T
    Bdot6223[shot,:] = data_pico_rest[11].T
    B3223[shot,:] = data_pico_rest[12]
    B4223[shot,:] = data_pico_rest[13]
    B5223[shot,:] = data_pico_rest[14]
    B6223[shot,:] = data_pico_rest[15]
    B3filt223[shot,:] = data_pico_rest[16]
    B4filt223[shot,:] = data_pico_rest[17]
    B5filt223[shot,:] = data_pico_rest[18]
    B6filt223[shot,:] = data_pico_rest[19]

z = raw.create_dataset('z', data=Bdot3raw223)
z = bdot.create_dataset('z', data=Bdot3223)
z = b.create_dataset('z', data=B3223)
z = Bfilt.create_dataset('z', data=B3filt223)

#I don't think the 223 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data223.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw223)
theta = raw.create_dataset('theta', data=Bdot5raw223)
z = raw.create_dataset('z', data=Bdot6raw223)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4223)
theta = bdot.create_dataset('theta', data=Bdot5223)
z = bdot.create_dataset('z', data=Bdot6223)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4223)
theta = b.create_dataset('theta', data=B5223)
z = b.create_dataset('z', data=B6223)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt223)
theta = bfilt.create_dataset('theta', data=B5filt223)
z = bfilt.create_dataset('z', data=B6filt223)


# Picoscope 3
for shot in np.arange(1, shotnum223):
    data_pico_rest = loadrest_223.load_picoscope_rest_02232022(shot, scopenum=3)

    Bdot3raw223[shot,:] = data_pico_rest[4].T
    Bdot4raw223[shot,:] = data_pico_rest[5].T
    Bdot5raw223[shot,:] = data_pico_rest[6].T
    Bdot6raw223[shot,:] = data_pico_rest[7].T
    Bdot3223[shot,:] = data_pico_rest[8].T
    Bdot4223[shot,:] = data_pico_rest[9].T
    Bdot5223[shot,:] = data_pico_rest[10].T
    Bdot6223[shot,:] = data_pico_rest[11].T
    B3223[shot,:] = data_pico_rest[12]
    B4223[shot,:] = data_pico_rest[13]
    B5223[shot,:] = data_pico_rest[14]
    B6223[shot,:] = data_pico_rest[15]
    B3filt223[shot,:] = data_pico_rest[16]
    B4filt223[shot,:] = data_pico_rest[17]
    B5filt223[shot,:] = data_pico_rest[18]
    B6filt223[shot,:] = data_pico_rest[19]

pos = data223.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw223)
theta = raw.create_dataset('theta', data=Bdot4raw223)
z = raw.create_dataset('z', data=Bdot5raw223)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3223)
theta = bdot.create_dataset('theta', data=Bdot4223)
z = bdot.create_dataset('z', data=Bdot5223)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3223)
theta = b.create_dataset('theta', data=B4223)
z = b.create_dataset('z', data=B5223)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt223)
theta = Bfilt.create_dataset('theta', data=B4filt223)
z = Bfilt.create_dataset('z', data=B5filt223)

#Pico3 has positions 19 and 21 in it.
pos = data223.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw223)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6223)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6223)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt223)

# Picoscope 4
for shot in np.arange(1,shotnum223):
    data_pico_rest = loadrest_223.load_picoscope_rest_02232022(shot, scopenum=4)
    
    Bdot3raw223[shot,:] = data_pico_rest[4].T
    Bdot4raw223[shot,:] = data_pico_rest[5].T
    Bdot5raw223[shot,:] = data_pico_rest[6].T
    Bdot6raw223[shot,:] = data_pico_rest[7].T
    Bdot3223[shot,:] = data_pico_rest[8].T
    Bdot4223[shot,:] = data_pico_rest[9].T
    Bdot5223[shot,:] = data_pico_rest[10].T
    Bdot6223[shot,:] = data_pico_rest[11].T
    B3223[shot,:] = data_pico_rest[12]
    B4223[shot,:] = data_pico_rest[13]
    B5223[shot,:] = data_pico_rest[14]
    B6223[shot,:] = data_pico_rest[15]
    B3filt223[shot,:] = data_pico_rest[16]
    B4filt223[shot,:] = data_pico_rest[17]
    B5filt223[shot,:] = data_pico_rest[18]
    B6filt223[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw223)
z = raw.create_dataset('z', data=Bdot4raw223)

theta = bdot.create_dataset('theta', data=Bdot3223)
z = bdot.create_dataset('z', data=Bdot4223)

theta = b.create_dataset('theta', data=B3223)
z = b.create_dataset('z', data=B4223)

theta = Bfilt.create_dataset('theta', data=B3filt223)
z = Bfilt.create_dataset('z', data=B4filt223)

# New position: 33
pos = data223.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw223)
theta = raw.create_dataset('theta', data=Bdot6raw223)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5223)
theta = bdot.create_dataset('theta', data=Bdot6223)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5223)
theta = b.create_dataset('theta', data=B6223)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt223)
theta = Bfilt.create_dataset('theta', data=B6filt223)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum223):
    data_pico_rest = loadrest_223.load_picoscope_rest_02232022(shot, scopenum=5)

    Bdot3raw223[shot,:] = data_pico_rest[4].T
    Bdot4raw223[shot,:] = data_pico_rest[5].T
    Bdot5raw223[shot,:] = data_pico_rest[6].T
    Bdot6raw223[shot,:] = data_pico_rest[7].T
    Bdot3223[shot,:] = data_pico_rest[8].T
    Bdot4223[shot,:] = data_pico_rest[9].T
    Bdot5223[shot,:] = data_pico_rest[10].T
    Bdot6223[shot,:] = data_pico_rest[11].T
    B3223[shot,:] = data_pico_rest[12]
    B4223[shot,:] = data_pico_rest[13]
    B5223[shot,:] = data_pico_rest[14]
    B6223[shot,:] = data_pico_rest[15]
    B3filt223[shot,:] = data_pico_rest[16]
    B4filt223[shot,:] = data_pico_rest[17]
    B5filt223[shot,:] = data_pico_rest[18]
    B6filt223[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw223)
z = bdot.create_dataset('z', data=Bdot3223)
z = b.create_dataset('z', data=B3223)
z = Bfilt.create_dataset('z', data=B3filt223)

#Last position, 35:
pos = data223.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw223)
theta = raw.create_dataset('theta', data=Bdot5raw223)
z = raw.create_dataset('z', data=Bdot6raw223)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4223)
theta = bdot.create_dataset('theta', data=Bdot5223)
z = bdot.create_dataset('z', data=Bdot6223)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4223)
theta = b.create_dataset('theta', data=B5223)
z = b.create_dataset('z', data=B6223)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt223)
theta = Bfilt.create_dataset('theta', data=B5filt223)
z = Bfilt.create_dataset('z', data=B6filt223)


Marked_Line_for_3172022 = np.zeros([1])
###########################
""" For dataset 3172022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info317 = data317.create_group('Probe Info')
probe_info317.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info317.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info317.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info317.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info317.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info317.create_dataset('radial_location', data = 'r=0cm')
probe_info317.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

run_info317 = data317.create_group('Run Info')
run_info317.create_dataset('Discharge Voltage (kV)', data=2.5)
run_info317.create_dataset('Stuffing Delay (ms)', data=1.5)
run_info317.create_dataset('Nozzle Fields', data='0mT Nozzle shots 1-52, 100mT Nozzle Shots 53-180, ~50mT Nozle Shots 181-253')

data317_pico1 = loadone_317.load_picoscope_one_03172022(1, scopenum=1)
data317_pico_rest = loadrest_317.load_picoscope_rest_03172022(1, scopenum=2)
shotnum317 = 254
time_s, time_us, timeB_s, timeB_us = data317_pico1[0],data317_pico1[1],data317_pico1[2],data317_pico1[3]

time317 = data317.create_group("Time")
time317.create_dataset('time_s', data=time_s)
time317.create_dataset('time_us', data = time_us)
time317.create_dataset('timeB_s', data = timeB_s)
time317.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw317 = np.zeros([shotnum317,datalength])
HighV_raw317 = np.zeros([shotnum317,datalength])
DisCurrent317 = np.zeros([shotnum317,datalength])
HighV317 = np.zeros([shotnum317,datalength])
Bdot1raw317 = np.zeros([shotnum317,datalength])
Bdot2raw317 = np.zeros([shotnum317,datalength])
Bdot3raw317 = np.zeros([shotnum317,datalength])
Bdot4raw317 = np.zeros([shotnum317,datalength])
Bdot5raw317 = np.zeros([shotnum317,datalength])
Bdot6raw317 = np.zeros([shotnum317,datalength])
Bdot1317 = np.zeros([shotnum317,datalength])
Bdot2317 = np.zeros([shotnum317,datalength])
Bdot3317 = np.zeros([shotnum317,datalength])
Bdot4317 = np.zeros([shotnum317,datalength])
Bdot5317 = np.zeros([shotnum317,datalength])
Bdot6317 = np.zeros([shotnum317,datalength])
B1317 = np.zeros([shotnum317,datalength-1])
B2317 = np.zeros([shotnum317,datalength-1])
B3317 = np.zeros([shotnum317,datalength-1])
B4317 = np.zeros([shotnum317,datalength-1])
B5317 = np.zeros([shotnum317,datalength-1])
B6317 = np.zeros([shotnum317,datalength-1])
B1filt317 = np.zeros([shotnum317,datalength-1])
B2filt317 = np.zeros([shotnum317,datalength-1])
B3filt317 = np.zeros([shotnum317,datalength-1])
B4filt317 = np.zeros([shotnum317,datalength-1])
B5filt317 = np.zeros([shotnum317,datalength-1])
B6filt317 = np.zeros([shotnum317,datalength-1])

#Time-Averaged Velocity Data
tavg_velocities_p57_3172022 = np.loadtxt(timeavg_data_directory_location+'Full3172022_ShotsandVelocities.txt', skiprows=1, unpack=True).T
vel317 = data317.create_group("Time-Averaged Velocities")
vel317.create_dataset('Data Key', data=r'Pos5Pos7[0] = Shot#, P5P7[1] = Delay (${\mu}$s), [2] = Velocity (km/s), [3] = Nozzle Strength (mT), [4] = Discharge Voltage (kV)')
vel317.create_dataset('Pos5Pos7', data=tavg_velocities_p57_3172022)

# Picoscope 1
for shot in np.arange (1, shotnum317):
    if shot==158:
        continue
    data317_pico1 = loadone_317.load_picoscope_one_03172022(shot, scopenum=1)
    DisCurrent_raw317[shot,:] = data_pico1[4].T
    HighV_raw317[shot,:] = data_pico1[5].T
    DisCurrent317[shot,:] = data_pico1[6]
    HighV317[shot,:] = data_pico1[7]
    Bdot1raw317[shot,:] = data_pico1[8].T
    Bdot2raw317[shot,:] = data_pico1[9].T
    Bdot1317[shot,:] = data_pico1[10].T
    Bdot2317[shot,:] = data_pico1[11].T
    B1317[shot,:] = data_pico1[12]
    B2317[shot,:] = data_pico1[13]
    B1filt317[shot,:] = data_pico1[14]
    B2filt317[shot,:] = data_pico1[15]

dis317 = data317.create_group("Discharge Data")
disraw317 = dis317.create_group("Discharge Raw")
DisCurrent_raw = disraw317.create_dataset("DisCurrent_raw", data=DisCurrent_raw317)
HighV_raw = disraw317.create_dataset("HighV_raw", data=HighV_raw317)
discharge317 = dis317.create_group("Discharge")
DisCurrent = discharge317.create_dataset("DisCurrent", data=DisCurrent317)
HighV = discharge317.create_dataset("HighV", data=HighV317)

pos = data317.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw317)
theta = raw.create_dataset('theta', data=Bdot2raw317)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1317)
theta = bdot.create_dataset('theta', data=Bdot2317)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1317)
theta = b.create_dataset('theta', data=B2317)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt317)
theta = Bfilt.create_dataset('theta', data=B2filt317)

# Picoscope 2
for shot in np.arange(1, shotnum317):
    if shot==158:
        continue
    data_pico_rest = loadrest_317.load_picoscope_rest_03172022(shot, scopenum=2)
    
    Bdot3raw317[shot,:] = data_pico_rest[4].T
    Bdot4raw317[shot,:] = data_pico_rest[5].T
    Bdot5raw317[shot,:] = data_pico_rest[6].T
    Bdot6raw317[shot,:] = data_pico_rest[7].T
    Bdot3317[shot,:] = data_pico_rest[8].T
    Bdot4317[shot,:] = data_pico_rest[9].T
    Bdot5317[shot,:] = data_pico_rest[10].T
    Bdot6317[shot,:] = data_pico_rest[11].T
    B3317[shot,:] = data_pico_rest[12]
    B4317[shot,:] = data_pico_rest[13]
    B5317[shot,:] = data_pico_rest[14]
    B6317[shot,:] = data_pico_rest[15]
    B3filt317[shot,:] = data_pico_rest[16]
    B4filt317[shot,:] = data_pico_rest[17]
    B5filt317[shot,:] = data_pico_rest[18]
    B6filt317[shot,:] = data_pico_rest[19]

z = raw.create_dataset('z', data=Bdot3raw317)
z = bdot.create_dataset('z', data=Bdot3317)
z = b.create_dataset('z', data=B3317)
z = Bfilt.create_dataset('z', data=B3filt317)

#I don't think the 317 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data317.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw317)
theta = raw.create_dataset('theta', data=Bdot5raw317)
z = raw.create_dataset('z', data=Bdot6raw317)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4317)
theta = bdot.create_dataset('theta', data=Bdot5317)
z = bdot.create_dataset('z', data=Bdot6317)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4317)
theta = b.create_dataset('theta', data=B5317)
z = b.create_dataset('z', data=B6317)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt317)
theta = bfilt.create_dataset('theta', data=B5filt317)
z = bfilt.create_dataset('z', data=B6filt317)

# Picoscope 3
for shot in np.arange(1, shotnum317):
    if shot==158:
        continue
    data_pico_rest = loadrest_317.load_picoscope_rest_03172022(shot, scopenum=3)
    
    Bdot3raw317[shot,:] = data_pico_rest[4].T
    Bdot4raw317[shot,:] = data_pico_rest[5].T
    Bdot5raw317[shot,:] = data_pico_rest[6].T
    Bdot6raw317[shot,:] = data_pico_rest[7].T
    Bdot3317[shot,:] = data_pico_rest[8].T
    Bdot4317[shot,:] = data_pico_rest[9].T
    Bdot5317[shot,:] = data_pico_rest[10].T
    Bdot6317[shot,:] = data_pico_rest[11].T
    B3317[shot,:] = data_pico_rest[12]
    B4317[shot,:] = data_pico_rest[13]
    B5317[shot,:] = data_pico_rest[14]
    B6317[shot,:] = data_pico_rest[15]
    B3filt317[shot,:] = data_pico_rest[16]
    B4filt317[shot,:] = data_pico_rest[17]
    B5filt317[shot,:] = data_pico_rest[18]
    B6filt317[shot,:] = data_pico_rest[19]
    
pos = data317.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw317)
theta = raw.create_dataset('theta', data=Bdot4raw317)
z = raw.create_dataset('z', data=Bdot5raw317)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3317)
theta = bdot.create_dataset('theta', data=Bdot4317)
z = bdot.create_dataset('z', data=Bdot5317)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3317)
theta = b.create_dataset('theta', data=B4317)
z = b.create_dataset('z', data=B5317)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt317)
theta = Bfilt.create_dataset('theta', data=B4filt317)
z = Bfilt.create_dataset('z', data=B5filt317)

#Pico3 has positions 19 and 21 in it.
pos = data317.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw317)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6317)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6317)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt317)

# Picoscope 4
for shot in np.arange(1,shotnum317):
    if shot==158:
        continue
    data_pico_rest = loadrest_317.load_picoscope_rest_03172022(shot, scopenum=4)
    
    Bdot3raw317[shot,:] = data_pico_rest[4].T
    Bdot4raw317[shot,:] = data_pico_rest[5].T
    Bdot5raw317[shot,:] = data_pico_rest[6].T
    Bdot6raw317[shot,:] = data_pico_rest[7].T
    Bdot3317[shot,:] = data_pico_rest[8].T
    Bdot4317[shot,:] = data_pico_rest[9].T
    Bdot5317[shot,:] = data_pico_rest[10].T
    Bdot6317[shot,:] = data_pico_rest[11].T
    B3317[shot,:] = data_pico_rest[12]
    B4317[shot,:] = data_pico_rest[13]
    B5317[shot,:] = data_pico_rest[14]
    B6317[shot,:] = data_pico_rest[15]
    B3filt317[shot,:] = data_pico_rest[16]
    B4filt317[shot,:] = data_pico_rest[17]
    B5filt317[shot,:] = data_pico_rest[18]
    B6filt317[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw317)
z = raw.create_dataset('z', data=Bdot4raw317)

theta = bdot.create_dataset('theta', data=Bdot3317)
z = bdot.create_dataset('z', data=Bdot4317)

theta = b.create_dataset('theta', data=B3317)
z = b.create_dataset('z', data=B4317)

theta = Bfilt.create_dataset('theta', data=B3filt317)
z = Bfilt.create_dataset('z', data=B4filt317)

# New position: 33
pos = data317.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw317)
theta = raw.create_dataset('theta', data=Bdot6raw317)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5317)
theta = bdot.create_dataset('theta', data=Bdot6317)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5317)
theta = b.create_dataset('theta', data=B6317)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt317)
theta = Bfilt.create_dataset('theta', data=B6filt317)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum317):
    if shot==158:
        continue
    data_pico_rest = loadrest_317.load_picoscope_rest_03172022(shot, scopenum=5)

    Bdot3raw317[shot,:] = data_pico_rest[4].T
    Bdot4raw317[shot,:] = data_pico_rest[5].T
    Bdot5raw317[shot,:] = data_pico_rest[6].T
    Bdot6raw317[shot,:] = data_pico_rest[7].T
    Bdot3317[shot,:] = data_pico_rest[8].T
    Bdot4317[shot,:] = data_pico_rest[9].T
    Bdot5317[shot,:] = data_pico_rest[10].T
    Bdot6317[shot,:] = data_pico_rest[11].T
    B3317[shot,:] = data_pico_rest[12]
    B4317[shot,:] = data_pico_rest[13]
    B5317[shot,:] = data_pico_rest[14]
    B6317[shot,:] = data_pico_rest[15]
    B3filt317[shot,:] = data_pico_rest[16]
    B4filt317[shot,:] = data_pico_rest[17]
    B5filt317[shot,:] = data_pico_rest[18]
    B6filt317[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw317)
z = bdot.create_dataset('z', data=Bdot3317)
z = b.create_dataset('z', data=B3317)
z = Bfilt.create_dataset('z', data=B3filt317)

#Last position, 35:
pos = data317.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw317)
theta = raw.create_dataset('theta', data=Bdot5raw317)
z = raw.create_dataset('z', data=Bdot6raw317)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4317)
theta = bdot.create_dataset('theta', data=Bdot5317)
z = bdot.create_dataset('z', data=Bdot6317)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4317)
theta = b.create_dataset('theta', data=B5317)
z = b.create_dataset('z', data=B6317)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt317)
theta = Bfilt.create_dataset('theta', data=B5filt317)
z = Bfilt.create_dataset('z', data=B6filt317)


Marked_Line_for_6202022 = np.zeros([1])
###########################
""" For dataset 6202022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info620 = data620.create_group('Probe Info')
probe_info620.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info620.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info620.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info620.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info620.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info620.create_dataset('radial_location', data = 'r=0cm')
probe_info620.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

disvolt_array = np.array([2.5, 3, 3.5, 4])
run_info620 = data620.create_group('Run Info')
run_info620.create_dataset('Discharge Voltage (kV)', data=disvolt_array)
run_info620.create_dataset('Stuffing Delay (ms)', data=1.5)
run_info620.create_dataset('Nozzle Fields', data='100mT Nozzle (39ms Coil Delay at 200V)')
run_info620.create_dataset('Shot Key', data='2.5kV Discharge Shots1-11, 3kV Discharge Shots11-30, 3.5kV Discharge Shots31-51, 4kV Discharge Shots52-54, (2.5kV Discharge) STEM Posse Shots55-56')

data620_pico1 = loadone_620.load_picoscope_one_06202022(1, scopenum=1)
data620_pico_rest = loadrest_620.load_picoscope_rest_06202022(1, scopenum=2)
shotnum620 = 57
time_s, time_us, timeB_s, timeB_us = data620_pico1[0],data620_pico1[1],data620_pico1[2],data620_pico1[3]

time620 = data620.create_group("Time")
time620.create_dataset('time_s', data=time_s)
time620.create_dataset('time_us', data = time_us)
time620.create_dataset('timeB_s', data = timeB_s)
time620.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw620 = np.zeros([shotnum620,datalength])
HighV_raw620 = np.zeros([shotnum620,datalength])
DisCurrent620 = np.zeros([shotnum620,datalength])
HighV620 = np.zeros([shotnum620,datalength])
Bdot1raw620 = np.zeros([shotnum620,datalength])
Bdot2raw620 = np.zeros([shotnum620,datalength])
Bdot3raw620 = np.zeros([shotnum620,datalength])
Bdot4raw620 = np.zeros([shotnum620,datalength])
Bdot5raw620 = np.zeros([shotnum620,datalength])
Bdot6raw620 = np.zeros([shotnum620,datalength])
Bdot1620 = np.zeros([shotnum620,datalength])
Bdot2620 = np.zeros([shotnum620,datalength])
Bdot3620 = np.zeros([shotnum620,datalength])
Bdot4620 = np.zeros([shotnum620,datalength])
Bdot5620 = np.zeros([shotnum620,datalength])
Bdot6620 = np.zeros([shotnum620,datalength])
B1620 = np.zeros([shotnum620,datalength-1])
B2620 = np.zeros([shotnum620,datalength-1])
B3620 = np.zeros([shotnum620,datalength-1])
B4620 = np.zeros([shotnum620,datalength-1])
B5620 = np.zeros([shotnum620,datalength-1])
B6620 = np.zeros([shotnum620,datalength-1])
B1filt620 = np.zeros([shotnum620,datalength-1])
B2filt620 = np.zeros([shotnum620,datalength-1])
B3filt620 = np.zeros([shotnum620,datalength-1])
B4filt620 = np.zeros([shotnum620,datalength-1])
B5filt620 = np.zeros([shotnum620,datalength-1])
B6filt620 = np.zeros([shotnum620,datalength-1])

#Time-Averaged Velocity Data
tavg_velocities_p57_6202022 = np.loadtxt(timeavg_data_directory_location+'Full6202022_ShotsandVelocities_P5P7.txt', skiprows=1, unpack=True).T
tavg_velocities_p19p21_6202022 = np.loadtxt(timeavg_data_directory_location+'Full6202022_ShotsandVelocities_P19P21.txt', skiprows=1, unpack=True).T
vel620 = data620.create_group("Time-Averaged Velocities")
vel620.create_dataset('Data Key', data=r'Pos5Pos7[0] = Shot#, P5P7[1] = Delay (${\mu}$s), [2] = Velocity (km/s), [3] = Nozzle Strength (mT), [4] = Discharge Voltage (kV)')
vel620.create_dataset('Pos5Pos7', data=tavg_velocities_p57_6202022)
vel620.create_dataset('Pos19Pos21', data=tavg_velocities_p19p21_6202022)

# Picoscope 1
for shot in np.arange (1, shotnum620):
    if shot==52:
        continue
    data620_pico1 = loadone_620.load_picoscope_one_06202022(shot, scopenum=1)
    
    DisCurrent_raw620[shot,:] = data_pico1[4].T
    HighV_raw620[shot,:] = data_pico1[5].T
    DisCurrent620[shot,:] = data_pico1[6]
    HighV620[shot,:] = data_pico1[7]
    Bdot1raw620[shot,:] = data_pico1[8].T
    Bdot2raw620[shot,:] = data_pico1[9].T
    Bdot1620[shot,:] = data_pico1[10].T
    Bdot2620[shot,:] = data_pico1[11].T
    B1620[shot,:] = data_pico1[12]
    B2620[shot,:] = data_pico1[13]
    B1filt620[shot,:] = data_pico1[14]
    B2filt620[shot,:] = data_pico1[15]

dis620 = data620.create_group("Discharge Data")
disraw620 = dis620.create_group("Discharge Raw")
DisCurrent_raw = disraw620.create_dataset("DisCurrent_raw", data=DisCurrent_raw620)
HighV_raw = disraw620.create_dataset("HighV_raw", data=HighV_raw620)
discharge620 = dis620.create_group("Discharge")
DisCurrent = discharge620.create_dataset("DisCurrent", data=DisCurrent620)
HighV = discharge620.create_dataset("HighV", data=HighV620)

pos = data620.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw620)
theta = raw.create_dataset('theta', data=Bdot2raw620)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1620)
theta = bdot.create_dataset('theta', data=Bdot2620)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1620)
theta = b.create_dataset('theta', data=B2620)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt620)
theta = Bfilt.create_dataset('theta', data=B2filt620)

# Picoscope 2
for shot in np.arange(1, shotnum620):
    if shot==45:
        continue
    data_pico_rest = loadrest_620.load_picoscope_rest_06202022(shot, scopenum=2)
    
    Bdot3raw620[shot,:] = data_pico_rest[4].T
    Bdot4raw620[shot,:] = data_pico_rest[5].T
    Bdot5raw620[shot,:] = data_pico_rest[6].T
    Bdot6raw620[shot,:] = data_pico_rest[7].T
    Bdot3620[shot,:] = data_pico_rest[8].T
    Bdot4620[shot,:] = data_pico_rest[9].T
    Bdot5620[shot,:] = data_pico_rest[10].T
    Bdot6620[shot,:] = data_pico_rest[11].T
    B3620[shot,:] = data_pico_rest[12]
    B4620[shot,:] = data_pico_rest[13]
    B5620[shot,:] = data_pico_rest[14]
    B6620[shot,:] = data_pico_rest[15]
    B3filt620[shot,:] = data_pico_rest[16]
    B4filt620[shot,:] = data_pico_rest[17]
    B5filt620[shot,:] = data_pico_rest[18]
    B6filt620[shot,:] = data_pico_rest[19]

z = raw.create_dataset('z', data=Bdot3raw620)
z = bdot.create_dataset('z', data=Bdot3620)
z = b.create_dataset('z', data=B3620)
z = Bfilt.create_dataset('z', data=B3filt620)

#I don't think the 620 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data620.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw620)
theta = raw.create_dataset('theta', data=Bdot5raw620)
z = raw.create_dataset('z', data=Bdot6raw620)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4620)
theta = bdot.create_dataset('theta', data=Bdot5620)
z = bdot.create_dataset('z', data=Bdot6620)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4620)
theta = b.create_dataset('theta', data=B5620)
z = b.create_dataset('z', data=B6620)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt620)
theta = bfilt.create_dataset('theta', data=B5filt620)
z = bfilt.create_dataset('z', data=B6filt620)

# Picoscope 3
for shot in np.arange(1, shotnum620):
    data_pico_rest = loadrest_620.load_picoscope_rest_06202022(shot, scopenum=3)

    Bdot3raw620[shot,:] = data_pico_rest[4].T
    Bdot4raw620[shot,:] = data_pico_rest[5].T
    Bdot5raw620[shot,:] = data_pico_rest[6].T
    Bdot6raw620[shot,:] = data_pico_rest[7].T
    Bdot3620[shot,:] = data_pico_rest[8].T
    Bdot4620[shot,:] = data_pico_rest[9].T
    Bdot5620[shot,:] = data_pico_rest[10].T
    Bdot6620[shot,:] = data_pico_rest[11].T
    B3620[shot,:] = data_pico_rest[12]
    B4620[shot,:] = data_pico_rest[13]
    B5620[shot,:] = data_pico_rest[14]
    B6620[shot,:] = data_pico_rest[15]
    B3filt620[shot,:] = data_pico_rest[16]
    B4filt620[shot,:] = data_pico_rest[17]
    B5filt620[shot,:] = data_pico_rest[18]
    B6filt620[shot,:] = data_pico_rest[19]

pos = data620.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw620)
theta = raw.create_dataset('theta', data=Bdot4raw620)
z = raw.create_dataset('z', data=Bdot5raw620)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3620)
theta = bdot.create_dataset('theta', data=Bdot4620)
z = bdot.create_dataset('z', data=Bdot5620)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3620)
theta = b.create_dataset('theta', data=B4620)
z = b.create_dataset('z', data=B5620)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt620)
theta = Bfilt.create_dataset('theta', data=B4filt620)
z = Bfilt.create_dataset('z', data=B5filt620)

#Pico3 has positions 19 and 21 in it.
pos = data620.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw620)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6620)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6620)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt620)

# Picoscope 4
for shot in np.arange(1,shotnum620):
    data_pico_rest = loadrest_620.load_picoscope_rest_06202022(shot, scopenum=4)
    
    Bdot3raw620[shot,:] = data_pico_rest[4].T
    Bdot4raw620[shot,:] = data_pico_rest[5].T
    Bdot5raw620[shot,:] = data_pico_rest[6].T
    Bdot6raw620[shot,:] = data_pico_rest[7].T
    Bdot3620[shot,:] = data_pico_rest[8].T
    Bdot4620[shot,:] = data_pico_rest[9].T
    Bdot5620[shot,:] = data_pico_rest[10].T
    Bdot6620[shot,:] = data_pico_rest[11].T
    B3620[shot,:] = data_pico_rest[12]
    B4620[shot,:] = data_pico_rest[13]
    B5620[shot,:] = data_pico_rest[14]
    B6620[shot,:] = data_pico_rest[15]
    B3filt620[shot,:] = data_pico_rest[16]
    B4filt620[shot,:] = data_pico_rest[17]
    B5filt620[shot,:] = data_pico_rest[18]
    B6filt620[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw620)
z = raw.create_dataset('z', data=Bdot4raw620)

theta = bdot.create_dataset('theta', data=Bdot3620)
z = bdot.create_dataset('z', data=Bdot4620)

theta = b.create_dataset('theta', data=B3620)
z = b.create_dataset('z', data=B4620)

theta = Bfilt.create_dataset('theta', data=B3filt620)
z = Bfilt.create_dataset('z', data=B4filt620)

# New position: 33
pos = data620.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw620)
theta = raw.create_dataset('theta', data=Bdot6raw620)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5620)
theta = bdot.create_dataset('theta', data=Bdot6620)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5620)
theta = b.create_dataset('theta', data=B6620)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt620)
theta = Bfilt.create_dataset('theta', data=B6filt620)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum620):
    data_pico_rest = loadrest_620.load_picoscope_rest_06202022(shot, scopenum=5)

    Bdot3raw620[shot,:] = data_pico_rest[4].T
    Bdot4raw620[shot,:] = data_pico_rest[5].T
    Bdot5raw620[shot,:] = data_pico_rest[6].T
    Bdot6raw620[shot,:] = data_pico_rest[7].T
    Bdot3620[shot,:] = data_pico_rest[8].T
    Bdot4620[shot,:] = data_pico_rest[9].T
    Bdot5620[shot,:] = data_pico_rest[10].T
    Bdot6620[shot,:] = data_pico_rest[11].T
    B3620[shot,:] = data_pico_rest[12]
    B4620[shot,:] = data_pico_rest[13]
    B5620[shot,:] = data_pico_rest[14]
    B6620[shot,:] = data_pico_rest[15]
    B3filt620[shot,:] = data_pico_rest[16]
    B4filt620[shot,:] = data_pico_rest[17]
    B5filt620[shot,:] = data_pico_rest[18]
    B6filt620[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw620)
z = bdot.create_dataset('z', data=Bdot3620)
z = b.create_dataset('z', data=B3620)
z = Bfilt.create_dataset('z', data=B3filt620)

#Last position, 35:
pos = data620.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw620)
theta = raw.create_dataset('theta', data=Bdot5raw620)
z = raw.create_dataset('z', data=Bdot6raw620)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4620)
theta = bdot.create_dataset('theta', data=Bdot5620)
z = bdot.create_dataset('z', data=Bdot6620)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4620)
theta = b.create_dataset('theta', data=B5620)
z = b.create_dataset('z', data=B6620)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt620)
theta = Bfilt.create_dataset('theta', data=B5filt620)
z = Bfilt.create_dataset('z', data=B6filt620)


Marked_Line_for_7152022 = np.zeros([1])
###########################
""" For dataset 7152022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info715 = data715.create_group('Probe Info')
probe_info715.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info715.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info715.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info715.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info715.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info715.create_dataset('radial_location', data = 'r=0cm')
probe_info715.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

run_info715 = data715.create_group('Run Info')
run_info715.create_dataset('Discharge Voltage (kV)', data=3.5)
run_info715.create_dataset('Stuffing Delay (ms)', data=1.5)
run_info715.create_dataset('Nozzle Fields', data='100mT Nozzle - All Shots')
run_info715.create_dataset('Shot Key', data='3.5kV Discharge, 100mT Nozzle Shots 1-51')

data715_pico1 = loadone_715.load_picoscope_one_07152022(1, scopenum=1)
data715_pico_rest = loadrest_715.load_picoscope_rest_07152022(1, scopenum=2)
shotnum715 = 52
time_s, time_us, timeB_s, timeB_us = data715_pico1[0],data715_pico1[1],data715_pico1[2],data715_pico1[3]

time715 = data715.create_group("Time")
time715.create_dataset('time_s', data=time_s)
time715.create_dataset('time_us', data = time_us)
time715.create_dataset('timeB_s', data = timeB_s)
time715.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw715 = np.zeros([shotnum715,datalength])
HighV_raw715 = np.zeros([shotnum715,datalength])
DisCurrent715 = np.zeros([shotnum715,datalength])
HighV715 = np.zeros([shotnum715,datalength])
Bdot1raw715 = np.zeros([shotnum715,datalength])
Bdot2raw715 = np.zeros([shotnum715,datalength])
Bdot3raw715 = np.zeros([shotnum715,datalength])
Bdot4raw715 = np.zeros([shotnum715,datalength])
Bdot5raw715 = np.zeros([shotnum715,datalength])
Bdot6raw715 = np.zeros([shotnum715,datalength])
Bdot1715 = np.zeros([shotnum715,datalength])
Bdot2715 = np.zeros([shotnum715,datalength])
Bdot3715 = np.zeros([shotnum715,datalength])
Bdot4715 = np.zeros([shotnum715,datalength])
Bdot5715 = np.zeros([shotnum715,datalength])
Bdot6715 = np.zeros([shotnum715,datalength])
B1715 = np.zeros([shotnum715,datalength-1])
B2715 = np.zeros([shotnum715,datalength-1])
B3715 = np.zeros([shotnum715,datalength-1])
B4715 = np.zeros([shotnum715,datalength-1])
B5715 = np.zeros([shotnum715,datalength-1])
B6715 = np.zeros([shotnum715,datalength-1])
B1filt715 = np.zeros([shotnum715,datalength-1])
B2filt715 = np.zeros([shotnum715,datalength-1])
B3filt715 = np.zeros([shotnum715,datalength-1])
B4filt715 = np.zeros([shotnum715,datalength-1])
B5filt715 = np.zeros([shotnum715,datalength-1])
B6filt715 = np.zeros([shotnum715,datalength-1])

#Time-Averaged Velocity Data
tavg_velocities_p57_7152022 = np.loadtxt(timeavg_data_directory_location+'Full7152022_ShotsandVelocities_P5P7.txt', skiprows=1, unpack=True).T
tavg_velocities_p19p21_7152022 = np.loadtxt(timeavg_data_directory_location+'Full7152022_ShotsandVelocities_P19P21.txt', skiprows=1, unpack=True).T
vel715 = data715.create_group("Time-Averaged Velocities")
vel715.create_dataset('Data Key', data=r'Pos5Pos7[0] = Shot#, P5P7[1] = Delay (${\mu}$s), [2] = Velocity (km/s), [3] = Nozzle Strength (mT), [4] = Discharge Voltage (kV)')
vel715.create_dataset('Pos5Pos7', data=tavg_velocities_p57_7152022)
vel715.create_dataset('Pos19Pos21', data=tavg_velocities_p19p21_7152022)

# Picoscope 1
for shot in np.arange (1, shotnum715):
    data715_pico1 = loadone_715.load_picoscope_one_07152022(shot, scopenum=1)
    
    DisCurrent_raw715[shot,:] = data_pico1[4].T
    HighV_raw715[shot,:] = data_pico1[5].T
    DisCurrent715[shot,:] = data_pico1[6]
    HighV715[shot,:] = data_pico1[7]
    Bdot1raw715[shot,:] = data_pico1[8].T
    Bdot2raw715[shot,:] = data_pico1[9].T
    Bdot1715[shot,:] = data_pico1[10].T
    Bdot2715[shot,:] = data_pico1[11].T
    B1715[shot,:] = data_pico1[12]
    B2715[shot,:] = data_pico1[13]
    B1filt715[shot,:] = data_pico1[14]
    B2filt715[shot,:] = data_pico1[15]

dis715 = data715.create_group("Discharge Data")
disraw715 = dis715.create_group("Discharge Raw")
DisCurrent_raw = disraw715.create_dataset("DisCurrent_raw", data=DisCurrent_raw715)
HighV_raw = disraw715.create_dataset("HighV_raw", data=HighV_raw715)
discharge715 = dis715.create_group("Discharge")
DisCurrent = discharge715.create_dataset("DisCurrent", data=DisCurrent715)
HighV = discharge715.create_dataset("HighV", data=HighV715)

pos = data715.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw715)
theta = raw.create_dataset('theta', data=Bdot2raw715)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1715)
theta = bdot.create_dataset('theta', data=Bdot2715)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1715)
theta = b.create_dataset('theta', data=B2715)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt715)
theta = Bfilt.create_dataset('theta', data=B2filt715)

# Picoscope 2
for shot in np.arange(1, shotnum715):
    #if shot==45:
    #    continue
    data_pico_rest = loadrest_715.load_picoscope_rest_07152022(shot, scopenum=2)
    
    Bdot3raw715[shot,:] = data_pico_rest[4].T
    Bdot4raw715[shot,:] = data_pico_rest[5].T
    Bdot5raw715[shot,:] = data_pico_rest[6].T
    Bdot6raw715[shot,:] = data_pico_rest[7].T
    Bdot3715[shot,:] = data_pico_rest[8].T
    Bdot4715[shot,:] = data_pico_rest[9].T
    Bdot5715[shot,:] = data_pico_rest[10].T
    Bdot6715[shot,:] = data_pico_rest[11].T
    B3715[shot,:] = data_pico_rest[12]
    B4715[shot,:] = data_pico_rest[13]
    B5715[shot,:] = data_pico_rest[14]
    B6715[shot,:] = data_pico_rest[15]
    B3filt715[shot,:] = data_pico_rest[16]
    B4filt715[shot,:] = data_pico_rest[17]
    B5filt715[shot,:] = data_pico_rest[18]
    B6filt715[shot,:] = data_pico_rest[19]

z = raw.create_dataset('z', data=Bdot3raw715)
z = bdot.create_dataset('z', data=Bdot3715)
z = b.create_dataset('z', data=B3715)
z = Bfilt.create_dataset('z', data=B3filt715)

#I don't think the 715 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data715.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw715)
theta = raw.create_dataset('theta', data=Bdot5raw715)
z = raw.create_dataset('z', data=Bdot6raw715)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4715)
theta = bdot.create_dataset('theta', data=Bdot5715)
z = bdot.create_dataset('z', data=Bdot6715)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4715)
theta = b.create_dataset('theta', data=B5715)
z = b.create_dataset('z', data=B6715)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt715)
theta = bfilt.create_dataset('theta', data=B5filt715)
z = bfilt.create_dataset('z', data=B6filt715)

# Picoscope 3
for shot in np.arange(1, shotnum715):
    data_pico_rest = loadrest_715.load_picoscope_rest_07152022(shot, scopenum=3)

    Bdot3raw715[shot,:] = data_pico_rest[4].T
    Bdot4raw715[shot,:] = data_pico_rest[5].T
    Bdot5raw715[shot,:] = data_pico_rest[6].T
    Bdot6raw715[shot,:] = data_pico_rest[7].T
    Bdot3715[shot,:] = data_pico_rest[8].T
    Bdot4715[shot,:] = data_pico_rest[9].T
    Bdot5715[shot,:] = data_pico_rest[10].T
    Bdot6715[shot,:] = data_pico_rest[11].T
    B3715[shot,:] = data_pico_rest[12]
    B4715[shot,:] = data_pico_rest[13]
    B5715[shot,:] = data_pico_rest[14]
    B6715[shot,:] = data_pico_rest[15]
    B3filt715[shot,:] = data_pico_rest[16]
    B4filt715[shot,:] = data_pico_rest[17]
    B5filt715[shot,:] = data_pico_rest[18]
    B6filt715[shot,:] = data_pico_rest[19]

pos = data715.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw715)
theta = raw.create_dataset('theta', data=Bdot4raw715)
z = raw.create_dataset('z', data=Bdot5raw715)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3715)
theta = bdot.create_dataset('theta', data=Bdot4715)
z = bdot.create_dataset('z', data=Bdot5715)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3715)
theta = b.create_dataset('theta', data=B4715)
z = b.create_dataset('z', data=B5715)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt715)
theta = Bfilt.create_dataset('theta', data=B4filt715)
z = Bfilt.create_dataset('z', data=B5filt715)

#Pico3 has positions 19 and 21 in it.
pos = data715.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw715)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6715)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6715)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt715)

# Picoscope 4
for shot in np.arange(1,shotnum715):
    data_pico_rest = loadrest_715.load_picoscope_rest_07152022(shot, scopenum=4)
    
    Bdot3raw715[shot,:] = data_pico_rest[4].T
    Bdot4raw715[shot,:] = data_pico_rest[5].T
    Bdot5raw715[shot,:] = data_pico_rest[6].T
    Bdot6raw715[shot,:] = data_pico_rest[7].T
    Bdot3715[shot,:] = data_pico_rest[8].T
    Bdot4715[shot,:] = data_pico_rest[9].T
    Bdot5715[shot,:] = data_pico_rest[10].T
    Bdot6715[shot,:] = data_pico_rest[11].T
    B3715[shot,:] = data_pico_rest[12]
    B4715[shot,:] = data_pico_rest[13]
    B5715[shot,:] = data_pico_rest[14]
    B6715[shot,:] = data_pico_rest[15]
    B3filt715[shot,:] = data_pico_rest[16]
    B4filt715[shot,:] = data_pico_rest[17]
    B5filt715[shot,:] = data_pico_rest[18]
    B6filt715[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw715)
z = raw.create_dataset('z', data=Bdot4raw715)

theta = bdot.create_dataset('theta', data=Bdot3715)
z = bdot.create_dataset('z', data=Bdot4715)

theta = b.create_dataset('theta', data=B3715)
z = b.create_dataset('z', data=B4715)

theta = Bfilt.create_dataset('theta', data=B3filt715)
z = Bfilt.create_dataset('z', data=B4filt715)

# New position: 33
pos = data715.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw715)
theta = raw.create_dataset('theta', data=Bdot6raw715)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5715)
theta = bdot.create_dataset('theta', data=Bdot6715)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5715)
theta = b.create_dataset('theta', data=B6715)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt715)
theta = Bfilt.create_dataset('theta', data=B6filt715)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum715):
    data_pico_rest = loadrest_715.load_picoscope_rest_07152022(shot, scopenum=5)

    Bdot3raw715[shot,:] = data_pico_rest[4].T
    Bdot4raw715[shot,:] = data_pico_rest[5].T
    Bdot5raw715[shot,:] = data_pico_rest[6].T
    Bdot6raw715[shot,:] = data_pico_rest[7].T
    Bdot3715[shot,:] = data_pico_rest[8].T
    Bdot4715[shot,:] = data_pico_rest[9].T
    Bdot5715[shot,:] = data_pico_rest[10].T
    Bdot6715[shot,:] = data_pico_rest[11].T
    B3715[shot,:] = data_pico_rest[12]
    B4715[shot,:] = data_pico_rest[13]
    B5715[shot,:] = data_pico_rest[14]
    B6715[shot,:] = data_pico_rest[15]
    B3filt715[shot,:] = data_pico_rest[16]
    B4filt715[shot,:] = data_pico_rest[17]
    B5filt715[shot,:] = data_pico_rest[18]
    B6filt715[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw715)
z = bdot.create_dataset('z', data=Bdot3715)
z = b.create_dataset('z', data=B3715)
z = Bfilt.create_dataset('z', data=B3filt715)

#Last position, 35:
pos = data715.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw715)
theta = raw.create_dataset('theta', data=Bdot5raw715)
z = raw.create_dataset('z', data=Bdot6raw715)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4715)
theta = bdot.create_dataset('theta', data=Bdot5715)
z = bdot.create_dataset('z', data=Bdot6715)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4715)
theta = b.create_dataset('theta', data=B5715)
z = b.create_dataset('z', data=B6715)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt715)
theta = Bfilt.create_dataset('theta', data=B5filt715)
z = Bfilt.create_dataset('z', data=B6filt715)


Marked_Line_for_7202022 = np.zeros([1])
###########################
""" For dataset 7202022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info720 = data720.create_group('Probe Info')
probe_info720.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info720.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info720.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info720.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info720.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info720.create_dataset('radial_location', data = 'r=0cm')
probe_info720.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos19/Pos21, Pos33/Pos35')

run_info720 = data720.create_group('Run Info')
run_info720.create_dataset('Discharge Voltage (kV)', data=3.5)
run_info720.create_dataset('Stuffing Delay (ms)', data=1.5)
run_info720.create_dataset('Nozzle Fields', data='0mT Nozzle shots 1-52, 100mT Nozzle Shots 53-180, ~50mT Nozle Shots 181-253')
run_info720.create_dataset('Shot Key', data='124mT Nozzle Shots1-61, 133mT Nozzle Shots62-95, ([Bad Data] Attempt at 5ms Nozzle (~0mT) shots96-203 - One coil was 39ms and the other at 5ms [Bad Data])')

data720_pico1 = loadone_720.load_picoscope_one_07202022(1, scopenum=1)
data720_pico_rest = loadrest_720.load_picoscope_rest_07202022(1, scopenum=2)
shotnum720 = 204
time_s, time_us, timeB_s, timeB_us = data720_pico1[0],data720_pico1[1],data720_pico1[2],data720_pico1[3]

time720 = data720.create_group("Time")
time720.create_dataset('time_s', data=time_s)
time720.create_dataset('time_us', data = time_us)
time720.create_dataset('timeB_s', data = timeB_s)
time720.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw720 = np.zeros([shotnum720,datalength])
HighV_raw720 = np.zeros([shotnum720,datalength])
DisCurrent720 = np.zeros([shotnum720,datalength])
HighV720 = np.zeros([shotnum720,datalength])
Bdot1raw720 = np.zeros([shotnum720,datalength])
Bdot2raw720 = np.zeros([shotnum720,datalength])
Bdot3raw720 = np.zeros([shotnum720,datalength])
Bdot4raw720 = np.zeros([shotnum720,datalength])
Bdot5raw720 = np.zeros([shotnum720,datalength])
Bdot6raw720 = np.zeros([shotnum720,datalength])
Bdot1720 = np.zeros([shotnum720,datalength])
Bdot2720 = np.zeros([shotnum720,datalength])
Bdot3720 = np.zeros([shotnum720,datalength])
Bdot4720 = np.zeros([shotnum720,datalength])
Bdot5720 = np.zeros([shotnum720,datalength])
Bdot6720 = np.zeros([shotnum720,datalength])
B1720 = np.zeros([shotnum720,datalength-1])
B2720 = np.zeros([shotnum720,datalength-1])
B3720 = np.zeros([shotnum720,datalength-1])
B4720 = np.zeros([shotnum720,datalength-1])
B5720 = np.zeros([shotnum720,datalength-1])
B6720 = np.zeros([shotnum720,datalength-1])
B1filt720 = np.zeros([shotnum720,datalength-1])
B2filt720 = np.zeros([shotnum720,datalength-1])
B3filt720 = np.zeros([shotnum720,datalength-1])
B4filt720 = np.zeros([shotnum720,datalength-1])
B5filt720 = np.zeros([shotnum720,datalength-1])
B6filt720 = np.zeros([shotnum720,datalength-1])

#Time-Averaged Velocity Data
tavg_velocities_p57_7202022 = np.loadtxt(timeavg_data_directory_location+'Full7202022_ShotsandVelocities_P5P7.txt', skiprows=1, unpack=True).T
tavg_velocities_p19p21_7202022 = np.loadtxt(timeavg_data_directory_location+'Full7202022_ShotsandVelocities_P19P21.txt', skiprows=1, unpack=True).T
vel720 = data720.create_group("Time-Averaged Velocities")
vel720.create_dataset('Data Key', data=r'Pos5Pos7[0] = Shot#, P5P7[1] = Delay (${\mu}$s), [2] = Velocity (km/s), [3] = Nozzle Strength (mT), [4] = Discharge Voltage (kV)')
vel720.create_dataset('Pos5Pos7', data=tavg_velocities_p57_7202022)
vel720.create_dataset('Pos19Pos21', data=tavg_velocities_p19p21_7202022)


# Picoscope 1
for shot in np.arange (1, shotnum720):
    data720_pico1 = loadone_720.load_picoscope_one_07202022(shot, scopenum=1)
    
    DisCurrent_raw720[shot,:] = data_pico1[4].T
    HighV_raw720[shot,:] = data_pico1[5].T
    DisCurrent720[shot,:] = data_pico1[6]
    HighV720[shot,:] = data_pico1[7]
    Bdot1raw720[shot,:] = data_pico1[8].T
    Bdot2raw720[shot,:] = data_pico1[9].T
    Bdot1720[shot,:] = data_pico1[10].T
    Bdot2720[shot,:] = data_pico1[11].T
    B1720[shot,:] = data_pico1[12]
    B2720[shot,:] = data_pico1[13]
    B1filt720[shot,:] = data_pico1[14]
    B2filt720[shot,:] = data_pico1[15]

dis720 = data720.create_group("Discharge Data")
disraw720 = dis720.create_group("Discharge Raw")
DisCurrent_raw = disraw720.create_dataset("DisCurrent_raw", data=DisCurrent_raw720)
HighV_raw = disraw720.create_dataset("HighV_raw", data=HighV_raw720)
discharge720 = dis720.create_group("Discharge")
DisCurrent = discharge720.create_dataset("DisCurrent", data=DisCurrent720)
HighV = discharge720.create_dataset("HighV", data=HighV720)

pos = data720.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw720)
theta = raw.create_dataset('theta', data=Bdot2raw720)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1720)
theta = bdot.create_dataset('theta', data=Bdot2720)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1720)
theta = b.create_dataset('theta', data=B2720)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt720)
theta = Bfilt.create_dataset('theta', data=B2filt720)

# Picoscope 2
for shot in np.arange(1, shotnum720):
    #if shot==45:
    #    continue
    data_pico_rest = loadrest_720.load_picoscope_rest_07202022(shot, scopenum=2)
    
    Bdot3raw720[shot,:] = data_pico_rest[4].T
    Bdot4raw720[shot,:] = data_pico_rest[5].T
    Bdot5raw720[shot,:] = data_pico_rest[6].T
    Bdot6raw720[shot,:] = data_pico_rest[7].T
    Bdot3720[shot,:] = data_pico_rest[8].T
    Bdot4720[shot,:] = data_pico_rest[9].T
    Bdot5720[shot,:] = data_pico_rest[10].T
    Bdot6720[shot,:] = data_pico_rest[11].T
    B3720[shot,:] = data_pico_rest[12]
    B4720[shot,:] = data_pico_rest[13]
    B5720[shot,:] = data_pico_rest[14]
    B6720[shot,:] = data_pico_rest[15]
    B3filt720[shot,:] = data_pico_rest[16]
    B4filt720[shot,:] = data_pico_rest[17]
    B5filt720[shot,:] = data_pico_rest[18]
    B6filt720[shot,:] = data_pico_rest[19]


z = raw.create_dataset('z', data=Bdot3raw720)
z = bdot.create_dataset('z', data=Bdot3720)
z = b.create_dataset('z', data=B3720)
z = Bfilt.create_dataset('z', data=B3filt720)

#I don't think the 720 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data720.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw720)
theta = raw.create_dataset('theta', data=Bdot5raw720)
z = raw.create_dataset('z', data=Bdot6raw720)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4720)
theta = bdot.create_dataset('theta', data=Bdot5720)
z = bdot.create_dataset('z', data=Bdot6720)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4720)
theta = b.create_dataset('theta', data=B5720)
z = b.create_dataset('z', data=B6720)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt720)
theta = bfilt.create_dataset('theta', data=B5filt720)
z = bfilt.create_dataset('z', data=B6filt720)

# Picoscope 3
for shot in np.arange(1, shotnum720):
    data_pico_rest = loadrest_720.load_picoscope_rest_07202022(shot, scopenum=3)

    Bdot3raw720[shot,:] = data_pico_rest[4].T
    Bdot4raw720[shot,:] = data_pico_rest[5].T
    Bdot5raw720[shot,:] = data_pico_rest[6].T
    Bdot6raw720[shot,:] = data_pico_rest[7].T
    Bdot3720[shot,:] = data_pico_rest[8].T
    Bdot4720[shot,:] = data_pico_rest[9].T
    Bdot5720[shot,:] = data_pico_rest[10].T
    Bdot6720[shot,:] = data_pico_rest[11].T
    B3720[shot,:] = data_pico_rest[12]
    B4720[shot,:] = data_pico_rest[13]
    B5720[shot,:] = data_pico_rest[14]
    B6720[shot,:] = data_pico_rest[15]
    B3filt720[shot,:] = data_pico_rest[16]
    B4filt720[shot,:] = data_pico_rest[17]
    B5filt720[shot,:] = data_pico_rest[18]
    B6filt720[shot,:] = data_pico_rest[19]


pos = data720.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw720)
theta = raw.create_dataset('theta', data=Bdot4raw720)
z = raw.create_dataset('z', data=Bdot5raw720)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3720)
theta = bdot.create_dataset('theta', data=Bdot4720)
z = bdot.create_dataset('z', data=Bdot5720)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3720)
theta = b.create_dataset('theta', data=B4720)
z = b.create_dataset('z', data=B5720)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt720)
theta = Bfilt.create_dataset('theta', data=B4filt720)
z = Bfilt.create_dataset('z', data=B5filt720)

#Pico3 has positions 19 and 21 in it.
pos = data720.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw720)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6720)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6720)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt720)

# Picoscope 4
for shot in np.arange(1,shotnum720):
    data_pico_rest = loadrest_720.load_picoscope_rest_07202022(shot, scopenum=4)
    
    Bdot3raw720[shot,:] = data_pico_rest[4].T
    Bdot4raw720[shot,:] = data_pico_rest[5].T
    Bdot5raw720[shot,:] = data_pico_rest[6].T
    Bdot6raw720[shot,:] = data_pico_rest[7].T
    Bdot3720[shot,:] = data_pico_rest[8].T
    Bdot4720[shot,:] = data_pico_rest[9].T
    Bdot5720[shot,:] = data_pico_rest[10].T
    Bdot6720[shot,:] = data_pico_rest[11].T
    B3720[shot,:] = data_pico_rest[12]
    B4720[shot,:] = data_pico_rest[13]
    B5720[shot,:] = data_pico_rest[14]
    B6720[shot,:] = data_pico_rest[15]
    B3filt720[shot,:] = data_pico_rest[16]
    B4filt720[shot,:] = data_pico_rest[17]
    B5filt720[shot,:] = data_pico_rest[18]
    B6filt720[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw720)
z = raw.create_dataset('z', data=Bdot4raw720)

theta = bdot.create_dataset('theta', data=Bdot3720)
z = bdot.create_dataset('z', data=Bdot4720)

theta = b.create_dataset('theta', data=B3720)
z = b.create_dataset('z', data=B4720)

theta = Bfilt.create_dataset('theta', data=B3filt720)
z = Bfilt.create_dataset('z', data=B4filt720)

# New position: 33
pos = data720.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw720)
theta = raw.create_dataset('theta', data=Bdot6raw720)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5720)
theta = bdot.create_dataset('theta', data=Bdot6720)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5720)
theta = b.create_dataset('theta', data=B6720)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt720)
theta = Bfilt.create_dataset('theta', data=B6filt720)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum720):
    data_pico_rest = loadrest_720.load_picoscope_rest_07202022(shot, scopenum=5)

    Bdot3raw720[shot,:] = data_pico_rest[4].T
    Bdot4raw720[shot,:] = data_pico_rest[5].T
    Bdot5raw720[shot,:] = data_pico_rest[6].T
    Bdot6raw720[shot,:] = data_pico_rest[7].T
    Bdot3720[shot,:] = data_pico_rest[8].T
    Bdot4720[shot,:] = data_pico_rest[9].T
    Bdot5720[shot,:] = data_pico_rest[10].T
    Bdot6720[shot,:] = data_pico_rest[11].T
    B3720[shot,:] = data_pico_rest[12]
    B4720[shot,:] = data_pico_rest[13]
    B5720[shot,:] = data_pico_rest[14]
    B6720[shot,:] = data_pico_rest[15]
    B3filt720[shot,:] = data_pico_rest[16]
    B4filt720[shot,:] = data_pico_rest[17]
    B5filt720[shot,:] = data_pico_rest[18]
    B6filt720[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw720)
z = bdot.create_dataset('z', data=Bdot3720)
z = b.create_dataset('z', data=B3720)
z = Bfilt.create_dataset('z', data=B3filt720)

#Last position, 35:
pos = data720.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw720)
theta = raw.create_dataset('theta', data=Bdot5raw720)
z = raw.create_dataset('z', data=Bdot6raw720)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4720)
theta = bdot.create_dataset('theta', data=Bdot5720)
z = bdot.create_dataset('z', data=Bdot6720)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4720)
theta = b.create_dataset('theta', data=B5720)
z = b.create_dataset('z', data=B6720)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt720)
theta = Bfilt.create_dataset('theta', data=B5filt720)
z = Bfilt.create_dataset('z', data=B6filt720)


Marked_Line_for_9232022 = np.zeros([1])
###########################
""" For dataset 9232022"""
###########################
probe_dia = 0.00158755  #meters, 1/16" probe
hole_sep = 0.001016 #meters (maybe check this, not sure what the difference is)
probe_info923 = data923.create_group('Probe Info')
probe_info923.create_dataset('Probe Stalk Diameter (m)', data=probe_dia)
probe_info923.create_dataset('Probe_Hole_Seperation(m)', data=hole_sep)
probe_info923.create_dataset(r'rloop_area ($m^2$)', data = np.pi*(probe_dia/2)**2)
probe_info923.create_dataset(r'tloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info923.create_dataset(r'zloop_area ($m^2$)', data = probe_dia*hole_sep)
probe_info923.create_dataset('radial_location', data = 'r=0cm')
probe_info923.create_dataset('Probe Locations', data = 'Pos5/Pos7, Pos6/Pos8, Davids New Probe')

disvolt_array = np.array([2, 3.5])
stuffing_array = np.array([1.5, 3])
gaspuff_array = np.array([2, 2.5, 5])
run_info923 = data923.create_group('Run Info')
run_info923.create_dataset('Discharge Voltage (kV)', data=disvolt_array)
run_info923.create_dataset('Stuffing Delay (ms)', data=stuffing_array)
run_info923.create_dataset('Gas Puff Delay (ms)', data=gaspuff_array)
run_info923.create_dataset('Shot Key', data='[2kV Discharge] 0mT Nozzle Shots37-57, 100mT Nozzle Shots58-88, 133mT Shot89 -Too Stuffed, 30ms Delay Shots90-99. [3.5kV Discharge] 0mT Nozzle Shots100-109 with misfires/delays, Shots110-113 [2kV] Phys121 Tour, Shots114-115 Delayed (NEW GASPUFF DELAY - 2.5ms), Shots116-120 Delayed (NEW GASPUFF DELAY - 5ms), Shot121 Delayed (NEW GASPUFF DELAY - 2ms) and 15mT Nozzle, Shots123-129 Misfires up to 45mT Nozzle, Shots 130-135 Misfires (Back to 0mT Nozzle), Shots136-138 (NEW STUFFING DELAY - 3ms (3mW Flux), Shots 139-141 Good, Shots 142-144 Misfire and PowerSupply Bang')


data923_pico1 = loadone_923.load_picoscope_one_09232022(1, scopenum=1)
data923_pico_rest = loadrest_923.load_picoscope_rest_09232022(1, scopenum=2)
shotnum923 = 145
time_s, time_us, timeB_s, timeB_us = data923_pico1[0],data923_pico1[1],data923_pico1[2],data923_pico1[3]

time923 = data923.create_group("Time")
time923.create_dataset('time_s', data=time_s)
time923.create_dataset('time_us', data = time_us)
time923.create_dataset('timeB_s', data = timeB_s)
time923.create_dataset('timeB_us', data = timeB_us)

#Discharge and Magnetic Data
DisCurrent_raw923 = np.zeros([shotnum923,datalength])
HighV_raw923 = np.zeros([shotnum923,datalength])
DisCurrent923 = np.zeros([shotnum923,datalength])
HighV923 = np.zeros([shotnum923,datalength])
Bdot1raw923 = np.zeros([shotnum923,datalength])
Bdot2raw923 = np.zeros([shotnum923,datalength])
Bdot3raw923 = np.zeros([shotnum923,datalength])
Bdot4raw923 = np.zeros([shotnum923,datalength])
Bdot5raw923 = np.zeros([shotnum923,datalength])
Bdot6raw923 = np.zeros([shotnum923,datalength])
Bdot1923 = np.zeros([shotnum923,datalength])
Bdot2923 = np.zeros([shotnum923,datalength])
Bdot3923 = np.zeros([shotnum923,datalength])
Bdot4923 = np.zeros([shotnum923,datalength])
Bdot5923 = np.zeros([shotnum923,datalength])
Bdot6923 = np.zeros([shotnum923,datalength])
B1923 = np.zeros([shotnum923,datalength-1])
B2923 = np.zeros([shotnum923,datalength-1])
B3923 = np.zeros([shotnum923,datalength-1])
B4923 = np.zeros([shotnum923,datalength-1])
B5923 = np.zeros([shotnum923,datalength-1])
B6923 = np.zeros([shotnum923,datalength-1])
B1filt923 = np.zeros([shotnum923,datalength-1])
B2filt923 = np.zeros([shotnum923,datalength-1])
B3filt923 = np.zeros([shotnum923,datalength-1])
B4filt923 = np.zeros([shotnum923,datalength-1])
B5filt923 = np.zeros([shotnum923,datalength-1])
B6filt923 = np.zeros([shotnum923,datalength-1])

#Time-Averaged Velocity Data
tavg_velocities_p57_9232022 = np.loadtxt(timeavg_data_directory_location + 'Full9232022_ShotsandVelocities_P5P7.txt', skiprows=1, unpack=True).T
tavg_velocities_p68_9232022 = np.loadtxt(timeavg_data_directory_location + 'Full9232022_ShotsandVelocities_P6P8.txt', skiprows=1, unpack=True).T
vel923 = data923.create_group("Time-Averaged Velocities")
vel923.create_dataset('Data Key', data=r'Pos5Pos7[0] = Shot#, P5P7[1] = Delay (${\mu}$s), [2] = Velocity (km/s), [3] = Nozzle Strength (mT), [4] = Discharge Voltage (kV)')
vel923.create_dataset('Pos5Pos7', data=tavg_velocities_p57_9232022)
vel923.create_dataset('Pos6Pos8', data=tavg_velocities_p68_9232022)

# Picoscope 1
for shot in np.arange (1, shotnum923):
    data923_pico1 = loadone_923.load_picoscope_one_09232022(shot, scopenum=1)
    
    DisCurrent_raw923[shot,:] = data_pico1[4].T
    HighV_raw923[shot,:] = data_pico1[5].T
    DisCurrent923[shot,:] = data_pico1[6]
    HighV923[shot,:] = data_pico1[7]
    Bdot1raw923[shot,:] = data_pico1[8].T
    Bdot2raw923[shot,:] = data_pico1[9].T
    Bdot1923[shot,:] = data_pico1[10].T
    Bdot2923[shot,:] = data_pico1[11].T
    B1923[shot,:] = data_pico1[12]
    B2923[shot,:] = data_pico1[13]
    B1filt923[shot,:] = data_pico1[14]
    B2filt923[shot,:] = data_pico1[15]

dis923 = data923.create_group("Discharge Data")
disraw923 = dis923.create_group("Discharge Raw")
DisCurrent_raw = disraw923.create_dataset("DisCurrent_raw", data=DisCurrent_raw923)
HighV_raw = disraw923.create_dataset("HighV_raw", data=HighV_raw923)
discharge923 = dis923.create_group("Discharge")
DisCurrent = discharge923.create_dataset("DisCurrent", data=DisCurrent923)
HighV = discharge923.create_dataset("HighV", data=HighV923)

pos = data923.create_group("pos5")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot1raw923)
theta = raw.create_dataset('theta', data=Bdot2raw923)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot1923)
theta = bdot.create_dataset('theta', data=Bdot2923)

b = pos.create_group('B')
r = b.create_dataset('r', data=B1923)
theta = b.create_dataset('theta', data=B2923)

Bfilt = pos.create_group('Bfilt')
r = Bfilt.create_dataset('r', data=B1filt923)
theta = Bfilt.create_dataset('theta', data=B2filt923)

# Picoscope 2
for shot in np.arange(1, shotnum923):
    #if shot==45:
    #    continue
    data_pico_rest = loadrest_923.load_picoscope_rest_09232022(shot, scopenum=2)
    
    Bdot3raw923[shot,:] = data_pico_rest[4].T
    Bdot4raw923[shot,:] = data_pico_rest[5].T
    Bdot5raw923[shot,:] = data_pico_rest[6].T
    Bdot6raw923[shot,:] = data_pico_rest[7].T
    Bdot3923[shot,:] = data_pico_rest[8].T
    Bdot4923[shot,:] = data_pico_rest[9].T
    Bdot5923[shot,:] = data_pico_rest[10].T
    Bdot6923[shot,:] = data_pico_rest[11].T
    B3923[shot,:] = data_pico_rest[12]
    B4923[shot,:] = data_pico_rest[13]
    B5923[shot,:] = data_pico_rest[14]
    B6923[shot,:] = data_pico_rest[15]
    B3filt923[shot,:] = data_pico_rest[16]
    B4filt923[shot,:] = data_pico_rest[17]
    B5filt923[shot,:] = data_pico_rest[18]
    B6filt923[shot,:] = data_pico_rest[19]


z = raw.create_dataset('z', data=Bdot3raw923)
z = bdot.create_dataset('z', data=Bdot3923)
z = b.create_dataset('z', data=B3923)
z = Bfilt.create_dataset('z', data=B3filt923)

#I don't think the 923 is necessary on the pos, since you're just coding the data into a location specified by the the right side of equals sign.
pos = data923.create_group('pos7')
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot4raw923)
theta = raw.create_dataset('theta', data=Bdot5raw923)
z = raw.create_dataset('z', data=Bdot6raw923)

bdot = pos.create_group('Bdot')
r = bdot.create_dataset('r', data=Bdot4923)
theta = bdot.create_dataset('theta', data=Bdot5923)
z = bdot.create_dataset('z', data=Bdot6923)

b = pos.create_group('B')
r = b.create_dataset('r', data=B4923)
theta = b.create_dataset('theta', data=B5923)
z = b.create_dataset('z', data=B6923)

bfilt = pos.create_group('Bfilt')
r = bfilt.create_dataset('r', data=B4filt923)
theta = bfilt.create_dataset('theta', data=B5filt923)
z = bfilt.create_dataset('z', data=B6filt923)

# Picoscope 3
for shot in np.arange(1, shotnum923):
    data_pico_rest = loadrest_923.load_picoscope_rest_09232022(shot, scopenum=3)

    Bdot3raw923[shot,:] = data_pico_rest[4].T
    Bdot4raw923[shot,:] = data_pico_rest[5].T
    Bdot5raw923[shot,:] = data_pico_rest[6].T
    Bdot6raw923[shot,:] = data_pico_rest[7].T
    Bdot3923[shot,:] = data_pico_rest[8].T
    Bdot4923[shot,:] = data_pico_rest[9].T
    Bdot5923[shot,:] = data_pico_rest[10].T
    Bdot6923[shot,:] = data_pico_rest[11].T
    B3923[shot,:] = data_pico_rest[12]
    B4923[shot,:] = data_pico_rest[13]
    B5923[shot,:] = data_pico_rest[14]
    B6923[shot,:] = data_pico_rest[15]
    B3filt923[shot,:] = data_pico_rest[16]
    B4filt923[shot,:] = data_pico_rest[17]
    B5filt923[shot,:] = data_pico_rest[18]
    B6filt923[shot,:] = data_pico_rest[19]


pos = data923.create_group("pos19")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot3raw923)
theta = raw.create_dataset('theta', data=Bdot4raw923)
z = raw.create_dataset('z', data=Bdot5raw923)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot3923)
theta = bdot.create_dataset('theta', data=Bdot4923)
z = bdot.create_dataset('z', data=Bdot5923)

b = pos.create_group("B")
r = b.create_dataset('r', data=B3923)
theta = b.create_dataset('theta', data=B4923)
z = b.create_dataset('z', data=B5923)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B3filt923)
theta = Bfilt.create_dataset('theta', data=B4filt923)
z = Bfilt.create_dataset('z', data=B5filt923)

#Pico3 has positions 19 and 21 in it.
pos = data923.create_group("pos21")
raw = pos.create_group('raw')
r = raw.create_dataset('r', data=Bdot6raw923)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot6923)

b = pos.create_group("B")
r = b.create_dataset('r', data=B6923)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B6filt923)

# Picoscope 4
for shot in np.arange(1,shotnum923):
    data_pico_rest = loadrest_923.load_picoscope_rest_09232022(shot, scopenum=4)
    
    Bdot3raw923[shot,:] = data_pico_rest[4].T
    Bdot4raw923[shot,:] = data_pico_rest[5].T
    Bdot5raw923[shot,:] = data_pico_rest[6].T
    Bdot6raw923[shot,:] = data_pico_rest[7].T
    Bdot3923[shot,:] = data_pico_rest[8].T
    Bdot4923[shot,:] = data_pico_rest[9].T
    Bdot5923[shot,:] = data_pico_rest[10].T
    Bdot6923[shot,:] = data_pico_rest[11].T
    B3923[shot,:] = data_pico_rest[12]
    B4923[shot,:] = data_pico_rest[13]
    B5923[shot,:] = data_pico_rest[14]
    B6923[shot,:] = data_pico_rest[15]
    B3filt923[shot,:] = data_pico_rest[16]
    B4filt923[shot,:] = data_pico_rest[17]
    B5filt923[shot,:] = data_pico_rest[18]
    B6filt923[shot,:] = data_pico_rest[19]

theta = raw.create_dataset('theta', data=Bdot3raw923)
z = raw.create_dataset('z', data=Bdot4raw923)

theta = bdot.create_dataset('theta', data=Bdot3923)
z = bdot.create_dataset('z', data=Bdot4923)

theta = b.create_dataset('theta', data=B3923)
z = b.create_dataset('z', data=B4923)

theta = Bfilt.create_dataset('theta', data=B3filt923)
z = Bfilt.create_dataset('z', data=B4filt923)

# New position: 33
pos = data923.create_group("pos33")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot5raw923)
theta = raw.create_dataset('theta', data=Bdot6raw923)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot5923)
theta = bdot.create_dataset('theta', data=Bdot6923)

b = pos.create_group("B")
r = b.create_dataset('r', data=B5923)
theta = b.create_dataset('theta', data=B6923)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B5filt923)
theta = Bfilt.create_dataset('theta', data=B6filt923)

# Picoscope 5, start with position 33
for shot in np.arange(1,shotnum923):
    data_pico_rest = loadrest_923.load_picoscope_rest_09232022(shot, scopenum=5)

    Bdot3raw923[shot,:] = data_pico_rest[4].T
    Bdot4raw923[shot,:] = data_pico_rest[5].T
    Bdot5raw923[shot,:] = data_pico_rest[6].T
    Bdot6raw923[shot,:] = data_pico_rest[7].T
    Bdot3923[shot,:] = data_pico_rest[8].T
    Bdot4923[shot,:] = data_pico_rest[9].T
    Bdot5923[shot,:] = data_pico_rest[10].T
    Bdot6923[shot,:] = data_pico_rest[11].T
    B3923[shot,:] = data_pico_rest[12]
    B4923[shot,:] = data_pico_rest[13]
    B5923[shot,:] = data_pico_rest[14]
    B6923[shot,:] = data_pico_rest[15]
    B3filt923[shot,:] = data_pico_rest[16]
    B4filt923[shot,:] = data_pico_rest[17]
    B5filt923[shot,:] = data_pico_rest[18]
    B6filt923[shot,:] = data_pico_rest[19]
    
z = raw.create_dataset('z', data=Bdot3raw923)
z = bdot.create_dataset('z', data=Bdot3923)
z = b.create_dataset('z', data=B3923)
z = Bfilt.create_dataset('z', data=B3filt923)

#Last position, 35:
pos = data923.create_group("pos35")
raw = pos.create_group("raw")
r = raw.create_dataset('r', data=Bdot4raw923)
theta = raw.create_dataset('theta', data=Bdot5raw923)
z = raw.create_dataset('z', data=Bdot6raw923)

bdot = pos.create_group("Bdot")
r = bdot.create_dataset('r', data=Bdot4923)
theta = bdot.create_dataset('theta', data=Bdot5923)
z = bdot.create_dataset('z', data=Bdot6923)

b = pos.create_group("B")
r = b.create_dataset('r', data=B4923)
theta = b.create_dataset('theta', data=B5923)
z = b.create_dataset('z', data=B6923)

Bfilt = pos.create_group("Bfilt")
r = Bfilt.create_dataset('r', data=B4filt923)
theta = Bfilt.create_dataset('theta', data=B5filt923)
z = Bfilt.create_dataset('z', data=B6filt923)
