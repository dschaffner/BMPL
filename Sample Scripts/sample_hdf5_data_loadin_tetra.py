# sample_hdf5_data_loadin.py

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 23:31:43 2019

@author: dschaffner
"""

import matplotlib.pylab as plt
from load_hdf5 import load_hdf5

#######################################################################
# Directory style depends on Mac vs PC. For PC, use a double backslash.
# for a Mac, use a single forward slash.

### PC Style ###
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2025\03272025\\'

### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################


# place the following file in the directory indicated above
datafilename = 'sample_2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'


# load hdf5 file
data = load_hdf5(data_directory_location+datafilename, verbose=True)

# the data storage array is organized by [']


# plot the 5th shot of z bdot for position 13 as a function of time in us
time = data['time']['time_us'][:]
shot = data['pos13']['bdot']['z'][4, :]
plt.plot(time, shot)
plt.xlabel('Time [us]')
plt.ylabel('Bdot (uncalibrated)')
