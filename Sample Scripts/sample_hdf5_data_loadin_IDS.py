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
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2026\\07082026\\'

### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################


# place the following file in the directory indicated above
datafilename = 'IDS_testshots.h5'


# load hdf5 file
data = load_hdf5(data_directory_location+datafilename, verbose=True)

# the data storage array is organized by [']


# plot the 5th shot of z bdot for position 13 as a function of time in us
plt.figure(1)
time = data['time']['time_us'][:]
shot_chan1 = data['ids']['ids_array'][0,0,:]
shot_chan9 = data['ids']['ids_array'][0,8,:]
shot_chan12 = data['ids']['ids_array'][0,11,:]
plt.plot(time, shot_chan1,label='Shot 1, Channel 1')
plt.plot(time, shot_chan9,label='Shot 1, Channel 9')
plt.plot(time, shot_chan12,label='Shot 1, Channel 12')
plt.xlabel('Time [us]')
plt.ylabel('IDS Signal [V]')
plt.xlim(-15,400)
plt.ylim(-1.0,1.0)
plt.title('2kV Discharge')
plt.legend(loc='upper right',frameon=False)

plt.figure(2)
time = data['time']['time_us'][:]
shot_chan1 = data['ids']['ids_array'][2,0,:]
shot_chan9 = data['ids']['ids_array'][2,8,:]
shot_chan12 = data['ids']['ids_array'][2,11,:]
plt.plot(time, shot_chan1,label='Shot 3, Channel 1')
plt.plot(time, shot_chan9,label='Shot 3, Channel 9')
plt.plot(time, shot_chan12,label='Shot 3, Channel 12')
plt.xlabel('Time [us]')
plt.ylabel('IDS Signal [V]')
plt.xlim(-15,400)
plt.ylim(-1.0,1.0)
plt.title('3kV Discharge')
plt.legend(loc='upper right',frameon=False)

data.close()
