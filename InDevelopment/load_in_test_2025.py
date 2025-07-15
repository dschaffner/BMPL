import numpy as np 
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5

data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2025\\06192025\\'

# place the following file in the directory indicated above
datafilename = 'TwoDens_twoBdot_2kV_p6msstuff_2mst1usgas_halfplateblocker_20shots.h5'


# load hdf5 file
data1 = load_hdf5(data_directory_location+datafilename, verbose=True)

isatdata=data1['isat probe']['isat']

time_us = data1['time']['time_us'][:]
timeB_us = data1['time']['time_us'][1:]

####figure 1
####probe 1 (near probe) 10th shot
#plt.plot(isatdata[0,9,:])
#plt.show()

####figure 2

####Aveage over shots

###create storage array
ave_isat1 = np.zeros([25004])
ave_isat2 = np.zeros([25004])

###loop over 20 shots
for loop in np.arange(20):
    print('Loop ',loop)
    ###select data from shot for probe 1
    isatshot = isatdata[0,loop,:]
    ###put data into stoarge array fir probe 1
    ave_isat1 = ave_isat1 + isatshot
    ###select data from shot for probe 2
    isatshot = isatdata[1,loop,:]
    ###put data into storage array for probe 2
    ave_isat2 = ave_isat2 + isatshot

ave_isat1 = ave_isat1/20.0
ave_isat2 = ave_isat2/20.0




#plt.plot(time_us,ave_isat1,label='Near Probe')
#plt.plot(time_us,ave_isat2,label='Far Probe')
#plt.xlabel('Time [us]')
#plt.ylabel('Isat [A]')
#plt.legend()
#plt.title('Average Isat over 20 shots')
#plt.show()





##practice magnetic data

ave_b = np.zeros([25003])
#Bdot Probe 1 The first shot
br = data1['mag_probe']['r']['b'][0,0,:]
bt = data1['mag_probe']['t']['b'][0,0,:]
bz = data1['mag_probe']['z']['b'][1,0,:]
#define Bmod (vector sum of all three component vectors)
Bmod = np.sqrt(br**2 + bt**2 + bz**2)

print(timeB_us.shape)
print(bz.shape)
#plt.plot(br)
#plt.plot(bt)
#plt.plot(bz)
#plt.plot(Bmod)
plt.plot(timeB_us,bz)







































