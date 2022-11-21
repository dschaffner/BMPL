# load_data_basic.py
import h5py
import numpy as np
import matplotlib.pylab as plt


def printname(name):
    print(name)


def load_hdf5(file, verbose=False):
    f = h5py.File(file, 'r')
    if verbose:
        print('All Groups Contained')
        f.visit(printname)
    return f

def tindex_min(timearr, timevalue):
    minval = np.min(np.abs((timearr)-timevalue))
    tind = np.where(np.abs((timearr)-(timevalue)) == minval)
    tind = tind[0][0]
    return tind


directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
datafilename='Dataset_01122022.h5'
data=load_hdf5(directory+datafilename,verbose=True)

"""
datastructure = data['mag_probe']['positions'][probe#][direction][datatype][shot,index_range]

options: 
    probe# = ['probe5','probe7','probe19','probe21','probe33','probe35']
    direction = ['r','t','z']
    datatype = ['bdot','b']
    shot = [0 to 70]
Example: We want magnetic field data from probe5 in the r direction from shot 30 from index 3000 to 5000
example_array = data['mag_probe']['positions']['probe5']['r']['b'][30,3000:5000]
"""

time_s = data['time']['time_s']
timeB_s = time_s[1:]
time_us = data['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = tindex_min(analysis_start_time,timeB_us)
end_time_index = tindex_min(analysis_end_time,timeB_us)

data1=data['mag_probe']['positions']['probe5']['r']['b'][30,start_time_index:end_time_index]
plt.plot(timeB_us[start_time_index:end_time_index],data1)