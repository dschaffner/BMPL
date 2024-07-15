
"""
Created on Tue Aug  2 22:18:16 2022

@author: dschaffner
"""

import numpy as np
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff

def generate_deltas_list(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries)
    initial_index = 0
    final_index = int(indexstep)
    initial=timeseries[initial_index]
    final=timeseries[final_index]
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        delta = final-initial
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial=timeseries[initial_index]
        final=timeseries[final_index]
    return deltas


### Compute Flatness of Distribution Using a Gaussian Distribution Signal Example ###
### A Gaussian Distribution Should have a Flatness/Kurtosis = 3 ###
    
#generate gaussian distribution
mu, sigma = 0, 5 # mean and standard deviation
noise_arr = np.random.normal(mu, sigma, 1000000)

#actually data



#generate list of deltas based on time interval
dt=1 #unit timestep
tau=3 #interval (must be an integer)
noise_deltas = generate_deltas_list(noise_arr,dt,tau)

#convert deltas list to array
noise_deltas_arr = np.array(noise_deltas)

#center deltas
noise_mean = np.mean(noise_deltas_arr)
noise_deltas_arr = noise_deltas_arr-noise_mean

#compute min and max of deltas (for histogramming)
noise_deltas_min = np.floor(np.min(noise_deltas_arr))
noise_deltas_max = np.ceil(np.max(noise_deltas_arr))

#calculate histogram
nbins = 300
noise_hist, noise_bins = np.histogram(noise_deltas_arr,nbins,range=[noise_deltas_min,noise_deltas_max])



#shift bins to be centered
noise_bins_centered = 0.5*(noise_bins[1:]+noise_bins[:-1])

#plot histogram
plt.figure(1)
plt.plot(noise_bins_centered,noise_hist/max(noise_hist),label='noise')
plt.title('Normal Plot')
plt.xlim(-30,30)
plt.figure(2)
plt.semilogy(noise_bins_centered,noise_hist/max(noise_hist),label='noise')
plt.title('Semilog Y Plot')

#compute normalized S2 moment
noise_S2 = np.sum(noise_hist*noise_bins_centered**2)/np.sum(noise_hist)
noise_S4 = np.sum(noise_hist*noise_bins_centered**4)/np.sum(noise_hist)

#compute Flatness/Kurtosis
noise_flatness = noise_S4/(noise_S2)**2
print('Noise Flatness/Kurtosis = ',noise_flatness)

#######################################################################
#######################################################################
#######################################################################
#intermittency of actual data
#######################################################################
# Directory style depends on Mac vs PC. For PC, use a double backslash.
# for a Mac, use a single forward slash.
### PC Style ###
data_directory_location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\04232019\\processed\\'
### Mac Style ###
#data_director_location = '/user/dschaffner/data/'
#######################################################################
# place the following file in the directory indicated above
datafilename = 'sample_2kV_oddpos1to13_2ms_stuffdelay_17shots_04232019.h5'
# load hdf5 file
data = load_hdf5(data_directory_location+datafilename, verbose=True)
timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us/1e6

# select time range for FFT (in us)
start_time = 50
end_time = 150
time_range = [start_time, end_time]
# compute indices from time range
start_time_index = iff.tindex_min(start_time, timeB_us)
end_time_index = iff.tindex_min(end_time, timeB_us)
# select shots to analyze
first_shot = 1
last_shot = 16
numshots = (last_shot-first_shot)+1
shot_range = [first_shot, last_shot]

#create empty list for deltas
data_deltas = []

#loop over all shots to build list of deltas
for shot in np.arange(first_shot, last_shot+1):
    data_arr=data['pos1']['b']['z'][shot, start_time_index:end_time_index]
    data_deltas = data_deltas+generate_deltas_list(data_arr,dt,tau)
    
#covert data_deltas list into array
data_deltas_arr = np.array(data_deltas)

#center deltas
data_mean = np.mean(data_deltas_arr)
data_deltas_arr = data_deltas_arr-data_mean

#compute min and max of deltas (for histogramming)
data_deltas_min = np.floor(np.min(data_deltas_arr))
data_deltas_max = np.ceil(np.max(data_deltas_arr))

#calculate histogram
nbins = 300
data_hist, data_bins = np.histogram(data_deltas_arr,nbins,range=[data_deltas_min,data_deltas_max])

#shift bins to be centered
data_bins_centered = 0.5*(data_bins[1:]+data_bins[:-1])

#plot histogram
plt.figure(1)
plt.plot(data_bins_centered,data_hist/max(data_hist),label='data')
plt.legend()
plt.figure(2)
plt.semilogy(data_bins_centered,data_hist/max(data_hist),label='data')
plt.legend()

#compute normalized S2 moment
data_S2 = np.sum(data_hist*data_bins_centered**2)/np.sum(data_hist)
data_S4 = np.sum(data_hist*data_bins_centered**4)/np.sum(data_hist)

#compute Flatness/Kurtosis
data_flatness = data_S4/(data_S2)**2
print('Data Flatness/Kurtosis = ',data_flatness)
