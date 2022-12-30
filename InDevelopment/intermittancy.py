
"""
Created on Tue Aug  2 22:18:16 2022

@author: dschaffner
"""

import numpy as np

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
mu, sigma = 0, 0.1 # mean and standard deviation
noise_arr = np.random.normal(mu, sigma, 100000)

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

#compute normalized S2 moment
noise_S2 = np.sum(noise_hist*noise_bins_centered**2)/np.sum(noise_hist)
noise_S4 = np.sum(noise_hist*noise_bins_centered**4)/np.sum(noise_hist)

#compute Flatness/Kurtosis
noise_flatness = noise_S4/(noise_S2)**2
print('Flatness = ',noise_flatness)

