
import scipy.io as spio
import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import indexfinderfuncs as iff
import numpy as np
import spectrum_wwind as spec

#load dataset 1
#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06192024\\'
#datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
#data1=load_hdf5(directory+datafilename,verbose=True)

#directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06202024\\'
#datafilename = 'Dataset_06192024_2kV_1p5stuff_centerprobes.h5'#with disk
#data2=load_hdf5(directory+datafilename,verbose=True)

"""
# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\06262024\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_nocurrent_20shots.h5'
data1=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_2kV_1p5stuff_wire_4kA_20shots.h5'#without disk
data2=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_1kV_1stuff_wire_nocurrent_20shots.h5'#without disk
data3=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset_06262024_1kV_1stuff_wire_4kA_20shots.h5'#without disk
data4=load_hdf5(directory+datafilename,verbose=True)
"""

# ***** CHANGE Pathway!! *****************************************************
directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\Dataset 4 Summer 2024 Wire Target\\'
# ***** CHANGE Pathway!! *****************************************************

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_0kAwire_40shots.h5'
data1=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_2kAwire_40shots.h5'#without disk
data2=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_4kAwire_40shots.h5'#without disk
data3=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1kV_1msstuff_6kAwire_40shots.h5'#without disk
data4=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_0kAwire_40shots.h5'#without disk
data5=load_hdf5(directory+datafilename,verbose=True)

datafilename= 'Dataset4_wiretarget_1p5kV_1msstuff_8kAwire_40shots.h5'#without disk
data6=load_hdf5(directory+datafilename,verbose=True)

#magnetic probe time
time_s = data1['time']['time_s']
timeB_s = time_s[1:]
time_us = data1['time']['time_us']
timeB_us = time_us[1:]
analysis_start_time = 75
analysis_end_time = 175
start_time_index = iff.tindex_min(analysis_start_time,timeB_us)
end_time_index = iff.tindex_min(analysis_end_time,timeB_us)
timerange_limit = 3e-6#s
port_sep = 0.0254#m

#Bdot Triplets
numshots=40
directions=3
probes = 3

#determine FFT size and generate an output array
fsize=int((data1['efield_probe']['efield'][0,start_time_index:end_time_index].shape[0]/2)+1)

aveE_spec_1kV_0kA = np.zeros([fsize])
aveE_spec_1kV_2kA = np.zeros([fsize])
aveE_spec_1kV_4kA = np.zeros([fsize])
aveE_spec_1kV_6kA = np.zeros([fsize])
aveE_spec_1p5kV_0kA = np.zeros([fsize])
aveE_spec_1p5kV_8kA = np.zeros([fsize])
#spec_frombdot_wdisk = np.zeros([probes,directions,numshots,fsize])
#spec_frombdot_wodisk = np.zeros([probes,directions,numshots,fsize])

direction_list = ['r','t','z']

#loop over shots to read in data and compute FFT

for shot in np.arange(numshots):    
    data1ts=data1['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1kV_0kA=aveE_spec_1kV_0kA+pwr
    
    data1ts=data2['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1kV_2kA=aveE_spec_1kV_2kA+pwr
    
    data1ts=data3['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1kV_4kA=aveE_spec_1kV_4kA+pwr
    
    data1ts=data4['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1kV_6kA=aveE_spec_1kV_6kA+pwr
    
    data1ts=data5['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1p5kV_0kA=aveE_spec_1p5kV_0kA+pwr
    
    data1ts=data6['efield_probe']['efield'][shot,:]
    f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data1ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    aveE_spec_1p5kV_8kA=aveE_spec_1p5kV_8kA+pwr
    
    #data2ts=data2['efield_probe']['efield'][shot,:]
    #f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data2ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    #aveE_spec_2kV_4kA=aveE_spec_2kV_4kA+pwr    
    
    #data3ts=data3['efield_probe']['efield'][shot,:]
    #f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data3ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    #aveE_spec_1kV_0kA=aveE_spec_1kV_0kA+pwr
    
    #data4ts=data4['efield_probe']['efield'][shot,:]
    #f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data4ts[start_time_index:end_time_index],time_s[start_time_index:end_time_index],window='hamming')
    #aveE_spec_1kV_4kA=aveE_spec_1kV_4kA+pwr


#savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\Dataset 4 Summer 2024 Wire Target\\Plots\\Raw_Efield\\1p5kV_8kA\\'
#for shot in np.arange(numshots):    
#    data1ts=data6['efield_probe']['efield'][shot,:]
#    plt.plot(time_us[:],data1ts)
#    plt.title('Shot '+str(shot))
#    savefile = 'Shot'+str(shot)+'.png'
#    plt.savefig(savedirectory+savefile, dpi=300, facecolor='w', edgecolor='k')
#    plt.clf()