#load picoscope

import numpy as np
import scipy.integrate as sp
from scipy.interpolate import interp1d
from scipy import signal
import matplotlib.pylab as plt


#data = '061615'
#shot = 1
#time_range = [20.0,80.0] #in us

def load_picoscope(shot_number,maxrange=5,scopenum=4,time_range=[-2.0,198.0],location='',plot=False):
    
    if scopenum == 1:
        scopename='pico1\\'
    if scopenum == 2:
        scopename='pico2\\'
    if scopenum == 3:
        scopename='pico3\\'
    if scopenum == 4:
        scopename='pico4\\'
    
    startintg_index=0#3000
    meancutoff = 1000
    #load file
    location = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2019\\Correlation Campaign\\Encoding Converted for PC\\10142019\\'
    filename = '20191014-0001 ('
    print(location+scopename+filename+str(shot_number)+').txt')
    data = np.loadtxt(location+scopename+filename+str(shot_number)+').txt',skiprows=2,unpack=True)

    #return data
    dataraw=data
    print(dataraw.shape)
    Discharge_raw=dataraw[2,:]
    HV_raw=dataraw[3,:]
    #data=data[:,startintg_index:]
    Discharge_raw[np.where(Discharge_raw==1)]=5.0
    Discharge_raw[np.where(Discharge_raw==-1)]=-5.0    
    
    HV_raw[np.where(HV_raw==1)]=5.0
    HV_raw[np.where(HV_raw==-1)]=-5.0
    
    time_ms = data[0,:]
    time_s = time_ms*1e-6
    
    Rogowski_gain = 10000.0
    Rogowski_dir = -1.0
    Rogowski_factor = 2.0
    HV_gain = 1000.0
    HV_dir = -1.0
    
    Discharge=Discharge_raw-np.mean(Discharge_raw[0:meancutoff])
    Discharge=Discharge*Rogowski_gain*Rogowski_dir*Rogowski_factor
    
    HV=HV_raw-np.mean(HV_raw[0:meancutoff])
    HV=HV*HV_gain*HV_dir
    
   
    return time_ms,time_s,Discharge_raw,HV_raw,Discharge,HV