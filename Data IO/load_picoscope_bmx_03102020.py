#load picoscope

import numpy as np
import scipy.integrate as sp
from scipy.interpolate import interp1d
from scipy import signal
import matplotlib.pylab as plt
import sys


#data = '061615'
#shot = 1
#time_range = [20.0,80.0] #in us


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

def load_picoscope(shot_number, maxrange = 1, scopenum = 4, time_range = [-2.0, 198.0], location = '', plot = False):
    
    def butter_highpass(cutoff, fs, order = 5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff/nyq
        b, a = signal.butter(order, normal_cutoff, btype = 'highpass', analog = False)
        return b, a

    def butter_highpass_filter(data, cutoff, fs, order = 5):
        b, a = butter_highpass(cutoff, fs, order = order)
        y = signal.filtfilt(b, a, data)
        return y
    
    
    if (type(scopenum) == int):
        if scopenum == 1:
            scopename = '03102020pico1/'
        elif scopenum == 2:
            scopename = '03102020pico2/'
        elif scopenum == 3:
            scopename = '03102020pico3/'
        elif scopenum == 4:
            scopename = '03102020pico4/'
        elif scopenum == 5:
            scopename = '03102020pico5/'
        elif scopenum == 6:
            scopename = '03102020pico6/'
        elif scopenum == 7:
            scopename = '03102020pico7/'
        else:
            scopename = '03102020pico8/'
    else:
        print(f'scopenum is not an int, {scopenum}')
        sys.exit()
    
    probe_dia = 0.003175    #m (1/8'' probe)
    probe_dia = 0.00158755  #m (1/16'' probe)
    ##hole_sep = 0.001016     #m (1/16''probe)  ## Aparently unused variable
    r_probe_area = np.pi*(probe_dia/2)**2
    #tz_probe_area = probe_dia*hole_sep  ## Aparently unused variable
    startintg_index = 0 #3000
    meancutoff = 1000
    ##### load file
    # The location and filename lines must be updated to your system.
    location = '/Volumes/CarFlor/Research/Data/2020/03102020/'
    filename = '20200310-0001 ('
    
    print(location + scopename + filename + str(shot_number) + ').txt')
    try:
        data = np.loadtxt(location + scopename + filename + str(shot_number) + ').txt', skiprows = 2, unpack = True)
    except NameError as err:
        print("Double check you have updated the location variable to your OS system; mac, pc: ", err)
    ##### return data
    dataraw = data
    
    print(dataraw.shape)
    Bdotraw1 = dataraw[1, :]
    Bdotraw2 = dataraw[2, :]
    Bdotraw3 = dataraw[3, :]
    #Bdotraw4 = dataraw[4, :]
    data = data[:, startintg_index:]
    
    time_ms = data[0, :]
    time_s = time_ms*1e-6
    timeB_s = time_s[1:]
    timeB_ms = time_ms[1:]
    timeraw = dataraw[0, :]
        
    Bdot1 = data[1,:] - np.mean(data[1, 0:meancutoff])
    neginfs = np.isneginf(Bdot1)
    Bdot1[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot1)
    Bdot1[np.where(posinfs)] = maxrange
    
    Bdot2 = data[2,:] - np.mean(data[2, 0:meancutoff])
    neginfs = np.isneginf(Bdot2)
    Bdot2[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot2)
    Bdot2[np.where(posinfs)] = maxrange
    
    Bdot3 = data[3,:] - np.mean(data[3, 0:meancutoff])
    neginfs = np.isneginf(Bdot3)
    Bdot3[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot3)
    Bdot3[np.where(posinfs)] = maxrange
    
    #### 03102020 does not use the fourth pico port for magnetic data
    """Bdot4 = data[4,:] - np.mean(data[4, 0:meancutoff])
    neginfs = np.isneginf(Bdot4)
    Bdot4[np.where(neginfs)] = -maxrange
    posinfs = np.isposinf(Bdot4)
    Bdot4[np.where(posinfs)] = maxrange"""
    
    B1 = sp.cumtrapz(Bdot1/r_probe_area,time_s)*1e4 #Gauss
    B2 = sp.cumtrapz(Bdot2/r_probe_area,time_s)*1e4 #Gauss
    B3 = sp.cumtrapz(Bdot3/r_probe_area,time_s)*1e4 #Gauss
    #B4 = sp.cumtrapz(Bdot4/r_probe_area,time_s)*1e4 #Gauss
    #Bt7 = 3.162*sp.cumtrapz(Btdot7/tz_probe_area,time_s)*1e4#Gauss
    #Bt9 = 3.162*sp.cumtrapz(Btdot9/tz_probe_area,time_s)*1e4#Gauss
    #Bz7 = sp.cumtrapz(Bzdot7/tz_probe_area,time_s)*1e4#Gauss
    #Bz9 = sp.cumtrapz(Bzdot9/tz_probe_area,time_s)*1e4#Gauss
    #filtering

    #fps = 30
    #sine_fq = 10 #Hz
    #duration = 10 #seconds
    #sine_5Hz = sine_generator(fps,sine_fq,duration)
    #sine_fq = 1 #Hz
    #duration = 10 #seconds
    #sine_1Hz = sine_generator(fps,sine_fq,duration)

    #sine = sine_5Hz + sine_1Hz

    #filtered_sine = butter_highpass_filter(sine.data,10,fps)
          
    
    #Integration and Calibration    
    #Bx =sp.cumtrapz(Bxdot/probe_area,time_s)
    #Bx = 3.162*Bx/1.192485591065652224e-03
        
    #By =sp.cumtrapz(Bydot/probe_area,time_s)
    #By = 3.162*By/1.784763055992550198e-03
        
    #Bz =sp.cumtrapz(Bzdot/probe_area,time_s)
    #Bz = 3.162*Bz/1.297485014039849059e-03
    #meanBx = np.mean(Bx)
    # Filtering
    B1filt = butter_highpass_filter(B1, 5e4, 125e6, order = 3)  
    B2filt = butter_highpass_filter(B2, 5e4, 125e6, order = 3)  
    B3filt = butter_highpass_filter(B3, 5e4, 125e6, order = 3) 
    #B4filt = butter_highpass_filter(B4, 5e4, 125e6, order = 3)  
    #Btot = np.sqrt(Bxfilt**2+Byfilt**2+Bzfilt**2)
    #Btotave=Btotave+Btot
    
    #if plot:
    #    plt.figure(1)
    #    plt.plot(time,data[1,:])
    #    plt.figure(2)
    #    plt.plot(time[1:],Btot)
    
    return time_ms, time_s, timeB_s, timeB_ms, Bdot1, Bdot2, Bdot3, B1, B2, B3, B1filt, B2filt, B3filt, Bdotraw1, Bdotraw2, Bdotraw3, timeraw
    
    #return time_ms, time_s, timeB_s, timeB_ms, Bdot1, Bdot2, Bdot3, Bdot4, B1, B2, B3, B4, B1filt, B2filt, B3filt, B4filt, Bdotraw1, Bdotraw2, Bdotraw3, Bdotraw4, timeraw
