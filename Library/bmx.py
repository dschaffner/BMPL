# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:36:19 2019
 This will be a growing list of functions.
@author: CACS
"""
import numpy as np
import scipy.signal as sps


def BMX_Pico_Read(filename):
    """
        Simple function which outputs the data as read.
    """
    data = np.loadtxt(filename, skiprows=3, unpack=True)

    return data


def find_Index(data, data_point):
    ##############################################################################################
    """ This function returns the index corresponding to a given time. """
    ##############################################################################################
    output_Index = 0
    for i in range(0, data.size - 1):
        if (data[i] <= data_point):
            output_Index = i
    return output_Index
    
    return None

def finding_Index_Time(time_Data, time_Interest):
    ##############################################################################################
    """ This function returns the index corresponding to a given time. """
    ##############################################################################################
    time = time_Data
    output_Time = 0
    for i in range(0, time.size - 1):
        if (time[i]*1e6 <= time_Interest):
            output_Time = i
    return output_Time


#|| DAQ picoscope has resolution of 100MHz ||#
def high_Pass_Filter(data, filter_Freq, fs=100e6, N=5):
    #########################################################
    """
        Beginnings of a filter function.
        fs -- sampling frequency
    """
    #########################################################
    Wn = 2.0*filter_Freq/fs
    B, A = sps.butter(N, Wn, btype='highpass')
    output = sps.filtfilt(B, A, data)

    return output
