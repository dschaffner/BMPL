# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:36:19 2019
 This will be a growing list of functions.
@author: Shadow
"""
#Github is dumb
import numpy as np
import matplotlib.pylab as plt
import scipy.integrate as sp
import scipy.signal as sps

def BMX_Pico_Read(filename):
    """
        Simple function which outputs the data as read.
    """
    data = np.loadtxt(filename, skiprows = 3, unpack = True)
    
    return data

def BMX_Single_Plot(data_Array, time_Array):
    return None

def BMX_Clear_Infs(input_Data, max_Range):
    """
        This functions changes the Inf values which plague the data.
    """
    neg_Infs = np.isneginf(input_Data)
    input_Data[np.where(neg_Infs)] = -max_Range
    pos_Infs = np.isposinf(input_Data)
    input_Data[np.where(pos_Infs)] = max_Range
    
    return input_Data

def BMX_Magnitude(V1,V2,V3):
    
    Vtot = np.sqrt(V1**2 + V2**2 + V3**2)
    
    return Vtot

def BMX_Correlation():
    return None

def bDot_Array_Name(indicator):
    """try:
        isinstance(indicator, str)
    except NameError:
        print("Provide a string for data_Structure variable.")"""
    
    if indicator.lower() == 'z':
        output = 'Bz_Dot'
    elif indicator.lower() == 'r':
        output = 'Br_Dot'
    elif indicator.lower() == 't':
        output = 'Bt_Dot'
    elif indicator.lower() == 'c':
        output = 'Current'
    elif indicator.lower() == 'v':
        output = 'Voltage'
    elif indicator.lower() == 'h':
        output = 'H-Alpha'
    
    return output

def BMX_Magnetic_BDOT(filename, data_Structure, output_Type = 'bdot', starting_Time = 0, mean_Subtraction = True, apply_Filter = True, cutoff_Region = False, max_Range = 1, ending_Time = -1, filter_Freq = 1e7):
    ####################################################################################################################
    """ This will be a general loading script. For now it will only contain barebones. """
    ####################################################################################################################
    def data_Struct():
        
        try:
            isinstance(data_Structure, str)
        except NameError:
            print("Provide a string for data_Structure variable.")
        
        ####################################################################################################################
        """ The two-component probes. """
        ####################################################################################################################
        if data_Structure.lower() == 'ztzt':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[0]]
            data_List2 = [data[3], data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'tztz':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[1], data[0]]
            data_List2 = [data[4], data[3], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), bDot_Array_Name(data_Structure[2]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'ztcv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[0]]
            data_List2 = [data[3], data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'zrzr':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[0]]
            data_List2 = [data[3], data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'rzrz':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[1], data[0]]
            data_List2 = [data[4], data[3], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), bDot_Array_Name(data_Structure[2]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'zrcv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[0]]
            data_List2 = [data[3], data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        else:
            pass
        ####################################################################################################################
        """ The three-component probes """
        ####################################################################################################################
            
        if data_Structure.lower() == 'zrtc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'ztrc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[3], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))

        elif data_Structure.lower() == 'rztc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[1], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'trzc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[2], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'tzrc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[3], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'rtzc':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[1], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'zrth':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'ztrh':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[3], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'rzth':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[1], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'trzh':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[2], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'tzrh':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[3], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'rtzh':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[1], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'zrtv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[2], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'ztrv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[1], data[3], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'rztv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[2], data[1], data[3], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[2]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
            
        elif data_Structure.lower() == 'trzv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[2], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'tzrv':
            """ I want data in form: z, r, t, time """          
            data_List1 = [data[2], data[3], data[1], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[1]), bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        elif data_Structure.lower() == 'rtzv':
            """ I want data in form: z, r, t, time """
            
            data_List1 = [data[3], data[1], data[2], data[0]]
            data_List2 = [data[4], data[0]]
            name_List1 = [bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]), 'Time']
            name_List2 = [bDot_Array_Name(data_Structure[3]), 'Time']
            
            my_Dict1 = dict(zip(name_List1, data_List1))
            my_Dict2 = dict(zip(name_List2, data_List2))
        
        else:
            pass
        
        two_Component = "ztzt, tztz, zrzr, rzrz, ztcv, zrcv"
        three_Component = ("zrtc, rztc, ztrc, rtzc, tzrc, trzc, "
                            + "zrth, rzth, ztrh, rtzh, tzrh, trzh, " 
                            + "zrtv, rztv, ztrv, rtzv, tzrv, trzv")
        try:
            my_Dict1
        except NameError:
            print("The provided data structure, %s,  is not included in the list: Two-Component %s; Three Component %s") %(data_Structure, two_Component, three_Component)
            print("Contact Carlos if you believe it needs to be added to the list.")
        return my_Dict1, my_Dict2
    
    data = BMX_Pico_Read(filename)
    my_Data1, my_Data2 = data_Struct()
    """
    my_Dict_1 = {'Time_B': timeB_Sec, 'Bz': Bz, 'Bt': Bt}
    """
    if output_Type.lower() == 'bdot':
        if len(my_Data1) == 3:
            name_List = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]),
                          bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3])]
            data1 = my_Data1['%s'%name_List[0]]
            data2 = my_Data1['%s'%name_List[1]]
            time = my_Data1['Time']
            data3 = my_Data2['%s'%name_List[2]]
            data4 = my_Data2['%s'%name_List[3]]
            seperator = ', '
            print("The output is of the form: " + seperator.join(name_List[0:2]) + ', Time, ' + seperator.join(name_List[2:4]))
            
        else:
            name_List = [bDot_Array_Name(data_Structure[0]), bDot_Array_Name(data_Structure[1]),
                          bDot_Array_Name(data_Structure[2]), bDot_Array_Name(data_Structure[3])]
            data1 = my_Data1['%s'%name_List[0]]
            data2 = my_Data1['%s'%name_List[1]]
            data3 = my_Data2['%s'%name_List[2]]
            time = my_Data1['Time']
            data4 = my_Data2['%s'%name_List[3]]
            seperator = ', '
            print("The output is of the form: " + seperator.join(name_List[0:2]) + ', Time, ' + seperator.join(name_List[2:4]))
        if mean_Subtraction:
            ending_Index = finding_Index_Time(time, 0)
            
            data1 -= np.mean(data1[:ending_Index])
            data2 -= np.mean(data2[:ending_Index])
            data3 -= np.mean(data3[:ending_Index])
            data4 -= np.mean(data4[:ending_Index])
        ####################################################################################################################
        """ This section of the code is focused on the time window. """
        ####################################################################################################################
        if cutoff_Region:
            ending_Index = finding_Index_Time(time, 0)
            
            time = time[:ending_Index]
            data1 = data1[:ending_Index]
            data2 = data2[:ending_Index]
            data3 = data3[:ending_Index]
            data4 = data4[:ending_Index]
        elif (starting_Time != 0 and ending_Time != -1):
            starting_Index = finding_Index_Time(time, starting_Time)
            ending_Index = finding_Index_Time(time, ending_Time)
            
            time = time[starting_Index:ending_Index]
            data1 = data1[starting_Index:ending_Index]
            data2 = data2[starting_Index:ending_Index]
            data3 = data3[starting_Index:ending_Index]
            data4 = data4[starting_Index:ending_Index]
        else:
            pass
        
        """elif (ending_Time == -1):
            starting_Index = finding_Index_Time(time, starting_Time)
            
            time = time[starting_Index:]
            data1 = data1[starting_Index:]
            data2 = data2[starting_Index:]
            data3 = data3[starting_Index:]
            data4 = data4[starting_Index:]
        """
    if apply_Filter:
        
        fs = 125e6
        N = 4
        Wn = 2.0*filter_Freq/fs
        B, A = sps.butter(N, Wn, output = 'ba')
        
        data1 = sps.filtfilt(B,A, data1)
        data2 = sps.filtfilt(B,A, data2)
        data3 = sps.filtfilt(B,A, data3)
        data4 = sps.filtfilt(B,A, data4)
    else:
        pass    
    
    return data1, data2, time, data3, data4

def BMX_Magnetic_Field_ZT_PM(filename, starting_Index = 0, mean_Cutoff = 0 , max_Range = 1, ending_Index = -1,
                        ztprobe_Dia = 0.00158755, hole_Sep = 0.001016):
    ################################################################################################
    """ The reason I am creating yet another load function has to do with a new mode of operation 
        pair mode (PM). Each picoscope, except the last one contains pairs of data corresponding 
        to a specific port. The first case neglects the r-component data, the data with the least 
        activity. """
    ################################################################################################

    data = BMX_Pico_Read(filename)
    if ending_Index == 0:
        ending_Index = data[0,:].size
    else:
        pass

    ################################################################################################
    """ Explicitely mean substracting data """
    ################################################################################################
    for i in range(1,5):
        data[i] -= np.mean(data[i][0:mean_Cutoff])
     
    data = data[:,starting_Index:ending_Index]
    ################################################################################################
    """ Setting the area variables """
    ################################################################################################
    tzprobe_Area = ztprobe_Dia*hole_Sep
    ################################################################################################
    """ Setting the time variables """
    ################################################################################################
    time_Microsec = data[0,:]
    time_Sec = time_Microsec*1e-6
    timeB_Sec = time_Sec[1:]
    ################################################################################################
    """ My standard is {Z, R, Theta} """
    ################################################################################################
    Bt_1dot = BMX_Clear_Infs(data[1], max_Range)
    Bz_1dot = BMX_Clear_Infs(data[2], max_Range)
    Bt_2dot = BMX_Clear_Infs(data[3], max_Range)
    Bz_2dot = BMX_Clear_Infs(data[4], max_Range)
    ################################################################################################
    ''' Performing the integration using the trap-rule (Guass units)'''
    ################################################################################################
    Bz_1 = sp.cumtrapz(Bz_1dot/(tzprobe_Area),time_Sec)*1e4
    Bt_1 = sp.cumtrapz(Bt_1dot/(tzprobe_Area),time_Sec)*1e4
    Bz_2 = sp.cumtrapz(Bz_2dot/(tzprobe_Area),time_Sec)*1e4
    Bt_2 = sp.cumtrapz(Bt_2dot/(tzprobe_Area),time_Sec)*1e4
    ################################################################################################
    ''' Detrending the field components '''
    ################################################################################################
    Bz_1 = sps.detrend(Bz_1)
    Bt_1 = sps.detrend(Bt_1)
    Bz_2 = sps.detrend(Bz_2)
    Bt_2 = sps.detrend(Bt_2)
    ################################################################################################
    ''' Mean subtracting and normalizing to mean :::: This should be in the beginning. '''
    ################################################################################################
    Bz_1mean = np.mean(Bz_1[:mean_Cutoff])
    Bt_1mean = np.mean(Bt_1[:mean_Cutoff])
    Bz_2mean = np.mean(Bz_2[:mean_Cutoff])
    Bt_2mean = np.mean(Bt_2[:mean_Cutoff])

    Bz_1 -= Bz_1mean
    Bt_1 -= Bt_1mean
    Bz_2 -= Bz_2mean
    Bt_2 -= Bt_2mean
    
    print(" Returning: timeB_Sec, Bz_1, Bt_1, Bz_2, Bt_2 ")
    
    return  timeB_Sec, Bz_1, Bt_1, Bz_2, Bt_2
    
def BMX_Bdot_ZT_PM(filename, starting_Index = 0, mean_Cutoff = 0, max_Range = 1, ending_Index = -1):
    ################################################################################################
    """ The reason I am creating yet another load function has to do with a new mode of operation 
        pair mode (PM). Each picoscope, except the last one contains pairs of data corresponding 
        to a specific port. The first case neglects the r-component data, the data with the least 
        activity. """
    ################################################################################################

    data = BMX_Pico_Read(filename)
    if ending_Index == 0:
        ending_Index = data[0,:].size
    else:
        pass
    
    ################################################################################################
    """ Explicitly mean substracting data """
    ################################################################################################
    for i in range(1,5):
     data[i] -= np.mean(data[i][0:mean_Cutoff])
    
    data = data[:,starting_Index:ending_Index]

    ################################################################################################
    """ Setting the time variables"""
    ################################################################################################
    time_Microsec = data[0,:]
    time_Sec = time_Microsec*1e-6
    timeB_Sec = time_Sec[1:]
    ################################################################################################
    """ My standard is {Z, R, Theta} """
    ################################################################################################
    
    Bt_1dot = BMX_Clear_Infs(data[1], max_Range)
    Bz_1dot = BMX_Clear_Infs(data[2], max_Range)
    Bt_2dot = BMX_Clear_Infs(data[3], max_Range)
    Bz_2dot = BMX_Clear_Infs(data[4], max_Range)
    ################################################################################################
    ''' Performing the integration using the trap-rule (Guass units)'''
    ################################################################################################
    ################################################################################################
    ''' Mean subtracting '''
    ################################################################################################
    """Bz_1dotmean = np.mean(Bz_1dot[:mean_Cutoff])
    Bt_1dotmean = np.mean(Bt_1dot[:mean_Cutoff])
    Bz_2dotmean = np.mean(Bz_2dot[:mean_Cutoff])
    Bt_2dotmean = np.mean(Bt_2dot[:mean_Cutoff])

    Bz_1dot -= Bz_1dotmean
    Bt_1dot -= Bt_1dotmean
    Bz_2dot -= Bz_2dotmean
    Bt_2dot -= Bt_2dotmean"""
    
    print(" Returning: time, Bz_1dot, Bt_1dot, Bz_2dot, Bt_2dot ")
    
    return  time_Sec, Bz_1dot, Bt_1dot, Bz_2dot, Bt_2dot

def BMX_Magnetic_Fieldf(filename, starting_Index = 0, mean_Cutoff = 0, max_Range = 1, ending_Index = -1,
                       rprobe_Dia = 0.003175, ztprobe_Dia = 0.00158755,
                       hole_Sep = 0.001016):
    """
        Outputs the magnetic field strengths, time and monitor value.
    """
    data = BMX_Pico_Read(filename)
    if ending_Index == 0:
        ending_Index = data[0,:].size
    else:
        pass
    ################################################################################################
    """ Explicitely mean subtracting data """
    ################################################################################################
    for i in range(1,4):
        data[i] -= np.mean(data[i][0:mean_Cutoff])
    data = data[:,starting_Index:ending_Index]
    ################################################################################################
    """ Setting the area variables """
    ################################################################################################
    rprobe_Area = np.pi*(rprobe_Dia/2)**2
    tzprobe_Area = ztprobe_Dia*hole_Sep
    ################################################################################################
    """ Setting the time variables"""
    ################################################################################################
    time_Microsec = data[0,:]
    time_Sec = time_Microsec*1e-6
    timeB_Sec = time_Sec[1:]
    ################################################################################################
    """ My standard is {Z, R, Theta} """
    ################################################################################################
    Bzdot = BMX_Clear_Infs(data[1], max_Range)
    Brdot = BMX_Clear_Infs(data[2], max_Range)
    Btdot = BMX_Clear_Infs(data[3], max_Range)
    ################################################################################################
    ''' Performing the integration using the trap-rule (Guass units): Start Index 0.'''
    ################################################################################################
    Bz = sp.cumtrapz(Bzdot/(tzprobe_Area),time_Sec)*1e4
    Br = sp.cumtrapz(Brdot/(rprobe_Area),time_Sec)*1e4
    Bt = sp.cumtrapz(Btdot/(tzprobe_Area),time_Sec)*1e4
    ################################################################################################
    ''' Detrending the field components ::::  '''
    ################################################################################################
    Bz = sps.detrend(Bz)
    Br = sps.detrend(Br)
    Bt = sps.detrend(Bt)
    ################################################################################################
    ''' Mean subtracting '''
    ################################################################################################
    Bz_mean = np.mean(Bz[:mean_Cutoff])
    Br_mean = np.mean(Br[:mean_Cutoff])
    Bt_mean = np.mean(Bt[:mean_Cutoff])

    Bz -= Bz_mean
    Br -= Br_mean
    Bt -= Bt_mean

    return  timeB_Sec, Bz, Br, Bt

def BMX_Magnetic_Fieldd(time, Bzdot, Brdot, Btdot, mean_Cutoff, max_Range,
                       rprobe_Dia = 0.003175, ztprobe_Dia = 0.00158755,
                       hole_Sep = 0.001016):
    """
        Outputs the magnetic field strengths, time and monitor value.
    """
   
    #Setting the area variables
    rprobe_Area = np.pi*(rprobe_Dia/2)**2
    tzprobe_Area = ztprobe_Dia*hole_Sep
    
    #Setting the time variables
    time_Microsec = time
    time_Sec = time_Microsec*1e-6
    timeB_Sec = time_Sec[1:]
    
    # My standard is {Z, R, Theta}
    Bzdot -= np.mean(Bzdot[0:mean_Cutoff])
    Brdot -= np.mean(Brdot[0:mean_Cutoff])
    Btdot -= np.mean(Btdot[0:mean_Cutoff])
    
    Bzdot = BMX_Clear_Infs(Bzdot, max_Range)
    Brdot = BMX_Clear_Infs(Brdot, max_Range)
    Btdot = BMX_Clear_Infs(Btdot, max_Range)
    
    # Performing the integration using the trap-rule (Guass units)
    Bz = sp.cumtrapz(Bzdot/(tzprobe_Area),time_Sec)*1e4
    Br = sp.cumtrapz(Brdot/(rprobe_Area),time_Sec)*1e4
    Bt = sp.cumtrapz(Btdot/(tzprobe_Area),time_Sec)*1e4
    
    Bz = sps.detrend(Bz)
    Br = sps.detrend(Br)
    Bt = sps.detrend(Bt)
    
    return  timeB_Sec, Bz, Br, Bt

def Windowed_Fluctuation_Amp(signal, domain, increment_Index, save = True, indicator = '(shot #)'):

    """
       This function is a window scan function. It calculates the mean value, unnormalized 
       fluctuation amplitude and a domain estimate for a specific domain range. The domain range is 
       specified by the "increment_Index" parameter. The data is usually saved as a dictionary. 

       If save = False, the function will output the data in this order: std_Values, mean_Values,
       time_Estimate, time_Range. ** In the code I have domain instead of time, I was feeling general 
       when writting the code.
    """
    ################################################################################################
    # The following lines are the variables
    ################################################################################################
    start_Index = 0
    end_Index = start_Index + increment_Index
    mean_Values = []
    std_Values = []
    domain_Ranges = []
    domain_Estimate = []
    ################################################################################################
    # The while loop contains the meat of the code.
    ################################################################################################
    while (end_Index <= (signal.size - 1)):
        mean_Values.append(np.abs(np.mean(signal[start_Index:end_Index])))
        std_Values.append(np.std(signal[start_Index:end_Index]))
        domain_Estimate.append(np.mean(domain[start_Index:end_Index]))
        domain_Ranges.append((domain[start_Index], domain[end_Index]))

        start_Index += 1
        end_Index += 1
    ################################################################################################
    # Converting the lists into numpy arrays.
    ################################################################################################
    std_Values = np.asarray(std_Values)
    mean_Values = np.asarray(mean_Values)
    domain_Ranges = np.asarray(domain_Ranges)
    domain_Estimate = np.asarray(domain_Estimate)
    ################################################################################################
    # The lines below are the dictionary definition.
    ################################################################################################
    if save :
        my_dict = {'STD' : std_Values, 'Mean' : mean_Values,
                     'Window Domain' : domain_Estimate, 'Window Range' : domain_Ranges}
        filesave = 'windo_Fluct_Amp_' + indicator + '.npy'
        print('Saving ' +  filesave)
        np.save(filesave, my_dict)
    else:
        return std_Values, mean_Values, domain_Estimate, domain_Ranges

def which_Indictor(pico_Num, shot_Num):
    ##############################################################################################
    """ The function creates the indicator for the general window_Scan """
    ##############################################################################################    
    if (pico_Num == 1):
        indicator_1 = '(Pos1S' + str(shot_Num) + ')'
        indicator_2 = '(Pos3S' + str(shot_Num) + ')'
    elif (pico_Num == 2):
        indicator_1 = '(Pos5S' + str(shot_Num) + ')'
        indicator_2 = '(Pos7S' + str(shot_Num) + ')'
    elif (pico_Num == 3):
        indicator_1 = '(Pos9S' + str(shot_Num) + ')'
        indicator_2 = '(Pos11S' + str(shot_Num) + ')'
    elif (pico_Num == 4):
        indicator_1 = '(Pos13S' + str(shot_Num) + ')'
        indicator_2 = '(Volt_Cur' + str(shot_Num) + ')'
    else:
        print('The pico number should be 1, 2, 3, or 4.')

    return indicator_1, indicator_2

def which_File(pos, shot):

        if (pos == 1 or pos == 3):
            pico = which_Pico(pos)
            indicator1, indicator2 = which_Indictor(pico, shot)
            if (pos == 1):
                indicator = indicator1
            else:
                indicator = indicator2
        elif (pos == 5 or pos == 7):
            pico = which_Pico(pos)
            indicator1, indicator2 = which_Indictor(pico, shot)
            if (pos == 5):
                indicator = indicator1
            else:
                indicator = indicator2
        elif (pos == 9 or pos == 11):
            pico = which_Pico(pos)
            indicator1, indicator2 = which_Indictor(pico, shot)
            if (pos == 9):
                indicator = indicator1
            else:
                indicator = indicator2
        else:
            pico = which_Pico(pos)
            indicator1, indicator2 = which_Indictor(pico, shot)
            indicator = indicator1
        
        return indicator

def which_Pico(pos):
        
        if (pos == 1 or pos == 3):
            pico = 1
        elif (pos == 5 or pos == 7):
            pico = 2
        elif (pos == 9 or pos == 11):
            pico = 3
        else:
            pico = 4
        
        return pico
    
def zero_Crossings(ydata, xdata):
    ##############################################################################################
    """ This functions name is self explainatory.
        ::::Both inputs are 1d-numpy arrays with the same size.    """
    ##############################################################################################    
    zero_Crossing = []
    index_Crossing = []
    for i in range(0, ydata.size - 1):
        if xdata[i] > 0.0:
            if ((ydata[i] > 0) and (ydata[i+1] < 0)) or ((ydata[i] < 0 ) and (ydata[i+1])):
                zero_Crossing.append(xdata[i])
                index_Crossing.append(i)
            else: 
                pass

        else:
            pass

    return zero_Crossing, index_Crossing

def finding_Index_Time(time_Data, time_Interest):
    ##############################################################################################
    """ This function returns the index corresponding to a given time. """
    ##############################################################################################
    """data_Path = '04232019/04232019pico1/'
    filename = data_Path + '20190423-(1).npy'
    print('Loading ' + filename)
    data_Dict = np.load(filename)"""
    """
    my_Dict = {'Time_B': timeB_Sec, 'Bz': Bz, 'Bt': Bt, 
                            'Bz_prime': Bz_prime, 'Bt_prime': Bt_prime}
    """
    time = time_Data
    output_Time = 0
    for i in range(0,time.size - 1):
        if (time[i]*1e6 <= time_Interest):
            output_Time = i

    #print('Outputing the index of interest')
    return output_Time