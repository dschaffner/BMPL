#data_conversion_04232019.py

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:34:02 2019

@author: dschaffner
"""
import numpy as np
import scipy.integrate as sp
import pandas as pd
import magnetic_pickup_probe as mpp
import matplotlib.pylab as plt
import astropy.units as u

datadirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2024\\07112024\\'

#numskipped=len(skipshots)
#numshots=maxshot-(startshot-1)-numskipped

#runfilename = 'Dataset4_wiretarget_1kV_1msstuff_6kAwire_40shots' #shots 1 to 63 skip 3,4,5,8,22,23,24,29,30,33,34,35,39,40,46,47,48,54,55,56,57,60,61
#runfilename = 'Dataset4_wiretarget_1p5kV_1msstuff_0kAwire_40shots' #shots 64 to 110 skip 70,71,91,92,95,96
runfilename = 'Dataset4_wiretarget_1p5kV_1msstuff_8kAwire_40shots' #shots 111 to 151 skip 128

"""
5 probes:
    Triplet Bdot at Port 21, 23, 25 all at axial center
    Squid Probe at Ports 27 and 29 (only one loop in z direction)
    Electric field probe at opposite side at port 25 at axial center (postive efield points away from gun)
    
    Bdots at 500mV/div on scope
"""

from scipy.signal import butter, sosfiltfilt

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y
    
# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 50e6

startintg_index=0#3000
meancutoff = 1000
bdotmaxrange=500#mV

probe_dia = 0.00158755#m (1/16'' probe)
hole_sep = 0.001016#m (1/16''probe)
rloop_area = np.pi*((probe_dia/2)**2)
rloop_area = rloop_area * u.m**2
tloop_area = probe_dia*hole_sep * u.m**2
zloop_area = probe_dia*hole_sep * u.m**2

testpicoshot=pd.read_csv(datadirectory+'\pico1\\20240711-0001 ('+str(6)+').csv',header=[0],skiprows=[1])

tstart = testpicoshot['Time'][0]*1e-6#pico reads out in us, convert to seconds
tinterval=np.abs((testpicoshot['Time'][1]-testpicoshot['Time'][0])*1e-6)
numsamples=int(np.shape(testpicoshot['Time'])[0])
time_s = tstart+(np.arange(numsamples)*tinterval)

time_s = time_s * u.s






###BEGIN Scan READIN#####
numshots=40#30
numprobes_trip = 3
numprobes_squid = 2
numpoints_squid = 6
Discharge_raw=np.zeros([numshots,numsamples])
HV_raw=np.zeros([numshots,numsamples])
Discharge=np.zeros([numshots,numsamples])
HV=np.zeros([numshots,numsamples])
#wirecurrent=np.zeros([num_rad_pos,num_currents,numshots,numsamples_wc])
Bdotr = np.zeros([numprobes_trip,numshots,numsamples])
Bdott = np.zeros([numprobes_trip,numshots,numsamples])
Bdotz = np.zeros([numprobes_trip,numshots,numsamples])
Br = np.zeros([numprobes_trip,numshots,numsamples-1])
Bt = np.zeros([numprobes_trip,numshots,numsamples-1])
Bz = np.zeros([numprobes_trip,numshots,numsamples-1])
Bsquid_Bdotz = np.zeros([numprobes_squid,numpoints_squid,numshots,numsamples])
Bsquid_Bz = np.zeros([numprobes_squid,numpoints_squid,numshots,numsamples-1])
efield = np.zeros([numshots,numsamples])
term50=2.0

pico=pd.read_csv(datadirectory+'\Pico1\\20240711-0001 ('+str(111)+').csv',header=[0],skiprows=[1])
pico.replace("-∞",-bdotmaxrange,inplace=True)
pico.replace("∞",bdotmaxrange,inplace=True)


bdottest1=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
bdottest2=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
btest1= sp.cumtrapz(bdottest1/rloop_area,time_s,initial=0)*1e4#Gauss
btest2= sp.cumtrapz(bdottest2/tloop_area,time_s,initial=0)*1e4#Gauss

test=np.array(pico['Channel C'],dtype=float)
bdottest1_units = test * term50 * u.millivolt
btest1f=mpp.compute_bfield(bdottest1_units.to(u.volt),time_s,rloop_area,num_loop=1,gain=1)


plt.plot(time_s,btest1)
plt.plot(time_s,btest1f.to('Gauss'))


"""
#shots 111 to 151 skip 128
startshot=111
maxshot=152
skipshots=[128]
### Extract PICO 1 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico1\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 1: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)

    HV_raw[savenumber,:]=np.array(pico['Channel A'],dtype=float)
    Discharge_raw[savenumber,:]=np.array(pico['Channel B'],dtype=float)
    Discharge[savenumber,:]=Discharge_raw[savenumber,:]*1000.0*2.0
    HV[savenumber,:]=HV_raw[savenumber,:]*1000.0
    
    Bdotr[0,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    Bdott[0,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bdotr[0,savenumber,:]=butter_bandpass_filter(Bdotr[0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[0,savenumber,:]=butter_bandpass_filter(Bdott[0,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br[0,savenumber,:]= sp.cumtrapz(Bdotr[0,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[0,savenumber,:]= sp.cumtrapz(Bdott[0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO1 READIN ##### 
   
### Extract PICO 2 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico2\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 2: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotz[0,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50*1e-3#to convert to volts
    Bdotr[1,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50*1e-3#to convert to volts
    Bdott[1,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    Bdotz[1,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bdotz[0,savenumber,:]=butter_bandpass_filter(Bdotz[0,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotr[1,savenumber,:]=butter_bandpass_filter(Bdotr[1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[1,savenumber,:]=butter_bandpass_filter(Bdott[1,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[1,savenumber,:]=butter_bandpass_filter(Bdotz[1,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bz[0,savenumber,:]= sp.cumtrapz(Bdotz[0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Br[1,savenumber,:]= sp.cumtrapz(Bdotr[1,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[1,savenumber,:]= sp.cumtrapz(Bdott[1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[1,savenumber,:]= sp.cumtrapz(Bdotz[1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO2 READIN #####
    
### Extract PICO 3 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico3\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 3: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bdotr[2,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50*1e-3#to convert to volts
    Bdott[2,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50*1e-3#to convert to volts
    Bdotz[2,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[0,0,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bdotr[2,savenumber,:]=butter_bandpass_filter(Bdotr[2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdott[2,savenumber,:]=butter_bandpass_filter(Bdott[2,savenumber,:],lowcut,highcut,fs,order=9)
    Bdotz[2,savenumber,:]=butter_bandpass_filter(Bdotz[2,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[0,0,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,0,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Br[2,savenumber,:]= sp.cumtrapz(Bdotr[2,savenumber,:]/rloop_area,time_s)*1e4#Gauss
    Bt[2,savenumber,:]= sp.cumtrapz(Bdott[2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bz[2,savenumber,:]= sp.cumtrapz(Bdotz[2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[0,0,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO3 READIN #####

### Extract PICO 4 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico4\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 4: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bsquid_Bdotz[0,1,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[0,2,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[0,3,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[0,4,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bsquid_Bdotz[0,1,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,1,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[0,2,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,2,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[0,3,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,3,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[0,4,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,4,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bsquid_Bz[0,1,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[0,2,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[0,3,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,3,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[0,4,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,4,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO4 READIN #####

### Extract PICO 5 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico5\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 5: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bsquid_Bdotz[0,5,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[1,0,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[1,1,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[1,2,savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bsquid_Bdotz[0,5,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[0,5,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[1,0,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,0,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[1,1,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,1,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[1,2,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,2,savenumber,:],lowcut,highcut,fs,order=9)
    
    #Compute Magnetic Field for Bdots
    Bsquid_Bz[0,5,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[0,5,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[1,0,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,0,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[1,1,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,1,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[1,2,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,2,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    #End of loop
    savenumber+=1
##END PICO5 READIN #####

### Extract PICO 6 ###
savenumber=0
for shot in np.arange(startshot,maxshot):
    if shot in skipshots:
        continue
    #Load Picoscope file
    pico=pd.read_csv(datadirectory+'\Pico6\\20240711-0001 ('+str(shot)+').csv',header=[0],skiprows=[1])
    print('Pico 6: loading shot number ', shot)
    
    #replace infinities (from clipping) with max range values
    pico.replace("-∞",-bdotmaxrange,inplace=True)
    pico.replace("∞",bdotmaxrange,inplace=True)
    
    Bsquid_Bdotz[1,3,savenumber,:]=np.array(pico['Channel A'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[1,4,savenumber,:]=np.array(pico['Channel B'],dtype=float)*term50*1e-3#to convert to volts
    Bsquid_Bdotz[1,5,savenumber,:]=np.array(pico['Channel C'],dtype=float)*term50*1e-3#to convert to volts
    
    #filter Bdot
    Bsquid_Bdotz[1,3,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,3,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[1,4,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,4,savenumber,:],lowcut,highcut,fs,order=9)
    Bsquid_Bdotz[1,5,savenumber,:]=butter_bandpass_filter(Bsquid_Bdotz[1,5,savenumber,:],lowcut,highcut,fs,order=9)
 
    #Compute Magnetic Field for Bdots
    Bsquid_Bz[1,3,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,3,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[1,4,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,4,savenumber,:]/tloop_area,time_s)*1e4#Gauss
    Bsquid_Bz[1,5,savenumber,:]= sp.cumtrapz(Bsquid_Bdotz[1,5,savenumber,:]/tloop_area,time_s)*1e4#Gauss

    efield[savenumber,:]=np.array(pico['Channel D'],dtype=float)*term50*200.0*(0.001)/0.0025#in V/m
    #End of loop
    savenumber+=1
##END PICO6 READIN #####






#create magnetic probe group
mag = h5f.create_group("mag_probe")
r=mag.create_group("r")
r.create_dataset("bdot",data=Bdotr)
r.create_dataset("b",data=Br)
t=mag.create_group("t")
t.create_dataset("bdot",data=Bdott)
t.create_dataset("b",data=Bt)
z=mag.create_group("z")
z.create_dataset("bdot",data=Bdotz)
z.create_dataset("b",data=Bz)
squid=mag.create_group("squid")
squid.create_dataset("bdot",data=Bsquid_Bdotz)
squid.create_dataset("b",data=Bsquid_Bz)

#create wirecurrent probe group
#wc=h5f.create_group("wirecurrent")
#wc.create_dataset("wirecurrent",data=wirecurrent)
 
#create floating potential probe group
efieldp=h5f.create_group("efield_probe")
efieldp.create_dataset("efield",data=efield)
                      
#create discharage group
dis=h5f.create_group("discharge")
discharge_current=dis.create_dataset('dis_I',data=Discharge)
discharge_voltage=dis.create_dataset('dis_V',data=HV)
discharge_current_raw=dis.create_dataset('dis_I_raw',data=Discharge_raw)
discharge_voltage_raw=dis.create_dataset('dis_V_raw',data=HV_raw)

h5f.close()

"""
