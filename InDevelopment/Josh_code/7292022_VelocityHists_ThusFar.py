# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:26:55 2022

@author: Josh0
"""

# This purpose of the code is to accumulate all velocity data thus far.  It will be organized by date.  There will be a key as well here below.

from scipy.stats import norm
import numpy as np
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from scipy import stats

#plt.style.use('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')
#plt.style.use('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')
#plt.style.use('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots_new.py')
plt.style.use('C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots_new.py')

arrival_data_directory_location = 'C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\Arrival_time_velocities\\'
timeavg_data_directory_location = 'C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\Time_avgd_velocities\\'
data_directory_location_9232022 = 'C:\\Users\\Josh0\\OneDrive\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Data\\2022\\Time_avgd_velocities\\9232022\\'
data_directory_location_9232022 = 'C:\\Users\\dschaffner\\Documents\\GitHub\\BMPL\\InDevelopment\\Josh_code\\Velocity_nozzle_disV_FullCatalogue\\'
# I don't need anything from the 2232022_3172022 .txt file, I split them up into individual .txt files.
# I'm also not doing anything with the 4kV data - more of a proof of concept than actual data.

"""
************************
Below are the datafilenames for Pos5 / Pos7.  We NOW HAVE datafiles for 2.5kV Discharge at 110mT fortunately.  The delay issue corrected itself at shot 210 in the 3172022 dataset.
************************
"""

tavg_data1_57 = '3172022_Shot208_253_Pos5and7_Tavgd'
tavg_data2_57 = '6202022_110mT_3kV_Pos5and7_Tavgd'
tavg_data3_57 = '6202022_110mT_3p5kV_Pos5and7_Tavgd'
tavg_data4_57 = '7152022_110mT_3p5kV_Pos5and7_Tavgd'
tavg_data5_57 = '7202022_120mT_3p5kV_Pos5and7_Tavgd'
tavg_data6_57 = '7202022_133mT_3p5kV_Pos5and7_Tavgd'
tavg_data7_57 = '7202022_17mT_3p5kV_Pos5and7_Tavgd'
tavg_data8_57_delayed = '3172022_Delayed_110mT_Shot57to179_Pos5and7_Tavgd'
tavg_data9_57 = '9232022_0mT_3p5kV_pos57'
tavg_data9_68 = '9232022_0mT_3p5kV_pos68'

tavg_pos57_317 = np.loadtxt(timeavg_data_directory_location + '3172022_Shot208_253_Pos5and7_Tavgd.txt', unpack=True)

#tavg_pos57_317 = np.loadtxt(timeavg_data_directory_location + tavg_data1_57 + '.txt', unpack=True)
tavg_pos57_317_delayed = np.loadtxt(timeavg_data_directory_location + tavg_data8_57_delayed + '.txt', unpack=True)
tavg_pos57_620a = np.loadtxt(timeavg_data_directory_location + tavg_data2_57 + '.txt', unpack=True)
tavg_pos57_620b = np.loadtxt(timeavg_data_directory_location + tavg_data3_57 + '.txt', unpack=True)
tavg_pos57_715 = np.loadtxt(timeavg_data_directory_location + tavg_data4_57 + '.txt', unpack=True)
tavg_pos57_720a = np.loadtxt(timeavg_data_directory_location + tavg_data5_57 + '.txt', unpack=True)
tavg_pos57_720b = np.loadtxt(timeavg_data_directory_location + tavg_data6_57 + '.txt', unpack=True)
tavg_pos57_720c = np.loadtxt(timeavg_data_directory_location + tavg_data7_57 + '.txt', unpack=True)
tavg_pos57_923 = np.loadtxt(data_directory_location_9232022 + tavg_data9_57 + '.txt', unpack=True)
tavg_pos68_923 = np.loadtxt(data_directory_location_9232022 + tavg_data9_68 + '.txt', unpack=True)

tavg_p57_317_55mT_2p5kV = np.array(tavg_pos57_317[:])
tavg_p57_317_110mT_2p5kV = np.array(tavg_pos57_317_delayed[:])
tavg_p57_620_110mT_3kV = np.array(tavg_pos57_620a[:])
tavg_p57_620_110mT_3p5kV = np.array(tavg_pos57_620b[:])
tavg_p57_715_110mT_3p5kV = np.array(tavg_pos57_715[:])
tavg_p57_720_120mT_3p5kV = np.array(tavg_pos57_720a[:])
tavg_p57_720_133mT_3p5kV = np.array(tavg_pos57_720b[:])
tavg_p57_720_17mT_3p5kV = np.array(tavg_pos57_720c[:])  # *********This is a unevenly nozzled; one coil at 5ms, other still at max*************** #
tavg_p57_923_0mT_3p5kV = np.array(tavg_pos57_923[:])
tavg_p68_923_0mT_3p5kV = np.array(tavg_pos68_923[:])

# Infinity removal if needed below.
"""
maxrange=1
posinfs = np.isposinf(arr_vel_data620_110mT_3kV)
arr_vel_data620_110mT_3kV[np.where(posinfs)] = maxrange
"""

(mu57_1, sigma57_1) = norm.fit(tavg_p57_317_55mT_2p5kV)
(mu57_2, sigma57_2) = norm.fit(tavg_p57_620_110mT_3kV)
(mu57_3, sigma57_3) = norm.fit(tavg_p57_620_110mT_3p5kV)
(mu57_4, sigma57_4) = norm.fit(tavg_p57_715_110mT_3p5kV)
(mu57_5, sigma57_5) = norm.fit(tavg_p57_720_120mT_3p5kV)
(mu57_6, sigma57_6) = norm.fit(tavg_p57_720_133mT_3p5kV)
(mu57_7, sigma57_7) = norm.fit(tavg_p57_720_17mT_3p5kV) # **** BAD DATA!****
(mu57_8, sigma57_8) = norm.fit(tavg_p57_317_110mT_2p5kV)
(mu57_9, sigma57_9) = norm.fit(tavg_p57_923_0mT_3p5kV)
(mu68, sigma68) = norm.fit(tavg_p68_923_0mT_3p5kV)

vavg57_1 = round(mu57_1)
vavg57_2 = round(mu57_2)
vavg57_3 = round(mu57_3)
vavg57_4 = round(mu57_4)
vavg57_5 = round(mu57_5)
vavg57_6 = round(mu57_6)
vavg57_7 = round(mu57_7) # **** BAD DATA!****
vavg57_8 = round(mu57_8)
vavg57_9 = round(mu57_9)
vavg68 = round(mu68)

errvavg57_1 = round(sigma57_1)
errvavg57_2 = round(sigma57_2)
errvavg57_3 = round(sigma57_3)
errvavg57_4 = round(sigma57_4)
errvavg57_5 = round(sigma57_5)
errvavg57_6 = round(sigma57_6)
errvavg57_7 = round(sigma57_7) # **** BAD DATA!****
errvavg57_8 = round(sigma57_8)
errvavg57_9 = round(sigma57_9)
errvavg68 = round(sigma68)


vavg57_array = np.array([vavg57_1, vavg57_2, vavg57_3, vavg57_4, vavg57_5, vavg57_6, vavg57_7, vavg57_1, vavg57_9])#vavg57_8, vavg57_9])
err57_array = np.array([errvavg57_1, errvavg57_2, errvavg57_3, errvavg57_4, errvavg57_5, errvavg57_6, errvavg57_7, errvavg57_8, errvavg57_9])

vavg68_array = np.array([vavg68])
err68_array = np.array([errvavg68])

tavg57_nozzle_velarray = ([vavg57_array[8], np.mean([vavg57_array[2],vavg57_array[3]]), vavg57_array[4], vavg57_array[5]])
err_tavg57_nozzle_array = ([err57_array[8], np.mean([err57_array[2], err57_array[3]]), err57_array[4], err57_array[5]])
#We have disV of 2.5kV, 3kV, and 3.5kV as of now, for shots at 110mT nozzle field.
tavg57_disvolt_velarray = ([vavg57_array[7], vavg57_array[1], np.mean([vavg57_array[2],vavg57_array[3]])])
err_tavg57_disvolt_array = ([err57_array[0], err57_array[1], np.mean([err57_array[2],err57_array[3]])])

nozzle_68 = np.array([0])
DisV_3p5kV_array = np.append(tavg_p57_620_110mT_3p5kV, tavg_p57_715_110mT_3p5kV)

(mu57_10, sigma57_10) = norm.fit(DisV_3p5kV_array)
vavg57_10 = round(mu57_10)
errvavg57_10 = round(sigma57_10)


"""
************************
The datafilenames below are all Pos19 / Pos21.  The Pos5 / Pos7 datafilenames are above.
************************
"""
arrival_datafilename1 = 'InnerCoilOnly_1122022'
arrival_datafilename2 = 'arrival_shot57to179_MaxNozzle_3172022'
arrival_datafilename3 = 'arrival_shot184to253_0p5MaxNozzle_3172022'
arrival_datafilename4 = '6202022_Shots31to51_3p5kV_Delays'
arrival_datafilename5 = '6202022_Shots11to30_3kV_Delays'
arrival_datafilename6 = '7152022_AllCoils3p5kV_ToFDelays'
arrival_datafilename7 = '7202022_225VNozzle_delays'
arrival_datafilename8 = '7202022_250VNozzle_delays'

# These are already in km/s and are not delays.
tavg_datafilename1 = '1122022_0mT_2kV_Discharge_Tavgd'
tavg_datafilename2 = '3172022_110mT_2p5kV_Discharge_Tavgd'
tavg_datafilename3 = '3172022_55mT_2p5kV_Discharge_Tavgd'
tavg_datafilename4 = '6202022_110mT_3kV_Discharge_Tavgd'
tavg_datafilename5 = '6202022_110mT_3p5kV_Discharge_Tavgd'
tavg_datafilename6 = '7152022_110mT_3p5kV_Discharge_Tavgd'
tavg_datafilename7 = '7202022_120mT_3p5kV_Discharge_Tavgd'
tavg_datafilename8 = '7202022_133mT_3p5kV_Discharge_Tavgd'
tavg_datafilename9 = '7202022_17mT_3p5kV_Discharge_Tavgd'

arr_data112 = np.loadtxt(arrival_data_directory_location + arrival_datafilename1 + '.txt', skiprows=2, unpack=True)
arr_data317a = np.loadtxt(arrival_data_directory_location + arrival_datafilename2 + '.txt', skiprows=2, unpack=True)
arr_data317b = np.loadtxt(arrival_data_directory_location + arrival_datafilename3 + '.txt', skiprows=2, unpack=True)
arr_data620a = np.loadtxt(arrival_data_directory_location + arrival_datafilename4 + '.txt', unpack=True)
arr_data620b = np.loadtxt(arrival_data_directory_location + arrival_datafilename5 + '.txt', unpack=True)
arr_data715 = np.loadtxt(arrival_data_directory_location + arrival_datafilename6 + '.txt', unpack=True)
arr_data720a = np.loadtxt(arrival_data_directory_location + arrival_datafilename7 + '.txt', unpack=True)
arr_data720b = np.loadtxt(arrival_data_directory_location + arrival_datafilename8 + '.txt', unpack=True)


tavg_data112 = np.loadtxt(timeavg_data_directory_location + tavg_datafilename1 + '.txt', unpack=True)
tavg_data317a = np.loadtxt(timeavg_data_directory_location + tavg_datafilename2 + '.txt', unpack=True)
tavg_data317b = np.loadtxt(timeavg_data_directory_location + tavg_datafilename3 + '.txt', unpack=True)
tavg_data620a = np.loadtxt(timeavg_data_directory_location + tavg_datafilename4 + '.txt', unpack=True)
tavg_data620b = np.loadtxt(timeavg_data_directory_location + tavg_datafilename5 + '.txt', unpack=True)
tavg_data715 = np.loadtxt(timeavg_data_directory_location + tavg_datafilename6 + '.txt', unpack=True)
tavg_data720a = np.loadtxt(timeavg_data_directory_location + tavg_datafilename7 + '.txt', unpack=True)
tavg_data720b = np.loadtxt(timeavg_data_directory_location + tavg_datafilename8 + '.txt', unpack=True)
tavg_data720c = np.loadtxt(timeavg_data_directory_location + tavg_datafilename9 + '.txt', unpack=True)

arr_data112_0mT = np.array(arr_data112[100:200])
arr_data317_55mT = np.array(arr_data317b[:])
arr_data317_110mT = np.array(arr_data317a[:])
arr_data620_110mT_3kV = np.array(arr_data620a[:])
arr_data620_110mT_3p5kV = np.array(arr_data620b[:])
arr_data715_110mT_3p5kV = np.array(arr_data715[:])
arr_data720_120mT_3p5kV = np.array(arr_data720a[:])
arr_data720_133mT_3p5kV = np.array(arr_data720b[:])

tavg_data112_0mT_2kV = np.array(tavg_data112[:])
tavg_data317_55mT_2p5kV = np.array(tavg_data317b[:])
tavg_data317_110mT_2p5kV = np.array(tavg_data317a[:])
tavg_data620_110mT_3kV = np.array(tavg_data620a[:])
tavg_data620_110mT_3p5kV = np.array(tavg_data620b[:])
tavg_data715_110mT_3p5kV = np.array(tavg_data715[:])
tavg_data720_120mT_3p5kV = np.array(tavg_data720a[:])
tavg_data720_133mT_3p5kV = np.array(tavg_data720b[:])
tavg_data720_5ms_3p5kV = np.array(tavg_data720c[:])


arr_vel_data112_0mT_2kV = arr_data112_0mT
arr_vel_data317_55mT_2p5kV = (0.026/arr_data317_55mT)*1000
arr_vel_data317_110mT_2p5kV = (0.026/arr_data317_110mT)*1000
arr_vel_data620_110mT_3kV = (0.026/arr_data620_110mT_3kV)*1000
arr_vel_data620_110mT_3p5kV = (0.026/arr_data620_110mT_3p5kV)*1000
arr_vel_data715_110mT_3p5kV = (0.026/arr_data715_110mT_3p5kV)*1000
arr_vel_data720_120mT_3p5kV = (0.026/arr_data720_120mT_3p5kV)*1000
arr_vel_data720_133mT_3p5kV = (0.026/arr_data720_133mT_3p5kV)*1000


maxrange=1
posinfs = np.isposinf(arr_vel_data620_110mT_3kV)
arr_vel_data620_110mT_3kV[np.where(posinfs)] = maxrange

(mu_112arr, sigma_112arr) = norm.fit(arr_vel_data112_0mT_2kV) #1
(mu_317arr1, sigma_317arr1) = norm.fit(arr_vel_data317_55mT_2p5kV) #2
(mu_317arr2, sigma_317arr2) = norm.fit(arr_vel_data317_110mT_2p5kV) #3
(mu_620arr1, sigma_620arr1) = norm.fit(arr_vel_data620_110mT_3kV) #4
(mu_620arr2, sigma_620arr2) = norm.fit(arr_vel_data620_110mT_3p5kV) #5
(mu_715arr, sigma_715arr) = norm.fit(arr_vel_data715_110mT_3p5kV) #6
(mu_720arr1, sigma_720arr1) = norm.fit(arr_vel_data720_120mT_3p5kV) #7
(mu_720arr2, sigma_720arr2) = norm.fit(arr_vel_data720_133mT_3p5kV) #8

(mu_112tavg, sigma_112tavg) = norm.fit(tavg_data112_0mT_2kV) #9
(mu_317tavg1, sigma_317tavg1) = norm.fit(tavg_data317_55mT_2p5kV) #10
(mu_317tavg2, sigma_317tavg2) = norm.fit(tavg_data317_110mT_2p5kV) #11
(mu_620tavg1, sigma_620tavg1) = norm.fit(tavg_data620_110mT_3kV) #12
(mu_620tavg2, sigma_620tavg2) = norm.fit(tavg_data620_110mT_3p5kV) #13
(mu_715tavg, sigma_715tavg) = norm.fit(tavg_data715_110mT_3p5kV) #14
(mu_720tavg1, sigma_720tavg1) = norm.fit(tavg_data720_120mT_3p5kV) #15
(mu_720tavg2, sigma_720tavg2) = norm.fit(tavg_data720_133mT_3p5kV) #16
(mu_720tavg3, sigma_720tavg3) = norm.fit(tavg_data720_5ms_3p5kV) #17


Vavg1 = round(mu_112arr)
Vavg2 = round(mu_317arr1)
Vavg3 = round(mu_317arr2)
Vavg4 = round(mu_620arr1)
Vavg5 = round(mu_620arr2)
Vavg6 = round(mu_715arr)
Vavg7 = round(mu_720arr1)
Vavg8 = round(mu_720arr2)
Vavg9 = round(mu_112tavg)
Vavg10 = round(mu_317tavg1)
Vavg11 = round(mu_317tavg2)
Vavg12 = round(mu_620tavg1)
Vavg13 = round(mu_620tavg2)
Vavg14 = round(mu_715tavg)
Vavg15 = round(mu_720tavg1)
Vavg16 = round(mu_720tavg2)
Vavg17 = round(mu_720tavg3)

err1 = round(sigma_112arr)
err2 = round(sigma_317arr1)
err3 = round(sigma_317arr2)
err4 = round(sigma_620arr1)
err5 = round(sigma_620arr2)
err6 = round(sigma_715arr)
err7 = round(sigma_720arr1)
err8 = round(sigma_720arr2)
err9 = round(sigma_112tavg)
err10 = round(sigma_317tavg1)
err11 = round(sigma_317tavg2)
err12 = round(sigma_620tavg1)
err13 = round(sigma_620tavg2)
err14 = round(sigma_715tavg)
err15 = round(sigma_720tavg1)
err16 = round(sigma_720tavg2)
err17 = round(sigma_720tavg3)

Vavg_array = np.zeros(17)

# I'll brute force the array:
Vavg_array = np.array([Vavg1, Vavg2, Vavg3, Vavg4, Vavg5, Vavg6, Vavg7, Vavg8, Vavg9, Vavg10, Vavg11, Vavg12, Vavg13, Vavg14, Vavg15, Vavg16, Vavg17])
Vavg_arrival = Vavg_array[0:8]    
Vavg_tavg = Vavg_array[8:18]

err_array = np.array([err9, err10, err11, err12, err13, err14, err15, err16, err17])

# We need to average #5 and #6 since they're the same dataset (and #13 and #14)

#print(vavg_arrival)
#print(vavg_tavg)

# How do I want to depict this data...

nozzle_array = np.array([0, 55, 110, 110, 110, 110, 120, 133])
discharge_voltages = np.array([2, 2.5, 2.5, 3, 3.5, 3.5, 3.5, 3.5])
arrival_bins = np.array(range(50, 500, 10))
tavg_bins = np.array(range(20, 100, 10))
bins1 = np.array(range(20, 250, 12))

dis_voltages1 = np.array([2.5, 3, 3.5])#kV #This is at constant 110mT nozzle field.
#dis_voltages2 = np.array([3, 3.5])
nozzle_fields1 = np.array([110, 120, 133])
nozzle_fields2 = np.array([0, 110, 120, 133])#mT #This is at consant 3.5kV discharge.

arrival_disvolt_VelocityArray = ([Vavg_arrival[2], Vavg_arrival[3], np.mean([Vavg_arrival[4],Vavg_arrival[5]])])
arrival_nozzle_VelocityArray = ([np.mean([Vavg_arrival[4],Vavg_arrival[5]]), Vavg_arrival[6], Vavg_arrival[7]])
tavg_disvolt_VelocityArray = ([Vavg_tavg[2], Vavg_tavg[3], np.mean([Vavg_tavg[4],Vavg_tavg[5]])])
tavg_nozzle_VelocityArray = ([Vavg_tavg[8], np.mean([Vavg_tavg[4],Vavg_tavg[5]]), Vavg_tavg[6], Vavg_tavg[7]])

err_tavg_disvolt1921 = ([err_array[2], err_array[3], np.mean([err_array[4],err_array[5]])])
err_tavg_nozzle1921 = ([err_array[8], np.mean([err_array[4],err_array[5]]), err_array[6], err_array[7]])

xticks1 = np.linspace(2, 4, 5)
xticks2 = np.array(range(5, 175, 25))
"""
fig1, (ax1, ax2) = plt.subplots(2)

ax1.hist(tavg_p57_720_17mT_3p5kV, bins=np.array(range(50, 300, 15)), align='left', edgecolor='black', color='green', label='3.5kV, 15mT Nozzle  --  $\mu$ = ' + str(vavg57_7) + ' km/s')
ax1.set_ylabel('Counts')
ax1.set_xlabel('Velocity (km/s)')
ax1.legend(loc='best', frameon=False)

ax2.hist(tavg_p57_720_133mT_3p5kV, bins=np.array(range(50, 300, 15)), align='left', edgecolor='black', color='red', label='3.5kV, 133mT Nozzle  --  $\mu$ = ' + str(vavg57_6) + ' km/s')
ax2.set_ylabel('Counts')
ax2.set_xlabel('Velocity (km/s)')
ax2.legend(loc='best', frameon=False)

#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\Pos57_3.5kV_15mTvs133mT_VelComparison.png', dpi=300)
"""

"""
fig1,(ax1, ax2) = plt.subplots(2)

ax1.scatter(dis_voltages1, tavg57_disvolt_velarray, label='Tavg Pos5/7', color='brown', s=75)
ax1.errorbar(dis_voltages1, tavg57_disvolt_velarray, yerr=err_tavg57_disvolt_array, xerr=False, capsize=3, capthick=1.5, ecolor='black', elinewidth=1, ls='none')
ax1.set_title('110mT Nozzle', fontsize='small')
ax1.set_ylabel('Velocity (km/s)')
#ax1.set_xlabel('Discharge Voltage (kV)')
ax1.set_ylim(30, 130)
ax1.legend(loc='best', frameon=False)
for i,j in zip(dis_voltages1, tavg57_disvolt_velarray):
    ax1.text(i+0.01,j+7, '{} km/s'.format(j), fontsize='x-small')

ax2.scatter(dis_voltages1, tavg_disvolt_VelocityArray, label='Tavg Pos19/21', color='blue', s=75)
ax2.errorbar(dis_voltages1, tavg_disvolt_VelocityArray, yerr=err_tavg_disvolt1921, xerr=False, capsize=3, capthick=1.5, ecolor='black', elinewidth=1, ls='none')
ax2.set_ylabel('Velocity (km/s)')
ax2.set_xlabel('Discharge Voltage (kV)')
ax2.legend(loc='best', frameon=False)
ax2.set_ylim(30, 130)
for i,j in zip(dis_voltages1, tavg_disvolt_VelocityArray):
    ax2.text(i+0.01,j+5, '{} km/s'.format(j), fontsize='x-small')

#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\11182022_disvolt_sweep_57_1921.png', dpi=300)
plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\EPS_filesforPaper\\11182022_disvolt_sweep_57_1921.eps', dpi=600)
"""

"""
#fig1, (ax1, ax2, ax3) = plt.subplots(3)

#ax1.hist(tavg_p57_317_110mT_2p5kV, bins=bins1, label='2.5kV Discharge', align='left', edgecolor='black', color='teal')
ax1.hist(tavg_p57_317_55mT_2p5kV, bins=bins1, label='2.5kV Discharge', align='left', edgecolor='black', color='teal')
#ax1.axvline(x=vavg57_8, color='black', label='Mean = '+str(vavg57_8) + ' km/s', linestyle='--', markersize='13')
ax1.axvline(x=vavg57_1, color='black', label='Mean = '+str(vavg57_1) + ' km/s', linestyle='--', markersize='13')

ax1.set_title('110mT Nozzle', fontsize='small')
ax1.set_ylabel('Counts')
#ax1.set_ylim(0, 8)
ax1.legend(loc='best', frameon=False)

ax2.hist(tavg_p57_620_110mT_3kV, bins=bins1, label='3kV Discharge', align='left', edgecolor='black', color='teal')# color='yellowgreen')
ax2.axvline(x=vavg57_2, color='black', label='Mean = '+str(vavg57_2)+ ' km/s', linestyle='--', markersize='13')
ax2.set_ylabel('Counts')
#ax2.set_ylim(0, 15)
ax2.legend(loc='best', frameon=False)

ax3.hist(DisV_3p5kV_array, bins=bins1, label='3.5kV Discharge', align='left', edgecolor='black', color='teal')# color='yellowgreen')
ax3.axvline(x=vavg57_10, color='black', label='Mean = 98 km/s', linestyle='--', markersize='13')
ax3.set_ylabel('Counts')
#ax2.set_ylim(0, 15)
ax3.legend(loc='best', frameon=False)
ax3.set_xlabel('Time-Averaged Velocity (km/s)')

plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\EPS_filesforPaper\\11232022_disvolt_Distributions_p57.eps', dpi=600)
"""


fig1, (ax2) = plt.subplots(1)
ax2.scatter(nozzle_fields2, tavg57_nozzle_velarray, label='Tavg Pos5/7', color='red', s=75)
ax2.errorbar(nozzle_fields2, tavg57_nozzle_velarray, yerr=err_tavg57_nozzle_array, xerr=False, capsize=3, capthick=1.5, ecolor='black', elinewidth=1, ls='none')
ax2.scatter(nozzle_68, vavg68_array, label='Tavg Pos6/8', color='green', s=75, alpha=0.5)
ax2.errorbar(nozzle_68, vavg68_array, yerr=err68_array, xerr=False, capsize=3, capthick=1.5, ecolor='black', elinewidth=1, ls='none')
#ax2.scatter(nozzle_fields2, tavg_nozzle_VelocityArray, label='Tavg Pos19/21', color='green', s=50)
#ax2.errorbar(nozzle_fields2, tavg_nozzle_VelocityArray, yerr=err_tavg_nozzle1921, xerr=False, capsize=3, capthick=1.5, ecolor='green', elinewidth=1, ls='none')
ax2.set_title('3.5kV Discharge', fontsize='small')
ax2.set_ylabel('Velocity (km/s)')
ax2.set_xlabel('Nozzle Field (mT)')
ax2.set_ylim(30, 175)
ax2.set_xlim(-25, 165)
ax2.legend(loc='best', frameon=False)
#for i,j in zip(nozzle_fields2, tavg_nozzle_VelocityArray):
#    ax2.text(i+1,j-10, '{} km/s'.format(j), fontsize='x-small')
for i,j in zip(nozzle_fields2, tavg57_nozzle_velarray):
    ax2.text(i+5,j-5, '{} km/s (Pos5/7)'.format(j), fontsize='x-small')
for i,j in zip(nozzle_68, vavg68_array):
    ax2.text(i+5,j+5, '{} km/s (Pos6/8)'.format(j), fontsize='x-small')

#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\9232022_nozzlefield_sweep_pos57_68.png', dpi=300)
#plt.savefig('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\EPS_filesforPaper\\9232022_nozzlefield_sweep_pos57_68.eps', dpi=600)





"""
This code below is a dual plot of both arrival and time averaged velocities, plotted as functions of Nozzle field and discharge voltage.

fig1, (ax1, ax2) = plt.subplots(2)


#ax1.scatter(dis_voltages, tavg_disvolt_VelocityArray, label='Time-Averaged Velocity (110mT Nozzle)', color='brown', s=50)
#ax1.scatter(dis_voltages, arrival_disvolt_VelocityArray, label='Plume Arrival Velocity (110mT Nozzle)', color='blue', s=50)
#ax1.set_ylabel('Velocity (km/s)')
#ax1.set_xlabel('Discharge Voltage (kV)')
#ax1.set_ylim(25, 350)
#ax1.set_xticks(xticks1)
#ax1.legend(loc='best', frameon=False)
#for i,j in zip(dis_voltages, tavg_disvolt_VelocityArray):
#    ax1.text(i,j+7, '{} km/s'.format(j))
#for i,j in zip(dis_voltages,arrival_disvolt_VelocityArray):
#    ax1.text(i,j+7, '{} km/s'.format(j))

ax1.scatter(nozzle_fields2, tavg_nozzle_VelocityArray, label='Time-Averaged Velocity (3.5kV)', color='red', s=50)
#ax2.scatter(nozzle_fields1, arrival_nozzle_VelocityArray, label='Plume Arrival Velocity (3.5kV)', color='green', s=50)
ax1.set_xlabel('Magnetic Nozzle Field (mT)')
ax1.set_ylabel('Velocity (km/s)')
ax1.set_xticks(xticks2)
ax1.set_ylim(60,85)
ax1.legend(loc='upper left', frameon=False)
for i,j in zip(nozzle_fields2, tavg_nozzle_VelocityArray):
    ax1.text(i+1,j+0.3, '({} mT, {} km/s)'.format(i,j), fontsize='small')
#for i,j in zip(nozzle_fields1, arrival_nozzle_VelocityArray):
#    ax2.text(i+1,j-40, '{} km/s'.format(j))
"""
#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\MagStructures_BulkVelocity\\ArrivalVsTimeAvgdVelocities_NozzleVsDisV.png', dpi=300)



