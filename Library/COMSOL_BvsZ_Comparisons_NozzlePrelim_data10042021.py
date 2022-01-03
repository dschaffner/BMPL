# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:01:40 2021

@author: Josh0
"""

# Various B vs Z plots for different coil positions.  

import numpy as np
import matplotlib.pylab as plt

plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_presentationplots.py')
#plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_tallplots.py')


comsol_location = 'C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\COMSOL_SimulationData\\Data\\10042021\\'
comsol_inner1 = 'BvsZ_btwn_outercoils_innercoil_2p935cm'
comsol_inner2 = 'BvsZ_btwn_outercoils_innercoil_9p385cm'
comsol_inner3 = 'BvsZ_btwn_outercoils_innercoil_15p835cm'
comsol_inner4 = 'BvsZ_btwn_outercoils_innercoil_22p285cm'
comsol_inner5 = 'BvsZ_btwn_outercoils_innercoil_28p735cm'
comsol_inner6 = 'BvsZ_btwn_outercoils_innercoil_32p935cm_currentPhysicalLocation'

comsol_configuration1 = 'BvsZ_OuterCoils_Closest_to_DiagSection'
comsol_configuration1_200V = '200V_BvsZ_OuterCoils_Closest_to_DiagSection_NoInnerCoil'
comsol_configuration1_300V = '300V_BvsZ_OuterCoils_Closest_to_DiagSection_NoInnerCoil'
comsol_configuration1_350V = '350V_BvsZ_OuterCoils_Closest_to_DiagSection_NoInnerCoil'
comsol_configuration2 = 'BvsZ_OuterCoil(single)_Close_to_DiagSection'
comsol_configuration2_225turns = 'BvsZ_OuterCoil(single_225turns)_Close_to_DiagSection'
comsol_configuration2_225turns_long = 'BvsZ_OuterCoil(single_225turns_DoubleLength)_Close_to_DiagSection'
comsol_configuration3_350turns = 'BvsZ_OuterCoil(single_350turns)_Close_to_DiagSection'


"""
Inner Coil Locations, centimeters into the inner electrode.  **This is with all coils firing**.
Simulation to see B vs. Z changes in electrode gap as inner coil location changes but outer coils stay the same.  
    1. 2.935cm
    2. 9.385cm
    3. 15.835cm
    4. 22.285cm
    5. 28.735cm
    6. 32.935cm - Current Physical location
    
Configuration 1: Outer Coils (300V) pushed to closest approach with each other to the Diagnostic section.  Inner coil sits at Gas Puff Valves.    
"""

com_BvZ_inner1 = np.loadtxt(comsol_location+comsol_inner1+'.txt', skiprows=8, unpack=True)
com_BvZ_inner2 = np.loadtxt(comsol_location+comsol_inner2+'.txt', skiprows=8, unpack=True)
com_BvZ_inner3 = np.loadtxt(comsol_location+comsol_inner3+'.txt', skiprows=8, unpack=True)
com_BvZ_inner4 = np.loadtxt(comsol_location+comsol_inner4+'.txt', skiprows=8, unpack=True)
com_BvZ_inner5 = np.loadtxt(comsol_location+comsol_inner5+'.txt', skiprows=8, unpack=True)
com_BvZ_inner6 = np.loadtxt(comsol_location+comsol_inner6+'.txt', skiprows=8, unpack=True)
com_Config1 = np.loadtxt(comsol_location+comsol_configuration1+'.txt', skiprows=8, unpack=True)
com_Config1_200V = np.loadtxt(comsol_location+comsol_configuration1_200V+'.txt', skiprows=8, unpack=True)
com_Config1_300V = np.loadtxt(comsol_location+comsol_configuration1_300V+'.txt', skiprows=8, unpack=True)
com_Config1_350V = np.loadtxt(comsol_location+comsol_configuration1_350V+'.txt', skiprows=8, unpack=True)
com_Config2 = np.loadtxt(comsol_location+comsol_configuration2+'.txt', skiprows=8, unpack=True)
com_Config2_225turns = np.loadtxt(comsol_location+comsol_configuration2_225turns+'.txt', skiprows=8, unpack=True)
com_Config2_225turns_long = np.loadtxt(comsol_location+comsol_configuration2_225turns_long+'.txt', skiprows=8, unpack=True)
com_Config3_350turns = np.loadtxt(comsol_location+comsol_configuration3_350turns+'.txt', skiprows=8, unpack=True)

#The position column from COMSOL plots
comsol_Z_1 = com_BvZ_inner1[0]*1e2
comsol_Z_2 = com_BvZ_inner2[0]*1e2
comsol_Z_3 = com_BvZ_inner3[0]*1e2
comsol_Z_4 = com_BvZ_inner4[0]*1e2
comsol_Z_5 = com_BvZ_inner5[0]*1e2
comsol_Z_6 = com_BvZ_inner6[0]*1e2

#The B column data from COMSOL plots
com_BvZ_1 = com_BvZ_inner1[1]*1e3
com_BvZ_2 = com_BvZ_inner2[1]*1e3
com_BvZ_3 = com_BvZ_inner3[1]*1e3
com_BvZ_4 = com_BvZ_inner4[1]*1e3
com_BvZ_5 = com_BvZ_inner5[1]*1e3
com_BvZ_6 = com_BvZ_inner6[1]*1e3

# Below are possible physical orientations for the Outer coils (without Inner coil firing) at various voltages.  The preliminary Nozzle simulations.
comsol_Z_config1 = com_Config1[0]*1e2
comsol_Z_config1_200V = com_Config1_200V[0]
comsol_Z_config1_300V = com_Config1_300V[0]
comsol_Z_config1_350V = com_Config1_350V[0]
comsol_Z_config2 = com_Config2[0]*1e2
comsol_Z_config2_225 = com_Config2_225turns[0]*1e2
comsol_Z_config2_225_long = com_Config2_225turns_long[0]*1e2
comsol_Z_config3_350 = com_Config3_350turns[0]*1e2

com_BvZ_config1 = com_Config1[1]*1e3
com_BvZ_config1_200V = com_Config1_200V[1]*1e3
com_BvZ_config1_300V = com_Config1_300V[1]*1e3
com_BvZ_config1_350V = com_Config1_350V[1]*1e3
com_BvZ_config2 = com_Config2[1]*1e3
com_BvZ_config2_225 = com_Config2_225turns[1]*1e3
com_BvZ_config2_225_long = com_Config2_225turns_long[1]*1e3
com_BvZ_config3_350 = com_Config3_350turns[1]*1e3

fig1, (ax1) = plt.subplots(1)

#These plots show field from Inner coil as its position is shifted within the inner electrode in Z (ALL COILS FIRING)
ax1.plot(comsol_Z_1, com_BvZ_1, color='red', label='Pos1')
ax1.plot(comsol_Z_2, com_BvZ_2, color='blue', label='Pos2')
ax1.plot(comsol_Z_3, com_BvZ_3, color='green', label='Pos3')
ax1.plot(comsol_Z_4, com_BvZ_4, color='black', label='Pos4')
ax1.plot(comsol_Z_5, com_BvZ_5, color='orange', label='Pos5')
ax1.plot(comsol_Z_6, com_BvZ_6, color='purple', label='Pos6 - Current')


# These plots below are for possible configurations of the Inner Coil with other outer coils, not built presently.
#ax1.plot(comsol_Z_config1, com_BvZ_config1, color='teal', label='Config1 - 300V OuterCoil Discharge')
#ax1.plot(comsol_Z_config2, com_BvZ_config2, color='brown', label='Config2 - 300V Single Outer Coil Discharge')
#ax1.plot(comsol_Z_config2_225, com_BvZ_config2_225, color='grey', label='Config2 - 300V Single (225Turns) Outer Coil')
#ax1.plot(comsol_Z_config2_225_long, com_BvZ_config2_225_long, color='pink', label='Config2 - 300V Single (225Turns, Longer), Outer Coil')
#ax1.plot(comsol_Z_config3_350, com_BvZ_config3_350, color='limegreen', label='Config3 - 300V Single (350Turns), Outer Coil')
#ax1.plot(comsol_Z_config1_200V, com_BvZ_config1_200V, color='black', label='200V, 170A')
#ax1.plot(comsol_Z_config1_300V, com_BvZ_config1_300V, color='red', label='300V, 260A')
#ax1.plot(comsol_Z_config1_350V, com_BvZ_config1_350V, color='blue', label='350V, 300A')
ax1.set_ylabel(r'B\_z (mT)')
ax1.set_xlabel(r'z (m)')
#ax1.set_ylim(-0.1, 200, 50)
#ax1.set_xlim(-0.01, 0.79)
ax1.legend(loc='best', frameon=False)

#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\COMSOL_SimulationData\\Config1_Comparison_PinchFields.pdf', dpi=600)
