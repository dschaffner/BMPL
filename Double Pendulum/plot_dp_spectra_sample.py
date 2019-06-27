#generate_dp_timeseries_sample.py

import numpy as np
import matplotlib.pyplot as plt
import loadnpzfile as ld
import spectrum_wwind as spec

datadir = 'savefiles/'
filename = 'dp_thetas.npz'
dp_thetas = ld.loadnpzfile(datadir+filename)
theta1=dp_thetas['theta1']
theta2=dp_thetas['theta2']
time=dp_thetas['time']

freq1,freq_alt1,comp1,pwr1,mag1,phase1,cos_phase1,dt1=spec.spectrum_wwind(np.cos(theta1),time,window='hanning')
freq2,freq_alt2,comp2,pwr2,mag2,phase2,cos_phase2,dt2=spec.spectrum_wwind(np.cos(theta2),time,window='hanning')

plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2,markeredgewidth=0.0)
fig=plt.figure(num=1,figsize=(6,3.5),dpi=300,facecolor='w',edgecolor='k')
left  = 0.2  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.17  # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.semilogx(freq1[1:],pwr1[1:],color='blue',label=r'FFT[$\cos(\theta_{1})(t)$]')
plt.semilogx(freq2[1:],pwr2[1:],color='red',label=r'FFT[$\cos(\theta_{2})(t)$]')

plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlabel('Freq [Hz]',fontsize=9)
plt.ylabel('Power (arb)',fontsize=9)
plt.xlim(1e-1,3e0)
plt.ylim(0,12)
plt.legend(loc='upper right',fontsize=5,frameon=False,handlelength=5)

