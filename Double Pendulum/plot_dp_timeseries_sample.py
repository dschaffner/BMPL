# generate_dp_timeseries_sample.py

import numpy as np
import matplotlib.pyplot as plt
import loadnpzfile as ld

datadir = 'savefiles/'
filename = 'dp_thetas.npz'
dp_thetas = ld.loadnpzfile(datadir+filename)
theta1 = dp_thetas['theta1']
theta2 = dp_thetas['theta2']
time = dp_thetas['time']

plt.rc('axes', linewidth=0.75)
plt.rc('xtick.major', width=0.75)
plt.rc('ytick.major', width=0.75)
plt.rc('xtick.minor', width=0.75)
plt.rc('ytick.minor', width=0.75)
plt.rc('lines', markersize=2, markeredgewidth=0.0)
fig = plt.figure(num=1, figsize=(6, 3.5), dpi=300,
                 facecolor='w', edgecolor='k')
left = 0.2  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.17  # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right,
                    top=top, wspace=wspace, hspace=hspace)
ax1 = plt.subplot(1, 1, 1)
plt.plot(time, np.cos(theta1), color='blue', label=r'$\cos(\theta_{1})(t)$')
plt.plot(time, np.cos(theta2), color='red', label=r'$\cos(\theta_{2})(t)$')

plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlabel('Time [s]', fontsize=9)
plt.ylabel(r'$\cos(\theta)$', fontsize=9)
plt.xlim(0, 30.0)
plt.ylim(-1.3, 1.1)
plt.legend(loc='lower right', fontsize=5, frameon=False, handlelength=5)
