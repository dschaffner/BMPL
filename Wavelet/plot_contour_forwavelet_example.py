plotscript = 'plot_contour_Bwaveletspec_forPPFCpaper.py'

import loadfile as lf
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import ssx_functions as ssxf
import MLE as mle
import numpy as np
import os
plt.ioff()#Turn interactive mode off so that plots don't appear automatically without plt.show()

day = '081413'
shotrange = 'Shots36to36'
#Bfile=lf.loadfile('wavelet_spec\B_chan1_wvspec_fixtime_'+day+'_'+shotrange+'.npz')
Bfile=lf.loadfile('wavelet_spec\B_chan1_wvspec_fulltime_'+day+shotrange+'.npz')
Bpwr_count = Bfile['Bpwr_tot_count']
Bwvfreq = Bfile['wvfreq']
time=Bfile['time']

axis = 1
if axis == 0: axisname = 'r'
if axis == 1: axisname = 't'
if axis == 2: axisname = 'z'

Bpwr = (Bfile['Bpwr_full'][axis,:,:])/Bpwr_count

plt.rc('axes',linewidth=1.75)
plt.rc('xtick.major',width=1.75)
plt.rc('ytick.major',width=1.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)

fig=plt.figure(num=1,figsize=(6,3.5),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax=plt.axes([0.15,0.15,0.80,0.65])

time_lims = [0,120]
pwr_lims = [-14,0]
min_freq = 2e4

#ax.axhline(linewidth=4)
#ax.axvline(color='blue')
Bpwr = Bpwr/Bpwr.max() #normalize to max
plotcwt = Bpwr[0:-1]
logplotcwt = np.log10(np.flipud(plotcwt))
#plotcwt = np.flipud(plotcwt)

im=plt.imshow(logplotcwt,extent=[time_lims[0],time_lims[-1],Bwvfreq[0],Bwvfreq[-1]],
        vmin=pwr_lims[0],vmax=pwr_lims[1],aspect='auto')
#im=plt.imshow(plotcwt,extent=[time[0],time[-1],Bwvfreq[0],Bwvfreq[-1]],
    #norm=LogNorm(pwr_lims[0],pwr_lims[1]),aspect='auto')
#plt.pcolor(time,Bwvfreq,plotcwt,norm=LogNorm(pwr_lims[0],pwr_lims[1]),cmap='gist_rainbow')

#modify axis settings
ax.set_yscale('log')
plt.yticks(fontsize=16)
plt.xticks(fontsize=9)
plt.ylim(min_freq,Bwvfreq[0])

plt.xticks(np.arange(0,130,10),fontsize=14)
plt.xlabel(r't [$\mu$s]',fontsize=14)
plt.ylabel(r'$f$ [Hz]',fontsize=14)

plt.vlines(9,2e4,1e8,color='white',linestyles='dotted',linewidth=1.0)
plt.vlines(40,2e4,1e8,color='white',linestyles='dotted',linewidth=1.0)
plt.vlines(60,2e4,1e8,color='white',linestyles='dotted',linewidth=1.0)
plt.text(11,3e4,'Formation/\nSelective Decay '+r'$\rightarrow$',fontsize=8,horizontalalignment='left',color='white')
plt.text(50,3e4,r'$\leftarrow$ Equlib. $\rightarrow$',fontsize=8,horizontalalignment='center',color='white')
plt.text(61,3e4,r'$\leftarrow$ Dissipation',fontsize=8,horizontalalignment='left',color='white')


ax2=plt.axes([0.15,0.89,0.80,0.05])
cbar = plt.colorbar(im,cax=ax2,ticks=[0,-2,-4,-6,-8,-10,-12,-14],orientation='horizontal')
ax2.set_xticklabels([r'$1$',r'$10^{-2}$',r'$10^{-4}$',r'$10^{-6}$',r'$10^{-8}$',r'$10^{-10}$',r'$10^{-12}$',r'$10^{-14}$'],fontname='Bitstream Vera Sans',fontsize=5)
plt.xticks(fontsize=14)
plt.text(0.7,0.95,'Log Power (Norm.)',fontsize=10,transform=fig.transFigure)

#label plot
#plt.text(0.01,0.01,plotscript,fontsize=6,transform=fig.transFigure,horizontalalignment='left')


process_dir = 'C:/Users/David Schaffner/Documents/ssxpython/plots/PPFC_paper/'
#filename = 'Bdot_wvcont_Run'+day+shotrange+'.png'
filename = 'Bdot_wvcont_Run'+day+shotrange+'.eps'
savefile = os.path.normpath(process_dir+filename)
#save figure with facecolor=white, edgecolor = black
plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
plt.clf()
plt.close(fig)
