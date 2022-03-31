# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 12:02:16 2022

@author: dschaffner
"""

import matplotlib.pylab as plt
import numpy as np
import os
import indexfinderfuncs as iff
def loadnpzfile(filename, supress=False):

    savefile = filename
    savefile = os.path.normpath(savefile)

    file = np.load(savefile)
    if not supress:
        print('Arrays loaded: ')
        for arr in file.files:
            print(arr)
    return file

directory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
processed_wavelet = loadnpzfile(directory+'all_processed_wavlet_spectra_corrvel.npz')
mean_power_interp_60t160 = processed_wavelet['mean_power_interp_60t160']
mean_power_interp_50t150 = processed_wavelet['mean_power_interp_50t150']
mean_power_interp_50t125 = processed_wavelet['mean_power_interp_50t125']
mean_power_interp_50t100 = processed_wavelet['mean_power_interp_50t100']
mean_power_interp_100t150 = processed_wavelet['mean_power_interp_100t150']
mean_power_interp_60t110 = processed_wavelet['mean_power_interp_60t110']
mean_power_interp_60t135 = processed_wavelet['mean_power_interp_60t135']
mean_power_interp_60t140 = processed_wavelet['mean_power_interp_60t140']
mean_power_interp_60t160_trace = processed_wavelet['mean_power_interp_60t160_trace']
mean_power_interp_50t150_trace = processed_wavelet['mean_power_interp_50t150_trace']
mean_power_interp_50t125_trace = processed_wavelet['mean_power_interp_50t125_trace']
mean_power_interp_50t100_trace = processed_wavelet['mean_power_interp_50t100_trace']
mean_power_interp_100t150_trace = processed_wavelet['mean_power_interp_100t150_trace']
mean_power_interp_60t110_trace = processed_wavelet['mean_power_interp_60t110_trace']
mean_power_interp_60t135_trace = processed_wavelet['mean_power_interp_60t135_trace']
mean_power_interp_60t140_trace = processed_wavelet['mean_power_interp_60t140_trace']
mean_Bwv60t160 = processed_wavelet['mean_Bwv60t160']
mean_Bwv50t150 = processed_wavelet['mean_Bwv50t150']
mean_Bwv50t125 = processed_wavelet['mean_Bwv50t125']
mean_Bwv50t100 = processed_wavelet['mean_Bwv50t100']
mean_Bwv100t150 = processed_wavelet['mean_Bwv100t150']
mean_Bwv60t110 = processed_wavelet['mean_Bwv60t110']
mean_Bwv60t135 = processed_wavelet['mean_Bwv60t135']
mean_Bwv60t140 = processed_wavelet['mean_Bwv60t140']
mean_Bwv60t160_trace = processed_wavelet['mean_Bwv60t160_trace']
mean_Bwv50t150_trace = processed_wavelet['mean_Bwv50t150_trace']
mean_Bwv50t125_trace = processed_wavelet['mean_Bwv50t125_trace']
mean_Bwv50t100_trace = processed_wavelet['mean_Bwv50t100_trace']
mean_Bwv100t150_trace = processed_wavelet['mean_Bwv100t150_trace']
mean_Bwv60t110_trace = processed_wavelet['mean_Bwv60t110_trace']
mean_Bwv60t135_trace = processed_wavelet['mean_Bwv60t135_trace']
mean_Bwv60t140_trace = processed_wavelet['mean_Bwv60t140_trace']
wavenum = processed_wavelet['wavenum']
freq = processed_wavelet['freq']

savedirectory = directory+'paperPlots\\'

plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)


### Plot 1: Frequency Spectra of 1 position
fig=plt.figure(num=1,figsize=(5,3.5),dpi=600,facecolor='w',edgecolor='k')
left  = 0.17  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.21  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.loglog(freq,mean_Bwv60t160[0,0,:],color='blue',label=r'$B_{r}$')
plt.loglog(freq,mean_Bwv60t160[0,1,:],color='green',label=r'$B_{\theta}$')
plt.loglog(freq,mean_Bwv60t160[0,2,:],color='red',label=r'$B_{z}$')
plt.loglog(freq,1e4*freq**(-5/3),color='black',linestyle='dashed',label=r'$f^{-5/3}$')
plt.xlabel(r'$f [Hz]$',fontsize=14)
plt.ylabel(r'$Arb. Power$',fontsize=14)
plt.legend(loc='lower left',fontsize=10,frameon=False,handlelength=5)
plt.savefig(savedirectory+'freq_spectra_60t160_pos5_rtz.png',dpi=600,facecolor='white',edgecolor='black')
plt.clf()
plt.close()

### Plot 2: Freq vs Wavenum spectra for three positions
labelfontsize=12
legendfontsize=8
fig=plt.figure(num=2,figsize=(10,3.5),dpi=600,facecolor='w',edgecolor='k')
left  = 0.17  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.21  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.25   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,2,1)
plt.loglog(freq,mean_Bwv60t160_trace[0,:],color='blue',label='Probe 5')
plt.loglog(freq,mean_Bwv60t160_trace[2,:],color='green',label='Probe 19')
plt.loglog(freq,mean_Bwv60t160_trace[4,:],color='red',label='Probe 33')
plt.loglog(freq,1e4*freq**(-5/3),color='black',linestyle='dashed',label=r'$f^{-5/3}$')
plt.xlabel(r'$f\,$ [Hz]',fontsize=labelfontsize)
plt.ylabel(r'$Arb. Power$',fontsize=labelfontsize)
plt.legend(loc='lower left',fontsize=legendfontsize,frameon=False,handlelength=5)

ax2=plt.subplot(1,2,2)
plt.loglog(wavenum,mean_power_interp_60t160_trace[0,:],color='blue',label='Probe 5')
plt.loglog(wavenum,mean_power_interp_60t160_trace[2,:],color='green',label='Probe 19')
plt.loglog(wavenum,mean_power_interp_60t160_trace[4,:],color='red',label='Probe 33')
plt.loglog(wavenum,1e-4*wavenum**(-5/3),color='black',linestyle='dashed',label=r'$\left(\frac{k}{2\pi}\right)^{-5/3}$')
plt.xlabel(r'$\left(\frac{k_{v}}{2\pi}\right)$ [m$^{-1}$]',fontsize=labelfontsize)
plt.ylabel(r'$Arb. Power$',fontsize=labelfontsize)
plt.legend(loc='lower left',fontsize=legendfontsize,frameon=False,handlelength=5)

#plt.savefig(savedirectory+'freq_and_wavenum_spectra_60t160_pos5-19-33_trace_initialplume_v.png',dpi=600,facecolor='white',edgecolor='black')
plt.savefig(savedirectory+'freq_and_wavenum_spectra_60t160_pos5-19-33_trace_corr_v.png',dpi=600,facecolor='white',edgecolor='black')
plt.clf()
plt.close()


### Plot 3: Wavenumber spectra fits
labelfontsize=12
legendfontsize=6
fig=plt.figure(num=1,figsize=(4,5),dpi=300,facecolor='w',edgecolor='k')
left  = 0.2  # the left side of the subplots of the figure
right = 0.94    # the right side of the subplots of the figure
bottom = 0.23  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.loglog(wavenum,mean_power_interp_60t160[0,0,:],color='blue',label='Probe 5')
plt.loglog(wavenum,mean_power_interp_60t160_trace[2,:],color='green',label='Probe 19')
plt.loglog(wavenum,mean_power_interp_60t160_trace[4,:],color='red',label='Probe 33')
plt.xlabel(r'$\left(\frac{k_{v}}{2\pi}\right)$ [m$^{-1}$]',fontsize=labelfontsize)
plt.ylabel(r'$Arb. Power$',fontsize=labelfontsize)
plt.legend(loc='lower left',fontsize=legendfontsize,frameon=False,handlelength=5)


#range 1
import MLE as mle
wavenum_range = [2.7e0,7e0]
wvtindex1 = iff.tindex_min(wavenum,wavenum_range[0])
wvtindex2 = iff.tindex_min(wavenum,wavenum_range[1])
#position 5
max_alpha,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[0,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha))
scale_dat = mean_power_interp_60t160_trace[0,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='blue',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(7e0,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='blue')
print(max_alpha)


#position 19
max_alpha2,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[2,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha2))
scale_dat = mean_power_interp_60t160_trace[2,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='green',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(7e0,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha2,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='green')
print(max_alpha2)

#position 33
max_alpha3,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[4,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha3))
scale_dat = mean_power_interp_60t160_trace[4,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],0.5*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='red',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(7e0,0.5*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha3,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='red')

print(max_alpha3)


#range 2
import MLE as mle
wavenum_range = [8e0,4e1]
wvtindex1 = iff.tindex_min(wavenum,wavenum_range[0])
wvtindex2 = iff.tindex_min(wavenum,wavenum_range[1])
#position 5
max_alpha,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[0,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha))
scale_dat = mean_power_interp_60t160_trace[0,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='blue',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(3e1,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='blue')
print(max_alpha)
#position 19
max_alpha2,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[2,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha2))
scale_dat = mean_power_interp_60t160_trace[2,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='green',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(3e1,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha2,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='green')
print(max_alpha2)
#position 33
max_alpha3,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[4,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha3))
scale_dat = mean_power_interp_60t160_trace[4,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],0.5*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='red',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(3e1,0.5*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha3,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='red')
print(max_alpha3)

#range 3
import MLE as mle
wavenum_range = [5e1,1e2]
wvtindex1 = iff.tindex_min(wavenum,wavenum_range[0])
wvtindex2 = iff.tindex_min(wavenum,wavenum_range[1])
#position 5
max_alpha,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[0,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha))
scale_dat = mean_power_interp_60t160_trace[0,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='blue',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(2e2,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='blue')
print(max_alpha)
#position 19
max_alpha2,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[2,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha2))
scale_dat = mean_power_interp_60t160_trace[2,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],2.0*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='green',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(2e2,2.0*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha2,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='green')
print(max_alpha2)
#position 33
max_alpha3,alpha,L_vals,cdf_dat,cdf_mod,KS,sigma = mle.mle_werr(mean_power_interp_60t160_trace[4,:],wavenum,wvtindex1,wvtindex2)
func = (wavenum**(-max_alpha3))
scale_dat = mean_power_interp_60t160_trace[4,wvtindex1]
scale_fit = func[wvtindex1]
ratio = scale_dat/scale_fit
plt.loglog(wavenum[wvtindex1:wvtindex2],0.5*ratio*func[wvtindex1:wvtindex2],linestyle='dotted',color='red',linewidth=0.5)#,label='Fixed Range Fit: -'+str(max_alpha) )
sigma = np.round(sigma,decimals=2)
plt.text(2e2,0.5*ratio*func[wvtindex2],r'$\alpha = $'+str(np.round(max_alpha3,decimals=2))+r'$\pm$'+str(sigma),fontsize=4,color='red')
print(max_alpha3)
