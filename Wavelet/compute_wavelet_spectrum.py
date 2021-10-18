#Compute Wavelet spectrum 

import numpy as np
import Wavelets as wv
import sys
import os
import matplotlib.pylab as plt

#custom funcs
import spectrum_wwind
#
def compute_wavelet(array,time,mother='Morlet',maxscale=2,notes=256,order=6,scaling='log',Bfield=False):
    if mother == 'Morlet': wavelet=wv.Morlet
    if mother == 'Mexican': wavelet=wv.MexicanHat
    if mother == 'Paul': wavelet=wv.Paul4
    
    #timestep
    dt = time[1]-time[0]#assumed to be in seconds
    
    #Wavelet transform the data
    cw = wavelet(array,maxscale,notes,order=order,scaling=scaling)
    cwt = cw.getdata()
    
    #Extract power spectrum
    pwr = cw.getpower()
    scalespec=np.sum(pwr,axis=1)#sum total of power over time
    
    #Calculate fourier scales
    scales = cw.getscales()
    y=cw.fourierwl*scales
    wvfreq = 1.0/(y*dt)
    
    if Bfield:
        #Calculate B-field spectrum from Bdot spectrum
        Bscalespec = scalespec/(wvfreq**2)
        Bpwr = np.zeros([pwr.shape[0],pwr.shape[1]])
        for n in range(pwr.shape[1]):
            Bpwr[:,n] = pwr[:,n]/(wvfreq**2)
        pwr = Bpwr
        scalespec = Bscalespec
    #Calculate B-field spectrum from Bdot spectrum for FFT
    fftfreq,freq2,comp,fft_pw,fft_mag,phase,cos_phase,dt = spectrum_wwind.spectrum_wwind(array,time,window='hanning')
    if Bfield:
        Bfft = fft_pw/(fftfreq**2)
        fft_pw = Bfft
    return pwr,scalespec,wvfreq,fft_pw,fftfreq

def plot_wavelet(array,shotdata,Bpwr,Bscalespec,wvfreq,Bfft,fftfreq,time,timerange=[0.0,100.0],pwr_lims=[-40,-10],min_freq=500e3,showPlot=True,savePlot=True):
    plt.ioff()#Turn interactive mode off so that plots don't appear automatically without plt.show()
    
    #Create figure
    #facecolor = white, edgecolor = black
    fig=plt.figure(num=1,figsize=(9.25,4.6),dpi=200,facecolor='w',edgecolor='k')
    
    #Create axes for 2D plot
    #ax=plt.axes([fromleft,frombottom,width,height])
    ax=plt.axes([0.3,0.105,0.68,0.55])
    plt.xlabel(r't [$\mu$s]',fontsize=12)
    plt.ylabel(r'$f$ [Hz]',fontsize=12)
    
    #prepare 2D array for plotting (flip and take log)
    plotcwt=Bpwr[0:-1]
    logplotcwt = np.log10(np.flipud(plotcwt))
    
    #create 2D image

    if not pwr_lims:
        print 'Min max range'
        vmin=logplotcwt.min()
        vmax=logplotcwt.max()
        print vmin
        print vmax
        im=plt.imshow(logplotcwt,extent=[time[0],time[-1],wvfreq[0],wvfreq[-1]],
        vmin=vmin,vmax=vmax,aspect='auto')
    if pwr_lims:
        print 'user range'
        im=plt.imshow(logplotcwt,extent=[time[0],time[-1],wvfreq[0],wvfreq[-1]],
        vmin=pwr_lims[0],vmax=pwr_lims[1],aspect='auto')
    
    #modify axis settings
    ax.set_yscale('log')
    plt.yticks(fontsize=6)
    plt.ylim(min_freq,wvfreq[0])
    
    #Create axes for time series plot (with x axis shared with 2Dplot)
    ax2=plt.axes([0.3,0.7,0.68,0.25],sharex=ax)
    
    #create time series image
    plt.plot(time,array,linewidth=0.5,color='blue')
    
    #modify axis settings
    plt.ylabel('B-dot',fontsize=12)
    plt.yticks(fontsize=6)
    plt.xticks(np.arange(0,110,10),fontsize=6)
    plt.xlim(timerange[0],timerange[1])
    
    #create total power spectrum comparison axes
    ax3=plt.axes([0.03,0.105,0.2,0.75])
    
    #plot total power spectra
    plt.loglog(fftfreq,Bfft,color='orange',linewidth=0.5,label='FFT')
    plt.loglog(wvfreq,Bscalespec,color='red',linewidth=1,label='Wavelet')
    
    #modify axis settings
    plt.xlabel(r'$f$ [Hz]',fontsize=12)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.xlim(fftfreq[0],5.0*fftfreq[-1])
    
    #make simple legend
    plt.legend(loc='lower left',fontsize=7)
    
    #label plot
    day=shotdata[0]
    nn=shotdata[1]#shot number
    axis_label = shotdata[2]#axis (r,t,z)
    channel = shotdata[3]#probe channel number
    plt.text(0.15,0.92,"Shot "+str(nn)+axis_label+" Chan:"+str(channel),
        fontsize=10,bbox=dict(facecolor='green',alpha=0.2),
        transform = fig.transFigure,horizontalalignment='center')
    plt.text(0.30,0.01,"Data Date: "+day,fontsize=6,
        transform = fig.transFigure,horizontalalignment='center')
    
    
    if showPlot and not savePlot: plt.show()
    if savePlot:
        process_dir = 'C:/Users/David Schaffner/Documents/ssxpython/plots/WaveletOutputDatabase/'#run073013_1mwb_single/chan1/'
        filename = 'wavelet_'+day+'_shot'+str(nn)+'_B'+axis_label[5]+'_chan'+str(channel)+'.png'
        savefile = os.path.normpath(process_dir+filename)
        #save figure with facecolor=white, edgecolor = black
        plt.savefig(savefile,dpi=150,facecolor='w',edgecolor='k')
        plt.clf()
        plt.close(fig)
        
def wv_quickplot(Bpwr,time,wvfreq,timerange=[0.0,100.0],pwr_lims=[-40,-10],min_freq=500e3,showPlot=True,savePlot=True):
    plt.ioff()#Turn interactive mode off so that plots don't appear automatically without plt.show()
    
    #Create figure
    #facecolor = white, edgecolor = black
    fig=plt.figure(num=1,figsize=(9.25,4.6),dpi=200,facecolor='w',edgecolor='k')
    
    #Create axes for 2D plot
    #ax=plt.axes([fromleft,frombottom,width,height])
    ax=plt.axes([0.1,0.105,0.85,0.8])
    plt.xticks(np.arange(0,130,10),fontsize=8)
    plt.xlabel(r't [$\mu$s]',fontsize=12)
    plt.ylabel(r'$f$ [Hz]',fontsize=12)
    
    #prepare 2D array for plotting (flip and take log)
    plotcwt=Bpwr[0:-1]
    logplotcwt = np.log10(np.flipud(plotcwt))
    
    #create 2D image

    if not pwr_lims:
        print 'Min max range'
        vmin=logplotcwt.min()
        vmax=logplotcwt.max()
        print vmin
        print vmax
        im=plt.imshow(logplotcwt,extent=[time[0],time[-1],wvfreq[0],wvfreq[-1]],
        vmin=vmin,vmax=vmax,aspect='auto')
    if pwr_lims:
        print 'user range'
        im=plt.imshow(logplotcwt,extent=[time[0],time[-1],wvfreq[0],wvfreq[-1]],
        vmin=pwr_lims[0],vmax=pwr_lims[1],aspect='auto')
    
    #modify axis settings
    ax.set_yscale('log')
    plt.yticks(fontsize=6)
    plt.ylim(min_freq,wvfreq[0])
    
    if showPlot and not savePlot: plt.show()
    if savePlot:
        process_dir = 'C:/Users/David Schaffner/Documents/ssxpython/plots/WaveletOutputDatabase/'#run073013_1mwb_single/chan1/'
        filename = 'wavelet_test.png'
        savefile = os.path.normpath(process_dir+filename)
        #save figure with facecolor=white, edgecolor = black
        plt.savefig(savefile,dpi=150,facecolor='w',edgecolor='k')
        plt.clf()
        plt.close(fig)

def compute_spwavelet(array,spacing,mother='Morlet',maxscale=1,notes=256,order=6,scaling='log'):
    if mother == 'Morlet': wavelet=wv.Morlet
    
    #spacing
    dr = spacing
    probe = np.arange(16)*dr
    
    #Wavelet transform the data
    cw = wavelet(array,maxscale,notes,order=order,scaling=scaling)
    cwt = cw.getdata()
    
    #Extract power spectrum
    pwr = cw.getpower()
    Bscalespec=np.sum(pwr,axis=1)#sum total of power over time
    
    #Calculate fourier scales
    scales = cw.getscales()
    y=cw.fourierwl*scales
    wvk = 1.0/(y*dr)
    
    #Calculate B-field spectrum from Bdot spectrum
    #Bscalespec = scalespec/(wvfreq**2)
    #Bpwr = np.zeros([pwr.shape[0],pwr.shape[1]])
    #for n in range(pwr.shape[1]):
    #    Bpwr[:,n] = pwr[:,n]/(wvfreq**2)
        
    #Calculate spatial spectrum from Bdot spectrum for FFT
    #fftk,fft_pw,fft_mag = spectrum.spectrum(array,probe)
    #Bfft = fft_pw#/(fftfreq**2)
    
    return Bscalespec,wvk#,Bfft,fftk