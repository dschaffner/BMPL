#test_wavelet_run

import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import process_wavelet as pw
import compute_wavelet_spectrum as cw

#generate sample data
#x=np.arange(10000)*0.000001
#x=x[3500:]
#y=np.sin(2.0/x)
#plt.figure(1)
#plt.plot(x,y)
#plt.show()

#generate chirp
fs = 800
T = 10
t = np.linspace(0, T, T*fs, endpoint=False)
w = sp.signal.chirp(t, f0=0.001, f1=200, t1=10, method='logarithmic')
#w = sp.signal.chirp(t, f0=0.01, f1=100, t1=10, method='linear')
plt.plot(t,w)
plt.show

#compute wavelet
wvfreq,waveletpwr=pw.process_wavelet(w,t)

#plot wavelet
plt.figure(2)
im=plt.contourf(t,np.log(wvfreq),(waveletpwr))