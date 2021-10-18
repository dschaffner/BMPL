#calc wavelet spectra

import compute_wavelet as cw

def compute_mag_wavelet(arr,time):
    waveletpwr,wavtot,wvfreq,fft,fftfreq = cw.compute_wavelet(arr,time)
    return wvfreq,waveletpwr