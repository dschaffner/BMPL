#calc wavelet spectra

import compute_wavelet_spectrum as cw

def process_wavelet(arr,time):
    waveletpwr,wavtot,wvfreq,fft,fftfreq = cw.compute_wavelet(arr,time)
    return wvfreq,waveletpwr