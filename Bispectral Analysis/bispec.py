#bispec.py
#A function to compute the single shot bispectrum (i.e. unaveraged) as a function of frequency
#That is: output = FFT(f1)*FFT(f2)*complexconj(FFT(f1+f2 or f1-f2))

import numpy as np
import spectrum_wwind as sw
import indexfinderfuncs as iff

def bispec(timearray,array1,array2,array3,maxfreq,auto=True):
    freq,freq2,comp1,pwr,mag,phase2,cos_phase,dt = sw.spectrum_wwind(array1,timearray)
    if auto==False:
        freq,freq2,comp2,pwr,mag,phase2,cos_phase,dt = sw.spectrum_wwind(array2,timearray)
        
        freq,freq2,comp3,pwr,mag,phase2,cos_phase,dt = sw.spectrum_wwind(array3,timearray)
    if auto==True:
        comp2 = comp1
        comp3 = comp1
    
    #create list of positive and negative frequencies
    maxfreq_index = iff.tindex_min(freq2,maxfreq)
    print('For ',maxfreq,' Hz, index is ',maxfreq_index)
    posfreqs = freq2[0:maxfreq_index]#pos freqs start at element 0, take first x chunk
    negfreqs = freq2[-maxfreq_index:]#neg freqs count up from min at halfway point, thus take last x chunk
    #print(posfreqs.shape)
    #print(negfreqs.shape)
    frequencies = np.append(negfreqs,posfreqs)
    #print(frequencies.shape)
    
    #create bispectral storage arrays based on frequency array size
    bisp = np.zeros([frequencies.shape[0]//2,376],dtype=complex)    
    norm1=np.zeros([frequencies.shape[0]//2,376],dtype=complex)
    norm2=np.zeros([frequencies.shape[0]//2,376],dtype=complex)
    #halfindex=(frequencies.shape[0]//2)+1
    
    #determine scaling amount so that frequencies are intergers
    power=np.log10(freq2[1])
    if power < 0.0:
        powerceil=np.ceil(np.abs(power))
        decimalpower=10.0**(powerceil)
    else:
        decimalpower=1.0
    print('decimalpower ',decimalpower)
    scaled_frequencies = frequencies*decimalpower
    rounded_frequencies = np.around(scaled_frequencies)
    rounded_posfrequencies = np.around(posfreqs*decimalpower)
    
    for s in np.arange(1,posfreqs.shape[0]):
        #print('On frequency',s,'  ',rounded_frequencies[s])
        for t in np.arange(1,376):
            
            f1=rounded_posfrequencies[s]
            f2=rounded_frequencies[t]
            if f1==f2: 
                print(s,' and ',t)
                print(f1,' and ',f2)
                continue
            f3=f1+f2#find sum of rounded frequncies
            #f3arg=np.nonzero(rounded_frequencies==f3)#find index of frequncies array that matches the summed freq, f3
            f3arg=np.where(rounded_frequencies[251:]==f3)#find index of frequncies array that matches the summed freq, f3
            #if len(f3arg[0])==0:
                #print('f3arg is zero for s=',s,' and t=',t)
            if len(f3arg[0])!=0:#only compute bispectrum if there is a frequency in the frequencies array that matches f3
                #print(f3arg[0][0])
                #print('###################')
                #print(f1,'+',f2,'=',f3)
                #print(frequencies[s],'+',frequencies[t],'=',frequencies[f3arg[0][0]])
                #print(f3arg[0])
                #if f3>500:
                #if s==10 and t==10:
                #    print(f3)
                #    break
                bisp[s,t] = comp1[s]*comp2[t]*np.conj(comp3[f3arg[0][0]])
                norm1[s,t] = comp1[s]*comp2[t]
                norm2[s,t] = comp3[f3arg[0][0]]
            
    #print(freq2)
    return frequencies,posfreqs,bisp,norm1,norm2