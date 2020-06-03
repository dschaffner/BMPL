#array_to_ensemble.py
#This is a function that will take a large timeseries array and break it into smaller sections
import numpy as np

def array_to_ensemble(array,num_sections):
    arrelements=array.shape[0]
    sectionsize = arrelements//num_sections
    array_out = np.zeros([sectionsize,num_sections])
    for sec in np.arange(num_sections):
        array_out[:,sec]=array[sec*sectionsize:(sec+1)*sectionsize]
    return array_out