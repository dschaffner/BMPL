import numpy as np
""" These functions are from Schaffner + 2014 
MULTIFRACTAL AND MONOFRACTAL SCALING IN A LABORATORY MAGNETOHYDRODYNAMIC TURBULENCE EXPERIMENT
"""

def structure_function(pdf, power):
    """ Compute the "power"th order structure function.
    """
    sFunc = np.average(np.abs(pdf)**power)
    return sFunc

def normalized_structure_function(pdf, power):
    """ Compute the "power"th order normalized structure function.
    """
    sFunc = structure_function(pdf, power)
    norm = structure_function(pdf, 2)**(power/2)
    norm_sFunc = sFunc/norm
    
    return norm_sFunc