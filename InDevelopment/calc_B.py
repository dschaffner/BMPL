# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 12:55:29 2024

@author: dschaffner
"""

import astropy.constants as con
import numpy as np
import matplotlib.pylab as plt

def B_wire_G(I,r):
    mu=con.mu0
    B = (mu.value)*I/(2*np.pi*r)
    B=B*1e4#convert to gauss
    return B

I=2500.0
r=(np.arange(500)*0.001)+0.001#meters

Bwire=B_wire_G(I,r)

plt.plot(r,Bwire)

    