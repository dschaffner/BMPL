#list of useful functions for analysis

from numpy import where, abs, min

def tindex(timearr,timevalue,delta=2e-2):
    tind = where(abs((timearr)-(timevalue))<delta)
    tind = tind[0][0]
    return tind
    
def tindex_min(timearr,timevalue):
    minval = min(abs((timearr)-timevalue))
    tind = where(abs((timearr)-(timevalue)) == minval)
    tind = tind[0][0]
    return tind
    
def tindex_near(timearr,timevalue,threshold):
    tinds = where(abs(timearr-timevalue) < threshold)
    return tinds
    
def firstzero(timearr,timevalue):
    for n in range(len(timearr)):
        tind = n
        if timearr[n] < 0: break
    return tind
    
