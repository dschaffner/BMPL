#probablity distribution function based on changes in signal at 
#given time steps
# Adopted and updated from David A. Schaffner (old_delta_pdf.py)
from typing_extensions import final
import numpy as np


def delta_spatial():
    """ 
    """

    return None

def delta_time(time, dt, windowSize):
    """ Generate the PDF of increments given timestep with time coordinate.
    """
    def mean_time(time1, time2):
        mtime = (time1 + time2)/2
        return mtime
    timeHold = []
    num_step = np.round(windowSize/dt)
    total = len(time)
    initial = 0
    steps = int(num_step)
    while steps < total:
        
        #print 'Index: ',steps, ' Delta: ',delta
        timeHold.append(mean_time(time[initial], time[initial + int(num_step)]))

        initial = steps
        steps = steps+int(num_step)

    return timeHold

def delta_pdf(timeseries,dt,windowSize):
    """ Generate the PDF of increments given timestep
    """
    deltas=[]
    indexstep = np.round(windowSize/dt)
    total = len(timeseries)
    initial = 0
    steps = int(indexstep)
    while steps < total:
        
        delta = diff_algorithm(timeseries, initial, steps)
        #print 'Index: ',steps, ' Delta: ',delta
        deltas.append(delta)

        initial = steps
        steps = steps+int(indexstep)
        
    return deltas

def diff_algorithm(timeseries, initial_index, final_index):
    """
        This function finds the difference of the timeseries at initial_index and final_index.
    """

    initial=timeseries[initial_index]
    final=timeseries[final_index]

    delta = final-initial
    
    return np.asarray(delta)