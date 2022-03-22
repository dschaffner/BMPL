#probablity distribution function based on changes in signal at 
#given time steps

import numpy as np

def delta_pdf(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries)
    initial = timeseries[0]
    steps = int(indexstep)
    #print 'Initial Steps = ',steps
    numloops = 0
    while steps < total:
        numloops+=1
        delta = timeseries[steps]-initial
        #print 'Index: ',steps, ' Delta: ',delta
        deltas.append(delta)
        initial = timeseries[steps]
        steps = steps+int(indexstep)
        
    return deltas
    
def delta_pdf_V2(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries)
    initial_index = 0
    final_index = int(indexstep)
    initial=timeseries[initial_index]
    final=timeseries[final_index]
    steps = int(indexstep)
    #print 'Initial Steps = ',steps
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        delta = final-initial
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial=timeseries[initial_index]
        final=timeseries[final_index]
        #print 'Index diff = ',final_index-initial_index
    return deltas

#find angle between vectors
def delta_pdf_angle(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries[0,:])
    initial_index = 0
    final_index = int(indexstep)
    initial=timeseries[:,initial_index]
    final=timeseries[:,final_index]
    steps = int(indexstep)
    #print 'Initial Steps = ',steps
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        delta = (180/np.pi)*np.arctan2(np.linalg.norm(np.cross(initial,final)),np.inner(initial,final))
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial=timeseries[:,initial_index]
        final=timeseries[:,final_index]
        #print 'Index diff = ',final_index-initial_index
    return deltas

#determine magnitude of difference vector
def delta_pdf_vecdiff(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries[0,:])
    initial_index = 0
    final_index = int(indexstep)
    initial=timeseries[:,initial_index]
    final=timeseries[:,final_index]
    steps = int(indexstep)
    #print 'Initial Steps = ',steps
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        vecdiff = (final-initial)
        delta = np.linalg.norm(vecdiff)
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial=timeseries[:,initial_index]
        final=timeseries[:,final_index]
        #print 'Index diff = ',final_index-initial_index
    return deltas
    
#determine DeltaB over B (i.e |B2-B1|/|B1|
def delta_pdf_deltavec(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries[0,:])
    initial_index = 0
    final_index = int(indexstep)
    initial=timeseries[:,initial_index]
    final=timeseries[:,final_index]
    steps = int(indexstep)
    #print 'Initial Steps = ',steps
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        vecdiff = (final-initial)
        deltaup = np.linalg.norm(vecdiff)
        deltadown = np.linalg.norm(initial)
        deltas.append(deltaup/deltadown)
        initial_index+=1
        final_index+=1
        initial=timeseries[:,initial_index]
        final=timeseries[:,final_index]
        #print 'Index diff = ',final_index-initial_index
    return deltas