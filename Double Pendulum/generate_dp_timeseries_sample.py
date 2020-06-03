#generate_dp_timeseries_sample.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import double_pendulum_functions as dp

#gravity
g=9.81#m/s2

#Pendulum 1 (upper pendulum)
m1=1#mass
L1=1#length of rod
p1=np.array([L1,m1])

#Pendulum 2 (lower pendulum)
m2=1#mass
L2=1#length of rod
p2=np.array([L2,m2])

#Initial Conditions (in rad or rad/s)
#[theta1, dtheta1/dt, theta2, dtheta2/dt]
y0 = np.array([3*np.pi/7, 0, 3*np.pi/4, 0])

edrift = 0.05

tmax = 30#seconds
dt = 0.01#seconds

time,t1,t1dot,t2,t2dot = dp.double_pendulum_calc(tmax,dt,y0,p1,p2,g,edrift=edrift,savefile=True)

x1,y1,x2,y2 = dp.convert_dp_to_cart(t1,t1dot,t2,t2dot,p1,p2)