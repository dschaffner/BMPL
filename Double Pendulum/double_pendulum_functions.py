# double_pendulum_functions.py

import numpy as np
from scipy.integrate import odeint


def deriv(y, t, L1, L2, m1, m2):
    """Return the first derivatives of y = theta1, z1, theta2, z2."""
    theta1, z1, theta2, z2 = y

    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1dot = z1
    z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    theta2dot = z2
    z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) +
             m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    return theta1dot, z1dot, theta2dot, z2dot


def calc_E(y):
    """Return the total energy of the system."""

    th1, th1d, th2, th2d = y.T
    V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
    T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
                                      2*L1*L2*th1d*th2d*np.cos(th1-th2))
    return T + V


def double_pendulum_calc(tmax, dt, y0, p1, p2, g, edrift=0.05, savefile=False):
    t = np.arange(0, tmax+dt, dt)
    L1 = p1[0]
    L2 = p2[0]
    m1 = p1[1]
    m2 = p2[1]

    # Compute timeseries of pendulum by solving the differential equations
    y = odeint(deriv, y0, t, args=(L1, L2, m1, m2))

    # Check the accumulated error
    cumerr = np.cumsum(np.abs(calc_E(y)-calc_E(y0)))
    if np.any(cumerr > edrift) == True:
        passed_threshold = np.where(cumerr > edrift)[0][0]
        print('Accumulated Error Exceeds Threshold at timestep: ' +
              str(passed_threshold))

    # angle [rad] pendulum 1 makes with vertical as function of time
    theta1 = y[:, 0]
    # angle [rad] pendulum 2 makes with vertical as function of time
    theta2 = y[:, 2]
    thetadot1 = y[:, 1]  # angluar velocity [rad/s] of pendulum 1
    thetadot2 = y[:, 3]  # angluar velocity [rad/s] of pendulum 2

    if savefile == True:
        datadir = 'savefiles/'
        filename = 'dp_thetas.npz'
        np.savez(datadir+filename, theta1=theta1, theta2=theta2,
                 thetadot1=thetadot1, thetadot2=thetadot2,
                 m1=m1, m2=m2, L1=L1, L2=L2, ic=y0, g=g,
                 time=t, dt=dt, tmax=tmax, nsteps=np.shape(t)[0],
                 E=calc_E(y), E0=calc_E(y0), cumerr=cumerr, edrift=edrift)

    return t, theta1, thetadot1, theta2, thetadot2


def convert_dp_to_cart(theta1, thetadot1, theta2, thetadot2, p1, p2):
    L1 = p1[0]
    L2 = p2[0]
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)

    return x1, y1, x2, y2


# gravity
g = 9.81  # m/s2

# Pendulum 1 (upper pendulum)
m1 = 1  # mass
L1 = 1  # length of rod
p1 = np.array([L1, m1])

# Pendulum 2 (lower pendulum)
m2 = 1  # mass
L2 = 1  # length of rod
p2 = np.array([L2, m2])

# Initial Conditions
#[theta1, dtheta1/dt, theta2, dtheta2/dt]
y0 = np.array([3*np.pi/7, 0, 3*np.pi/4, 0])

edrift = 0.05

tmax = 30  # seconds
dt = 0.01  # seconds
