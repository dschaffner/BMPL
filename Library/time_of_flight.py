import numpy as np
import bmag_spatial_correlation as bsc
import bmx


def get_tau_spatial(tau, shot):
    """ Get the separation array based on the inputed time delay array.
    """
    velocity = get_bulk_speed(shot)
    spatial = tau*(1e-6)*velocity*1e3*100

    return spatial


def get_spatial_tau(spatial, shot):
    """ Get the time delay array based on the inputed separation array.
    """
    velocity = get_bulk_speed(shot)
    tau = spatial*1e6/(velocity*1e3*100)

    return tau


def compute_time_delay(refPos, pos, shot, window=[75, 150], norm='ref', normalized=False):
    """ Output the time delay corresponding to the maximum correlation.
    """
    magwCorr, tau = get_shot_correlation(
        refPos, pos, shot, window, norm, normalized)
    zeroDex = bmx.finding_Index_Time(tau*1e-6, -10)
    limitDex = bmx.finding_Index_Time(tau*1e-6, 10)

    magwCorr = magwCorr[zeroDex:limitDex]
    tau = tau[zeroDex:limitDex]
    maxIndex = np.where(magwCorr == np.max(magwCorr))[0][0]
    timeDelay = tau[maxIndex]

    return timeDelay


def get_shot_correlation(refPos, pos, shot, window=[75, 150], norm='ref', normalized=False):
    """ Output the correlation and time delay of the windowed bmags form 
    specified probes.
    """
    magwCorr, tau = bsc.gen_shot_pos_unfiltered_bmag_correlation(
        refPos, pos, shot, window, norm=norm, normalized=normalized)

    return magwCorr, tau


def get_bulk_speed(shot):
    """ Outputs the tabulated bulk velocity for specified shot
    time_delay_array = np.asarray([0.817, 0.4, 0.617, 0.55, 0.6, 0.37, 0.63, 0.52, 0.65, 0.3,
                               1.34, 1.12, 0.38, 1.05, 0.6, 0.55, 0.83, 0.87, 0.76, 0.67, 0.93, 1.00, 0.65, 0.58, 0.84])
    speed_array = 2.6e-2/(time_delay_array*1e-6)
    """
    time_delay_array = np.asarray([0.817, 0.4, 0.617, 0.55, 0.6, 0.37, 0.63, 0.52, 0.65, 0.3,
                                   1.34, 1.12, 0.38, 1.05, 0.6, 0.55, 0.83, 0.87, 0.76, 0.67, 0.93, 1.00, 0.65, 0.58, 0.84])
    speed_array = 2.6e-2/(time_delay_array*1e-6)*1e-3

    return speed_array[shot]
