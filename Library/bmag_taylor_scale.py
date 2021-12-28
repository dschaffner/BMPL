import tools
import numpy as np
import bmag_spatial_correlation as bsc
import bmag_temporal_correlation as btc
import richardson_extrapolation as re
import bmx
import time_of_flight as tof


def single_shot_bmag_biased_taylor_scales(refPos, shot, filt, window, norm='ref'):
    """ Outputs the biased Taylor scales from bmag correlations of shot referenced 
    to refPos windowed and filt."""
    wfMagCorr = bsc.gen_shot_bmag_spatial_correlation(
        refPos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tools.gen_separation_array(refPos=refPos)
    biasedTSArray, maxSeparationArray, uncertainty = biased_taylor_scales(
        wfMagCorr, separations)
    return biasedTSArray, maxSeparationArray, uncertainty


def single_shot_abs_bmag_biased_taylor_scales(refPos, shot, filt, window, norm='ref'):
    """ Outputs the biased Taylor scales from ABS-bmag correlations of shot referenced 
    to refPos windowed and filt."""
    wfMagCorr = bsc.gen_shot_bmag_spatial_correlation(
        refPos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tools.gen_separation_array(refPos=refPos)
    biasedTSArray, maxSeparationArray, uncertainty = biased_taylor_scales(
        np.abs(wfMagCorr), separations)
    return biasedTSArray, maxSeparationArray, uncertainty


def single_shot_bmag_biased_temporal_taylor_scales(refPos, shot, filt, window, norm='ref', window2=[0, 13]):
    wfMagCorr, tau = btc.gen_shot_pos_bmag_correlation(
        refPos=refPos, pos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tof.get_tau_spatial(tau, shot)
    startDex = bmx.finding_Index_Time(separations*1e-6, window2[0])
    endDex = bmx.finding_Index_Time(separations*1e-6, window2[1])
    biasedTSArray, maxSeparationArray, uncertainty = biased_taylor_scales(
        wfMagCorr[startDex:endDex], separations[startDex:endDex])
    return biasedTSArray, maxSeparationArray, uncertainty


def single_shot_abs_bmag_biased_temporal_taylor_scales(refPos, shot, filt, window, norm='ref', window2=[0, 13]):
    wfMagCorr, tau = btc.gen_shot_pos_bmag_correlation(
        refPos=refPos, pos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tof.get_tau_spatial(tau, shot)
    startDex = bmx.finding_Index_Time(separations*1e-6, window2[0])
    endDex = bmx.finding_Index_Time(separations*1e-6, window2[1])
    biasedTSArray, maxSeparationArray, uncertainty = biased_taylor_scales(
        np.abs(wfMagCorr[startDex:endDex]), separations[startDex:endDex])
    return biasedTSArray, maxSeparationArray, uncertainty


def single_shot_bmag_temporal_taylor_y_intercepts(refPos, shot, filt, window, norm='ref', window2=[0, 13]):
    wfMagCorr, tau = btc.gen_shot_pos_bmag_correlation(
        refPos=refPos, pos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tof.get_tau_spatial(tau, shot)
    startDex = bmx.finding_Index_Time(separations*1e-6, window2[0])
    endDex = bmx.finding_Index_Time(separations*1e-6, window2[1])
    yTSArray, maxSeparationArray, uncertainty, _ = y_intercepts(
        wfMagCorr[startDex:endDex], separations[startDex:endDex])
    return yTSArray, maxSeparationArray, uncertainty


def single_shot_abs_bmag_temporal_taylor_y_intercepts(refPos, shot, filt, window, norm='ref', window2=[0, 13]):
    """ Outputs the y_intercept array from the temporal correlations of the absolute Bmag with respect
    to shot.
    """
    wfMagCorr, tau = btc.gen_shot_pos_bmag_correlation(
        refPos=refPos, pos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tof.get_tau_spatial(tau, shot)
    startDex = bmx.finding_Index_Time(separations*1e-6, window2[0])
    endDex = bmx.finding_Index_Time(separations*1e-6, window2[1])
    yTSArray, maxSeparationArray, uncertainty, _ = y_intercepts(
        np.abs(wfMagCorr[startDex:endDex]), separations[startDex:endDex])
    return yTSArray, maxSeparationArray, uncertainty


def single_shot_abs_bmag_taylor_y_intercepts(refPos, shot, filt, window, norm='ref'):
    """ Outputs the TS y-intercepts from bmag correlations of shot referenced 
    to refPos windowed and filt."""
    wfMagCorr = bsc.gen_shot_bmag_spatial_correlation(
        refPos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tools.gen_separation_array(refPos=refPos)
    y_interceptsTSArray, maxSeparationArray, uncertainty, slopeArray = y_intercepts(
        np.abs(wfMagCorr), separations)
    return y_interceptsTSArray, maxSeparationArray, uncertainty, slopeArray


def single_shot_bmag_taylor_y_intercepts(refPos, shot, filt, window, norm='ref'):
    """ Outputs the biased Taylor scales from bmag correlations of shot referenced 
    to refPos windowed and filt."""
    wfMagCorr = bsc.gen_shot_bmag_spatial_correlation(
        refPos=refPos, shot=shot, filt=filt, window=window, norm=norm)
    separations = tools.gen_separation_array(refPos=refPos)
    y_interceptsTSArray, maxSeparationArray, uncertainty, slopeArray = y_intercepts(
        wfMagCorr, separations)
    return y_interceptsTSArray, maxSeparationArray, uncertainty, slopeArray


def biased_taylor_scales(correlation, separation, uncertainties=None):
    """ Outputs the biased Taylor scales, uncertainties and max separation array"""
    biasedTaylorArray, errorArray, maxSeparationArray, _ = re.gen_ordered_sequence(
        correlation, separation, uncertainties)
    return biasedTaylorArray, maxSeparationArray, errorArray[0]


def y_intercepts(correlation, separation, uncertainties=None):
    """ Outputs the y-intercepts from the linear fits on the biased Taylor scales."""
    yInterceptArray, errorArray, maxSeparationArray, slope_Array = re.richardson_extrapolation(
        correlation, separation, uncertainties)
    return yInterceptArray, maxSeparationArray, errorArray[0], slope_Array
