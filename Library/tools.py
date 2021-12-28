import numpy as np
import bmx as BMX
import bmx_data as data03
import h5py
import matplotlib.pylab as plt
import spectrum_wwind as sw

def get_windowed_bmag(pos, shot, window):
    """ Output the windowed bmag and time arrays of {shot} at {pos}."""
    bmag, time = get_bmag(pos, shot)

    start_dex = BMX.finding_Index_Time(time*1e-6, window[0])
    end_dex = BMX.finding_Index_Time(time*1e-6, window[1])

    wbmag = bmag[start_dex:end_dex]
    time = time[start_dex:end_dex]

    return wbmag, time


def get_windowed_filtered_bmag(pos, shot, window, filt):
    """ Output the windowed filtered bmag and time arrays of {shot} at {pos}."""
    filt_bmag, time = get_filtered_bmag(pos, shot, filt)

    start_dex = BMX.finding_Index_Time(time*1e-6, window[0])
    end_dex = BMX.finding_Index_Time(time*1e-6, window[1])

    wfilt_bmag = filt_bmag[start_dex:end_dex]
    time = time[start_dex:end_dex]

    return wfilt_bmag, time


def get_filtered_bmag(pos, shot, filt):
    """ Output the filtered bmag and time arrays of {shot} at {pos}."""
    bmag, time = get_bmag(pos, shot)
    filt_bmag = HPF(bmag, filt)
    return filt_bmag, time


def get_bmag(pos, shot):
    """ Output the bmag and time arrays of shot at pos."""
    br, bt, bz, time = data03.load_data(pos, shot)
    bmag = np.sqrt(br**2 + bt**2 + bz**2)
    return bmag, time


def get_bank_voltage_current():
    """ Outputs the discharge voltage, current and time."""
    filename = '/Volumes/cacs_resear/Data/2020/03102020/HDF5/bankData.h5'
    f = h5py.File(filename, 'r')
    current = f['/data/curr'][()]
    voltage = f['/data/volt'][()]
    time = f['/time/time'][()]
    f.close()
    return current, voltage, time


def HPF(data, filt):
    """ Outputs the filtered data """
    data = BMX.high_Pass_Filter(data, filter_Freq=filt)
    return data


def correlation(sing_1, sing_2, time, normalized=True, mode='same'):
    """ Normalization based on the number of point considered """

    corr = np.correlate(sing_2, sing_1, mode=mode)
    dt = time[1]-time[0]
    tau = dt*(np.arange(corr.size) - corr.size/2)

    if normalized:
        corr = corr/(time.size - 1)

    return corr, tau


def plot_components(pos, shot):
    """ Plot the raw magnetic data.
    """
    x, y, z, time = data03.load_data(pos, shot)
    fig = plt.figure()
    fig.suptitle('Raw Fields')
    gs = fig.add_gridspec(3,1)
    axx = fig.add_subplot(gs[0,0])
    axy = fig.add_subplot(gs[1,0])
    axz = fig.add_subplot(gs[2,0])

    axx.plot(time, x)
    axy.plot(time, y)
    axz.plot(time, z)

    plt.show()
    plt.close()

    return None


def gen_separation_array(refPos, max_probe_number=8):
    """ Outputs the separation array w.r.t the reference probe"""
    range = gen_max_probe_range(refPos, max_probe_number=max_probe_number)
    separation = np.arange(0, range)*2.6

    return separation


def gen_max_probe_range(refPos, max_probe_number=8):
    """ The probe number range"""
    startProbeNumber = (refPos - 1)/2
    range = max_probe_number - startProbeNumber
    return int(range)


def get_pos_index(refPos, pos):
    positions = np.arange(refPos, 16, 2)
    indices = np.where(positions == pos)
    posIndex = indices[0][0]
    return posIndex


def get_power_spectrum(data, time, alias = 'hanning'):
    """ Output the power spectrum and frequency
    """

    freq, _, _, pwr, _, _, _, _ = sw.spectrum_wwind(
                    data, time*1e-6, window=alias)

    return freq, pwr