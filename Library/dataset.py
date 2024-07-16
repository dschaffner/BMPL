# List of input and output functions for 03102020 dataset, this script will need to be update for each dataset.
import numpy as np
import h5py


def load_data(pos, shot):
    ####################################################################
    """ Obtaining the data from the HDF5 file."""
    ####################################################################
    filename = '/Volumes/cacs_resear/Data/2020/03102020/HDF5/03102020.h5'
    f = h5py.File(filename, 'r')
    posList = np.asarray([1, 3, 5, 7, 9, 11, 13, 15])
    keyList = ['pos1/b', 'pos3/b', 'pos5/b', 'pos7/b', 'pos9/b', 'pos11/b',
               'pos13/b', 'pos15/b']

    time = f['time/timeB_us'][()]
    ind = np.where(pos == posList)[0]
    posD = f[keyList[int(ind)]]
    Br = posD['r'][()]
    Bt = posD['theta'][()]
    Bz = posD['z'][()]
    f.close()
    return Br[shot], Bt[shot], Bz[shot], time


def load_bdot_data(pos, shot):
    ####################################################################
    """ Obtaining the data from the HDF5 file."""
    ####################################################################
    filename = '/Volumes/cacs_resear/Data/2020/03102020/HDF5/03102020.h5'
    f = h5py.File(filename, 'r')
    posList = np.asarray([1, 3, 5, 7, 9, 11, 13, 15])
    keyList = ['pos1/bdot', 'pos3/bdot', 'pos5/bdot', 'pos7/bdot', 'pos9/bdot', 'pos11/bdot',
               'pos13/bdot', 'pos15/bdot']

    time = f['time/time_us'][()]
    ind = np.where(pos == posList)[0]
    posD = f[keyList[int(ind)]]
    Brdot = posD['r'][()]
    Btdot = posD['theta'][()]
    Bzdot = posD['z'][()]
    f.close()
    return Brdot[shot], Btdot[shot], Bzdot[shot], time


def get_bank_voltage_current():
    """ Outputs the discharge voltage, current and time."""
    filename = '/Volumes/cacs_resear/Data/2020/03102020/HDF5/bankData.h5'
    f = h5py.File(filename, 'r')
    current = f['/data/curr'][()]
    voltage = f['/data/volt'][()]
    time = f['/time/time'][()]
    f.close()
    return current, voltage, time