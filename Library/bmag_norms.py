import numpy as np
import tools
import bmx as BMX


def which_norm(norm, refPos, window, filt, max_probe_number=8):
    """ A decider function that chooses which normalization to return."""
    if isinstance(norm, str):
        if norm.lower() == 'ref':
            normArray = shot_bmag_refnorm(refPos, window, filt)
        elif norm.lower() == 'std':
            normArray = shot_bmag_stdnorm(
                refPos, window, filt, max_probe_number=max_probe_number)
        else:
            raise TypeError(f'{norm} must be either ref or std')
    else:
        raise TypeError('The norm variable must be a string')

    return normArray


def shot_bmag_refnorm(refPos, window, filt):
    """ Outputs the Bmag normalization per shot based on the reference probe."""
    normArray = np.zeros(25)
    for shot in range(0, 25):

        wfBmag, time = tools.get_windowed_filtered_bmag(
            pos=refPos, shot=shot, window=window, filt=filt)

        norm_Mag, _ = tools.correlation(wfBmag, wfBmag, time, mode='valid')

        normArray[shot] = norm_Mag
    return normArray


def shot_bmag_stdnorm(refPos, window, filt, max_probe_number=8):
    """ Outputs the tools.correlation norm array w.r.t shot. The norm is product of 
    the standard deviations of the two windowed and filtered bmag signals. """

    maxRange = tools.gen_max_probe_range(
        refPos, max_probe_number=max_probe_number)
    type(maxRange)
    normWFArray = np.zeros([maxRange, 25])
    posArray = [*range(refPos, 16, 2)]

    for posIndex, pos in enumerate(posArray):
        for shot in range(0, 25):

            Br_1, Bt_1, Bz_1, time = tools.load_data(refPos, shot)
            Br_2, Bt_2, Bz_2, time = tools.load_data(pos, shot)

            start_Dex = BMX.finding_Index_Time(time*1e-6, window[0])
            end_Dex = BMX.finding_Index_Time(time*1e-6, window[1])

            Bmag1 = np.sqrt(Br_1**2 + Bt_1**2 + Bz_1**2)
            Bmag2 = np.sqrt(Br_2**2 + Bt_2**2 + Bz_2**2)
            fBmag1 = tools.HPF(Bmag1, filt)
            fBmag2 = tools.HPF(Bmag2, filt)

            wfBmag1 = fBmag1[start_Dex:end_Dex]
            wfBmag2 = fBmag2[start_Dex:end_Dex]

            normWFBmag = np.sqrt(np.var(wfBmag1)*np.var(wfBmag2))

            normWFArray[posIndex, shot] = normWFBmag

    return normWFArray


def which_norm_tof(norm, refPos, max_probe_number=8):
    """ A decider function that chooses which normalization to return."""
    if isinstance(norm, str):
        if norm.lower() == 'ref':
            normArray = shot_bmag_refnorm_tof(refPos)
        elif norm.lower() == 'std':
            normArray = shot_bmag_stdnorm_tof(
                refPos, max_probe_number=max_probe_number)
        else:
            raise TypeError(f'{norm} must be either ref or std')
    else:
        raise TypeError('The norm variable must be a string')

    return normArray


def shot_bmag_refnorm_tof(refPos):
    """ Outputs the Bmag normalization per shot based on the reference probe."""
    normArray = np.zeros(25)
    for shot in range(0, 25):

        bmag, time = tools.get_bmag(
            pos=refPos, shot=shot)

        norm_Mag, _ = tools.correlation(bmag, bmag, time, mode='valid')

        normArray[shot] = norm_Mag
    return normArray


def shot_bmag_stdnorm_tof(refPos, max_probe_number=8):
    """ Outputs the tools.correlation norm array w.r.t shot. The norm is product of 
    the standard deviations of the two windowed and filtered bmag signals. """

    maxRange = tools.gen_max_probe_range(
        refPos, max_probe_number=max_probe_number)
    type(maxRange)
    normArray = np.zeros([maxRange, 25])
    posArray = [*range(refPos, 16, 2)]

    for posIndex, pos in enumerate(posArray):
        for shot in range(0, 25):

            refBmag, _ = tools.get_bmag(pos=refPos, shot=shot)
            bmag, _ = tools.get_bmag(pos=pos, shot=shot)

            normBmag = np.sqrt(np.var(refBmag)*np.var(bmag))

            normArray[posIndex, shot] = normBmag

    return normArray


def which_windowed_norm(norm, refPos, window, max_probe_number=8):
    """ A decider function that chooses which normalization to return."""
    if isinstance(norm, str):
        if norm.lower() == 'ref':
            normArray = shot_windowed_bmag_refnorm(refPos, window)
        elif norm.lower() == 'std':
            normArray = shot_windowed_bmag_stdnorm(
                refPos, window=window, max_probe_number=max_probe_number)
        else:
            raise TypeError(f'{norm} must be either ref or std')
    else:
        raise TypeError('The norm variable must be a string')

    return normArray


def shot_windowed_bmag_stdnorm(refPos, window, max_probe_number=8):
    """ Outputs the tools.correlation norm array w.r.t shot. The norm is product of 
    the standard deviations of the two windowed bmag signals. """

    maxRange = tools.gen_max_probe_range(
        refPos, max_probe_number=max_probe_number)
    type(maxRange)
    normArray = np.zeros([maxRange, 25])
    posArray = [*range(refPos, 16, 2)]

    for posIndex, pos in enumerate(posArray):
        for shot in range(0, 25):

            refBmag, _ = tools.get_windowed_bmag(
                pos=refPos, shot=shot, window=window)
            bmag, _ = tools.get_windowed_bmag(
                pos=pos, shot=shot, window=window)

            normBmag = np.sqrt(np.var(refBmag)*np.var(bmag))

            normArray[posIndex, shot] = normBmag

    return normArray


def shot_windowed_bmag_refnorm(refPos, window):
    """ Outputs the windowed Bmag normalization per shot based on the reference probe."""
    normArray = np.zeros(25)
    for shot in range(0, 25):

        bmag, time = tools.get_windowed_bmag(
            pos=refPos, shot=shot, window=window)

        norm_Mag, _ = tools.correlation(bmag, bmag, time, mode='valid')

        normArray[shot] = norm_Mag
    return normArray
