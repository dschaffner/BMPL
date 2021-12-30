import numpy as np
import matplotlib.pylab as plt
import tools
import bmx_data as data03
import bmag_norms as norms
import bmx as BMX


def gen_bmag_spatial_correlation_plot(refPos, filt, window, norm='ref'):
    """ Plot the bmag spatial correlation plot based on specified normalization. """
    fig = plt.figure()
    gs = fig.add_gridspec(1, 1)
    axmag = fig.add_subplot(gs[0, 0])

    plt.title(f'Starting Probe {refPos}')

    separation = tools.gen_separation_array(refPos)
    magShotArray = gen_bmag_spatial_correlation(
        refPos, filt, window, norm)

    axmag.plot(separation, magShotArray, 'o', color='red')
    # axmag.set_ylim(0, 1.1)
    axmag.set_ylabel(r'$R(s)$')
    axmag.set_xlabel(r'Separation (cm)')
    plt.show()
    # plt.savefig(f'Figure/magSpaCorr/magSpaCorrelationP5.png')
    plt.close()

    return None


def save_bmag_spatial_correlation_plot(filename, refPos, filt, window, norm='ref'):
    """ Plot the bmag spatial correlation plot based on specified normalization. """
    fig = plt.figure()
    gs = fig.add_gridspec(1, 1)
    axmag = fig.add_subplot(gs[0, 0])

    plt.title(f'Starting Probe {refPos} Standard Deviation Normalization')

    separation = tools.gen_separation_array(refPos)
    magShotArray = gen_bmag_spatial_correlation(
        refPos, filt, window, norm)

    axmag.plot(separation, magShotArray, 'o', color='red')
    # axmag.set_ylim(0, 1.1)
    axmag.set_ylabel(r'$R(s)$')
    axmag.set_xlabel(r'Separation (cm)')
    # plt.show()
    plt.savefig(filename)
    plt.close()

    return None


def gen_shot_bmag_spatial_correlation_plot(refPos, shot, filt, window, norm='ref'):
    """ Plot the bmag spatial correlation plot of shot based on specified normalization. """
    fig = plt.figure()
    gs = fig.add_gridspec(1, 1)
    axmag = fig.add_subplot(gs[0, 0])
    plt.title(f'Shot {shot} Spatial Correlation Normalized by {norm}')

    plt.title(f'Starting Probe {refPos}')

    separation = tools.gen_separation_array(refPos)
    magShotArray = gen_shot_bmag_spatial_correlation(
        refPos, shot, filt, window, norm)

    axmag.plot(separation, magShotArray, 'o', color='red')
    # axmag.set_ylim(0, 1.1)
    axmag.set_ylabel(r'$R(s)$')
    axmag.set_xlabel(r'Separation (cm)')
    plt.show()
    # plt.savefig(f'Figure/magSpaCorr/magSpacorrelationP{refPos}shot{shot}.png')
    plt.close()

    return None


def gen_bmag_spatial_correlation(refPos, filt, window, norm='ref'):
    """ Outputs the spatial correlation array based on the bmag data. """
    normArray = norms.which_norm(norm, refPos, window, filt)

    def norm_value(norm, posIndex, shot):
        """ Wrapper function to make choosing the norm value easier"""
        if norm.lower() == 'ref':
            normValue = normArray[shot]
        else:
            normValue = normArray[posIndex, shot]

        return normValue

    magShotArray = np.asarray([])
    for pos2Index, pos2 in enumerate([*range(refPos, 16, 2)]):
        for shot in range(0, 25):
            Bmag1, _ = tools.get_windowed_filtered_bmag(
                refPos, shot, window, filt)
            Bmag2, time = tools.get_windowed_filtered_bmag(
                pos2, shot, window, filt)
            corr_MagHold, _ = tools.correlation(
                Bmag1, Bmag2, time, mode='valid')
            corr_MagHold = np.abs(corr_MagHold)
            normValue = norm_value(norm, pos2Index, shot)
            if (shot == 0):
                corr_Mag = corr_MagHold/(25*normValue)
            else:
                corr_Mag += corr_MagHold/(25*normValue)

        magShotArray = np.append(
            magShotArray, corr_Mag)
    return magShotArray


def gen_shot_bmag_spatial_correlation(refPos, shot, filt=3e4, window=[75, 150], norm='ref'):

    normArray = norms.which_norm(norm, refPos, window, filt)

    def norm_value(norm, posIndex, shot):
        """ Wrapper function to make choosing the norm value easier"""
        if norm.lower() == 'ref':
            normValue = normArray[shot]
        else:
            normValue = normArray[posIndex, shot]

        return normValue
    magShotArray = np.asarray([])

    for pos2Index, pos2 in enumerate([*range(refPos, 16, 2)]):
        Br_1, Bt_1, Bz_1, time = data03.load_data(refPos, shot)
        Br_2, Bt_2, Bz_2, time = data03.load_data(pos2, shot)

        if (pos2 == refPos):
            start_Dex = BMX.finding_Index_Time(time*1e-6, window[0])
            if (window[1] != 190):
                end_Dex = BMX.finding_Index_Time(time*1e-6, window[1])
            else:
                end_Dex = -1

        Bmag1 = np.sqrt(Br_1**2 + Bt_1**2 + Bz_1**2)
        Bmag2 = np.sqrt(Br_2**2 + Bt_2**2 + Bz_2**2)

        Bmag1 = tools.HPF(Bmag1, filt)
        Bmag2 = tools.HPF(Bmag2, filt)

        time = time[start_Dex:end_Dex]

        Bmag1 = Bmag1[start_Dex:end_Dex]
        Bmag2 = Bmag2[start_Dex:end_Dex]

        corrMagHold, _ = tools.correlation(
            Bmag1, Bmag2, time, mode='valid')
        # corr_MagHold = np.abs(corr_MagHold)
        normValue = norm_value(norm=norm, posIndex=pos2Index, shot=shot)
        normCorrMagHold = corrMagHold/normValue

        magShotArray = np.append(
            magShotArray, normCorrMagHold)
    return magShotArray


def gen_shot_pos_bmag_correlation(refPos, pos, shot, filt=3e4, window=[75, 150], norm='ref', normalized=True):
    """ This function outputs the bmag correlation of shot at pos
        and the corresponding time delay."""
    normArray = norms.which_norm(norm, refPos, window, filt)

    def norm_value(norm, posIndex, shot):
        """ Wrapper function to make choosing the norm value easier"""
        if norm.lower() == 'ref':
            normValue = normArray[shot]
        else:
            normValue = normArray[posIndex, shot]

        return normValue

    bmagRef, time = tools.get_windowed_filtered_bmag(
        pos=refPos, shot=shot, window=window, filt=filt)
    bmagPos, _ = tools.get_windowed_filtered_bmag(pos, shot, window, filt)

    magCorr, tau = tools.correlation(bmagRef, bmagPos, time)

    if normalized:
        posIndex = tools.get_pos_index(refPos, pos)
        normMagCorr = magCorr/norm_value(norm, posIndex, shot)
    else:
        normMagCorr = magCorr

    return normMagCorr, tau


def gen_shot_pos_unproccessed_bmag_correlation(refPos, pos, shot, norm='ref', normalized=True):
    """ Outputs the two-time correlation and time delay of Bmags from refPos and pos for shot.
    This correlation is for the time-of-flight analysis and will not be windowed or filtered.
    """
    normArray = norms.which_norm_tof(norm, refPos)

    def norm_value(norm, posIndex, shot):
        """ Wrapper function to make choosing the norm value easier"""
        if norm.lower() == 'ref':
            normValue = normArray[shot]
        else:
            normValue = normArray[posIndex, shot]

        return normValue
    refBmag, _ = tools.get_bmag(pos=refPos, shot=shot)
    bmag, time = tools.get_bmag(pos=pos, shot=shot)

    magCorr, tau = tools.correlation(refBmag, bmag, time)

    if normalized:
        posIndex = tools.get_pos_index(refPos, pos)
        normMagCorr = magCorr/norm_value(norm, posIndex, shot)
    else:
        normMagCorr = magCorr

    return normMagCorr, tau


def gen_shot_pos_unfiltered_bmag_correlation(refPos, pos, shot, window=[75, 150], norm='ref', normalized=False):
    """ Output the correlation and time-delay of windowed Bmag correlations between the specified
    probes.
    """

    normArray = norms.which_windowed_norm(norm, refPos, window)

    def norm_value(norm, posIndex, shot):
        """ Wrapper function to make choosing the norm value easier"""
        if norm.lower() == 'ref':
            normValue = normArray[shot]
        else:
            normValue = normArray[posIndex, shot]

        return normValue
    refBmag, _ = tools.get_windowed_bmag(pos=refPos, shot=shot, window=window)
    bmag, time = tools.get_windowed_bmag(pos=pos, shot=shot, window=window)
    # create new function and delete
    #refBmag -= np.mean(refBmag)
    #bmag -= np.mean(bmag)

    magCorr, tau = tools.correlation(
        refBmag, bmag, time, normalized=normalized)

    if normalized:
        posIndex = tools.get_pos_index(refPos, pos)
        normMagCorr = magCorr/norm_value(norm, posIndex, shot)
    else:
        normMagCorr = magCorr

    return normMagCorr, tau


# save_bmag_spatial_correlation_plot(filename='Figure/magSpaCorr/startP5/magCorrelationstdNormAvg.png',
#                                   refPos=5, filt=3e4, window=[75, 150], norm='std')
"""gen_shot_bmag_spatial_correlation_plot(
    refPos=5, shot=1, filt=3e4, window=[75, 150], norm='ref')"""
