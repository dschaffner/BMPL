import numpy as np
import tools
import matplotlib.pylab as plt
import bmag_norms as norms


def gen_shot_bmag_temporal_correlation_plot(Pos, shot, filt, window, norm='ref'):
    """ Plot the bmag spatial correlation plot of shot based on specified normalization. """
    fig = plt.figure()
    gs = fig.add_gridspec(1, 1)
    axmag = fig.add_subplot(gs[0, 0])
    plt.title(f'Shot {shot} Temporal Correlation Normalized by {norm}')

    plt.title(f'Corr5-{Pos}')

    magShotArray, tau = gen_shot_pos_bmag_correlation(
        5, Pos, shot, filt, window, norm)
    axmag.axhline(y=0, color='black', alpha=0.5)
    axmag.plot(tau, magShotArray, color='red')
    # axmag.set_ylim(0, 1.1)
    axmag.set_ylabel(r'$R(s)$')
    axmag.set_xlabel(r'Separation (cm)')
    plt.show()
    # plt.savefig(f'Figure/magSpaCorr/magSpacorrelationP{refPos}shot{shot}.png')
    plt.close()

    return None


def gen_bmag_temporal_correlation(refPos, filt, window, norm='ref'):
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
    for shot in range(0, 25):
        Bmag1, time = tools.get_windowed_filtered_bmag(
            refPos, shot, window, filt)
        corr_MagHold, tau = tools.correlation(
            Bmag1, Bmag1, time, mode='same')
        corr_MagHold = np.abs(corr_MagHold)
        normValue = norm_value(norm, 0, shot)
        if (shot == 0):
            corr_Mag = corr_MagHold/(25*normValue)
        else:
            corr_Mag += corr_MagHold/(25*normValue)

    magShotArray = np.append(
        magShotArray, corr_Mag)
    return magShotArray, tau


def gen_shot_pos_bmag_correlation(refPos, pos, shot, filt, window, norm='ref', normalized=True):
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
