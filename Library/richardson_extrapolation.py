import numpy as np
from scipy.optimize import curve_fit


def gen_ordered_sequence(correlations, separations, uncertainties, fit_type='parabolic'):
    """ Outputs sequence of fits"""

    taylor_List = []
    slope_List = []
    error_List = []
    num_List = []
    step_size = separations[1] - separations[0]

    def parabolic_Fit(z, L):
        ###################################################
        """
            The functional form of the natural log of the
            gaussian, the form is not general but specific
            to this situation.
        """
        ###################################################
        nt_Corr = 1 - z**2/(2*L**2)

        return nt_Corr

    def linear_Fit(z, m, b):

        linear = m*z + b

        return linear

    if fit_type.lower() == 'parabolic':
        for i in range(1, separations.size - 1):
            popt, pcov = curve_fit(
                parabolic_Fit, separations[:i + 2], correlations[:i + 2])
            taylor_List.append(popt[0])
            error_List.append(pcov[0])
            num_List.append(separations[i]/step_size + 1)
    elif fit_type.lower() == 'linear':
        for i in range(1, separations.size):
            popt, pcov = curve_fit(
                linear_Fit, separations[:i+1], correlations[:i+1], sigma=uncertainties[:i+1])
            taylor_List.append(popt[1])
            slope_List.append(popt[0])
            error_List.append(pcov[1, 1])

            num_List.append(separations[i])

    taylor_Array = np.asarray(taylor_List)
    num_Array = np.asarray(num_List)
    error_Array = np.transpose(np.sqrt(np.asarray(error_List)))
    slope_Array = np.asarray(slope_List)

    return taylor_Array, error_Array, num_Array, slope_Array


def richardson_extrapolation(correlations, separations, uncertainties):
    """ Outputs the final ordered sequence in the richardson extrapolation technique, with uncertainties
    and the ordering array."""
    taylor_array, error_array, max_separation_array, _ = gen_ordered_sequence(
        correlations, separations, uncertainties, fit_type='parabolic')

    final_taylor_array, final_error_array, final_max_separation_array, slope_Array = gen_ordered_sequence(
        taylor_array, max_separation_array, error_array[0], fit_type='linear')

    return final_taylor_array, final_error_array, final_max_separation_array, slope_Array
