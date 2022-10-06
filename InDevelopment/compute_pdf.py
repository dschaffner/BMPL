

def get_increment(signal1, signal2, delay_index):
    """ This function takes two signals and compute the difference between the corresponding field values with respect to delay_index.
    Inputs:
        delay_index: integer
        signal1: N 1D numpy array 
        signal2: N 1D numpy array
    Outputs:
        increment_array: N - delay_index 1D numpy array 
    """
    return [signal2[:-delay_index] - signal1[delay_index:]]


def get_delay_index(delay, dt=1/(125e6)):
    """ Function that gets the delay index from delay time value
        Inputs:
            delay: float (seconds)
            dt: sampling time
        Outputs:
            delay_index: integer
    """
    return [int(delay/dt)]
