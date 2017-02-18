import numpy as np

def mad_clipper(arr,mad_num=None):
    """
    Data clipper that removes values from an array of values using the Minimum Absolute Deviation
    as the centering function instead of standard deviation

    If mad_num = None then default is mad_num = 3

    Function returns masking array the same size as input arr
    """
    arr = np.array(arr)

    arr_median = np.median(arr)

    if mad_num == None:
        mad_num = 3.


    """
    This factor of 1.4826 is the typical correction factor that allows the MAD value for a sample to 
    be used as the standard deviation of the population.
    """
    MAD = 1.4826*np.median(np.abs(arr-arr_median))
    mad_values = np.abs(arr-arr_median)
    num_mads = mad_values / MAD

    bound_cond = (num_mads <= 3)
    masked_arr = bound_cond

    return masked_arr
