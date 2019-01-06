import numpy as np
import scipy as sp

def remove_trend(Y_rm,detrend_option='linear'):
    mean_pixel = Y_rm.mean(axis=1, keepdims=True)
    Y_rm2 = Y_rm - mean_pixel
    # Detrend
    if detrend_option=='linear':
        detr_data = sp.signal.detrend(Y_rm2,axis=1,type='l')
    #elif detrend_option=='quad':
        #detr_data = detrend(Y_rm)
    else:
        print('Add option')
    Y_det = detr_data + mean_pixel
    offset = Y_rm - Y_det
    return Y_det, offset


def unpad(x):
    """
    Given padded matrix with nan
    Get rid of all nan in order (row, col)

    Parameters:
    ----------
    x:          np.array
                array to unpad (all nan values)

    Outputs:
    -------
    x:          np.array
                unpaded array (will not contain nan values)
                dimension might be different from input array
    """
    x = x[:, ~np.isnan(x).all(0)]
    x = x[~np.isnan(x).all(1)]
    return x


def pad(array, reference_shape, offsets, array_type=np.nan):
    """
    Pad array wrt reference_shape exlcluding offsets with dtype=array_type

    Parameters:
    ----------
    array:          np.array
                    array to be padded
    reference_shape:tuple
                    size of narray to create
    offsets:        tuple
                    list of offsets (number of elements must be equal
                    to the dimension of the array)
                    will throw a ValueError if offsets is too big and the
                    reference_shape cannot handle the offsets
    array_type:     dtype
                    data type to pad array with.

    Outputs:
    -------
    result:         np.array (reference_shape)
                    padded array given input
    """

    # Create an array of zeros with the reference shape
    result = np.ones(reference_shape) * array_type
    # Create a list of slices from offset to offset + shape in each dimension
    insertHere = [slice(offsets[dim], offsets[dim] + array.shape[dim])
                  for dim in range(array.ndim)]
    # Insert the array in the result at the specified offsets
    result[insertHere] = array
    return result

def nextpow2(value):
    """
    Extracted from
caiman.source_extraction.cnmf.deconvolution import axcov

    Find exponent such that 2^exponent is >= abs(value).

    Parameters:
    ----------
    value : int

    Returns:
    -------
    exponent : int
    """

    exponent = 0
    avalue = np.abs(value)
    while avalue > np.power(2, exponent):
        exponent += 1
    return exponent


def axcov(data, maxlag=10):
    """
    Edited from cnmf.deconvolution
    Compute the autocovariance of data at lag = -maxlag:0:maxlag

    Parameters:
    ----------
    data : array
        Array containing fluorescence data

    maxlag : int
        Number of lags to use in autocovariance calculation

    Output:
    -------
    axcov : array
        Autocovariances computed from -maxlag:0:maxlag
    """

    data = data - np.mean(data)
    T = len(data)
    bins = np.size(data)
    xcov = np.fft.fft(data, np.power(2, nextpow2(2 * bins - 1)))
    xcov = np.fft.ifft(np.square(np.abs(xcov)))
    xcov = np.concatenate([xcov[np.arange(xcov.size - maxlag, xcov.size)],
                           xcov[np.arange(0, maxlag + 1)]])
    return np.real(np.divide(xcov, T))