import numpy as np
import scipy.signal as sps
import scipy.ndimage as ndi
import scipy.sparse as spr


def butter_highpass(signals, cutoff, fs, order=6, axis=-1):
    """ Forward-backward filter inpute signals with butterworth kernel"""
    return sps.filtfilt(*sps.butter(order,
                                    cutoff / (0.5 * fs),
                                    btype='high',
                                    analog=False),
                        signals,
                        axis=axis)


def gaussian_bandpass(signals, kernel_len, kernel_sigma, axis=-1):
    """ Convolve inputs signals with 0 mean gaussian kernel"""
    return ndi.convolve1d(signals,
                          gaussian_kernel(kernel_len,
                                          kernel_sigma),
                          axis=axis)


def gaussian_kernel(length, sigma):
    """
    Returns a 1D gaussian filter of specified length and stdv for convolution
    """
    n = (length - 1.) / 2.
    kernel = np.exp((-1 * np.power(np.arange(-n, n + 1), 2)) / (2. * sigma**2))
    kernel[kernel < np.finfo(kernel.dtype).eps * kernel.max()] = 0
    sumh = kernel.sum()
    if sumh != 0:
        kernel /= sumh
    return kernel - kernel.mean()


def detrend(signals, degree=2):
    """ Substracts the order k trend from each input pixel """
    _, len_signal = signals.shape
    X = np.array([np.power(np.arange(len_signal), k)
                  for k in range(degree + 1)]).T
    Hat = X.dot(np.linalg.inv(X.T.dot(X))).dot(X.T)
    return signals - signals.dot(Hat)


def HALS4activity(Yr, A, C, K, nb, iters=2):
    U = A.T.dot(Yr)
    V = A.T.dot(A)
    for _ in range(iters):
        for m in range(K):  # neurons
            C[m] = np.clip(C[m] + (U[m] - V[m].dot(C)) /
                           V[m, m], -np.inf, np.inf)
        for m in range(nb):  # background
            C[K + m] = np.clip(C[K + m] + (U[K + m] - V[K + m].dot(C)) /
                               V[K + m, K + m], -np.inf, np.inf)
    return C


def HALS4shape(Yr, A, C, K, nb, ind_A, iters=2):
    U = C.dot(Yr.T)
    V = C.dot(C.T)
    for _ in range(iters):
        for m in range(K):  # neurons
            ind_pixels = np.squeeze(ind_A[:, m].toarray())
            A[ind_pixels, m] = np.clip(
                A[ind_pixels, m] + ((U[m, ind_pixels] - V[m].dot(A[ind_pixels].T)) / V[m, m]), 0, np.inf)
        for m in range(nb):  # background
            A[:, K + m] = np.clip(A[:, K + m] + ((U[K + m] -
                                                  V[K + m].dot(A.T)) / V[K + m, K + m]), 0, np.inf)
    return A


def hals(Y, A, C, b, f, bSiz=3, maxIter=5):
    """ Hierarchical alternating least square method for solving NMF problem
    Y = A*C + b*f
    input:
    ------
       Y:      d1 X d2 [X d3] X T, raw data.
        It will be reshaped to (d1*d2[*d3]) X T in this
       function
       A:      (d1*d2[*d3]) X K, initial value of spatial components
       C:      K X T, initial value of temporal components
       b:      (d1*d2[*d3]) X nb, initial value of background spatial component
       f:      nb X T, initial value of background temporal component
       bSiz:   int or tuple of int
        blur size. A box kernel (bSiz X bSiz [X bSiz]) (if int) or bSiz (if tuple) will
        be convolved with each neuron's initial spatial component, then all nonzero
       pixels will be picked as pixels to be updated, and the rest will be
       forced to be 0.
       maxIter: maximum iteration of iterating HALS.
    output:
    -------
        the updated A, C, b, f
    @Author: Johannes Friedrich, Andrea Giovannucci
    See Also:
        http://proceedings.mlr.press/v39/kimura14.pdf
    """

    # smooth the components
    dims, T = np.shape(Y)[:-1], np.shape(Y)[-1]
    K = A.shape[1]  # number of neurons
    nb = b.shape[1]  # number of background components
    if isinstance(bSiz, (int, float)):
        bSiz = [bSiz] * len(dims)
    ind_A = ndi.uniform_filter(np.reshape(
        A, dims + (K,), order='F'), size=bSiz + [0])
    ind_A = np.reshape(ind_A > 1e-10, (np.prod(dims), K), order='F')
    ind_A = spr.csc_matrix(ind_A)  # indicator of nonnero pixels

    Ab = np.c_[A, b]
    Cf = np.r_[C, f.reshape(nb, -1)]
    for _ in range(maxIter):
        Cf = HALS4activity(np.reshape(Y, (np.prod(dims), T), order='F'),
                           Ab, Cf, K, nb)
        Ab = HALS4shape(np.reshape(Y, (np.prod(dims), T), order='F'),
                        Ab, Cf, K, nb, ind_A)

    return Ab[:, :-nb], Cf[:-nb], Ab[:, -nb:], Cf[-nb:].reshape(nb, -1)

######## Some general utils