# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: nonecheck=False

import numpy as np
cimport numpy as np
cimport cython

from cpython cimport bool

np.import_array()

cdef extern from "C/l1tf.c":
    int call_l1tf "l1tf" (const int n,
                       const double *y,
                       const double lambda_,
                       double *x,
                       const int MAXITER,
                       const int verbose) nogil

    double call_l1tf_lambdamax "l1tf_lambdamax"(const int n,
                                                double *y,
                                                const int verbose) nogil

def l1_tf(double[::1] data_,
          double lambda_,
          bool max_lambda,
          int MAXITER,
          int verbose):
    """
    Solve L1 trend filter via primal-dual interior point method

    minimize rescale*||data_-x||_2 ^2 + lambda ||z||_1
    subject to  z=Dx

    Parameters:
    ----------
    data_:      np.array (L,)
                data sequence 
    lambda_:    float
                regularization parameter
    max_lambda: boolean
                find the upper bound for regularization parameter
                scale this value by given lambda_
    rescale:    float
                scaling parameter for L2 norm
    verbose:    int {0,1}
                flag
    Output:
    -------
    data_hat :  np.array (L,)
                data sequence
    """
    
    cdef np.intp_t data_length = data_.shape[0]
    cdef np.double_t[::1] data_hat = np.empty(data_length, dtype=np.double)
    cdef double lambda_max
    cdef int iter_status

    if max_lambda:
        lambda_max = call_l1tf_lambdamax(data_length,
                                         &data_[0],
                                         verbose)
        if lambda_max < 0:
            raise RuntimeError("Lambda < 0")
        lambda_ *= lambda_max

    with nogil:
        iter_status = call_l1tf(data_length,
                                &data_[0],
                                lambda_,
                                &data_hat[0],
                                MAXITER,
                                verbose)

    if iter_status < 0:
        raise RuntimeError("Error")

    return data_hat


def l1tf_lambda_max(double[::1] data_,
                    int verbose):
    """
    Upper bound for regularization parameter

    Parameters:
    ----------
    data_:      np.array (L,)
                data sequence 
    verbose:    int {0,1}
                flag

    Outputs:
    --------
    lambda_max: float
                lambda upper bound
    """

    cdef np.intp_t data_length = data_.shape[0]
    cdef double lambda_max

    lambda_max = call_l1tf_lambdamax(data_length,
                                     &data_[0],
                                     verbose)

    if lambda_max < 0:
        raise RuntimeError("Lambda < 0")

    return lambda_max
