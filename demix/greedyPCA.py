import numpy
import numpy as np
import scipy as sp

from scipy.stats import norm
from sklearn.utils.extmath import randomized_svd
from sklearn import preprocessing

import concurrent
import cvxpy as cp
import multiprocessing
import itertools
import time

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

# see how this changes with get_noise_fft
#import spatial_filtering as sp_filters

import trefide
import util_plot as uplot
import tools as tools_
import noise_estimator
from l1_trend_filter.l1_tf_C.c_l1_tf import l1_tf# cython

def kurto_one(x):
    """
    Computer kurtosis of rows of x with 2D (d x T)
    Outputs: idx of components with kurtosis > population mean + std
    """
    # Compute kurtosis
    kurt = sp.stats.kurtosis(x, axis=1, fisher=True, bias=True)
    # keep component wrt mean + std
    kurt_th = kurt.mean() + kurt.std()
    keep = np.where(kurt > kurt_th)[0].tolist()
    return keep

def stimulus_segmentation(T,
                          stim_knots=None,
                          stim_delta=10):
    """
    Returns boolean areas to ignore in frame
    """
    if stim_knots is None:
        return np.zeros(T).astype('bool')

    stim_knots = np.where(np.abs(np.diff(stim_knots))>1)[0]

    ignore_frames =np.zeros(T)
    for knot in stim_knots:
        arg_np = (np.ones(stim_delta).cumsum()-stim_delta//2+knot).astype('int')
        ignore_frames[arg_np]=1
    ignore_frames=ignore_frames.astype('bool')
    return ignore_frames


def mean_confidence_interval(data,
                            confidence=0.90,
                            one_sided=False):
    """
    Compute mean confidence interval (CI)
    for a normally distributed population
    _____________

    Parameters:
    ___________
    data:       np.array  (L,)
                input vector from which to calculate mean CI
                assumes gaussian-like distribution
    confidence: float
                confidence level for test statistic
    one_sided:  boolean
                enforce a one-sided test

    Outputs:
    _______
    th:         float
                threshold for mean value at CI
    """
    if one_sided:
        confidence = 1 - 2*(confidence-1)

    _, th = sp.stats.norm.interval(confidence,
                                loc =np.mean(data),
                                scale=data.std())

    return th


def choose_rank(Vt,
                maxlag=10,
                iterate=False,
                confidence=0.90,
                corr=True,
                kurto=False,
                mean_th=None,
                mean_th_factor=1.,
                min_rank=0,
               enforce_both=False):
    """
    Select rank vectors in Vt which pass test statistic(s) enabled
    (e.g. axcov and/or kurtosis)

    __________
    Parameters:

    Vt:         np.array (k x T)
                array of k temporal components lasting T samples
    maxlag:     int
                max correlation lag for correlation null hypothesis in samples
                (e.g. indicator decay in samples)
    iterate:    boolean
                flag to include correlated components iteratively
    confidence: float
                confidence interval (CI) for correlation null hypothesis
    corr:       boolean
                flag to include components which pass correlation null hypothesis
    kurto:      boolean
                flag to include components which pass kurtosis null hypothesis
    mean_th:    float
                threshold employed to reject components according to
                correlation null hypothesis
    mean_th_factor: float
                factor to scale mean_th
                typically set to 2 if greedy=True and mean_th=None or
                if mean_th has not been scaled yet.
    min_rank:   int
                minimum number of components to include in output
                even if no components of Vt pass any test

    Outputs:
    -------

    vtid:       np.array (3,d)
                indicator 3D matrix (corr-kurto-reject) which points which statistic
                a given component passed and thus it is included.
                can vary according to min_rank

    """
    if enforce_both:
        corr= True
        #kurto = True

    n, L = Vt.shape
    vtid = np.zeros(shape=(3, n)) * np.nan

    # Null hypothesis: white noise ACF
    if corr is True:
        if mean_th is None:
            mean_th = wnoise_acov_CI(L,confidence=confidence,maxlag=maxlag)
        mean_th*= mean_th_factor
        keep1 = vector_acov(Vt,
                            mean_th = mean_th,
                            maxlag=maxlag,
                            iterate=iterate,
                            min_rank=min_rank)
    else:
        keep1 = []
    if kurto is True:
        keep2 = kurto_one(Vt)
    else:
        keep2 = []

    keep = list(set(keep1 + keep2))
    loose = np.setdiff1d(np.arange(n),keep)
    loose = list(loose)

    if enforce_both:
        keep1 = np.intersect1d(keep1,keep2)
        keep2 = keep1
        print(len(keep1))


    vtid[0, keep1] = 1  # components stored due to cov
    vtid[1, keep2] = 1  # components stored due to kurto
    vtid[2, loose] = 1  # extra components ignored
    # print('rank cov {} and rank kurto {}'.format(len(keep1),len(keep)-len(keep1)))
    return vtid


def wnoise_acov_CI(L,
                    confidence=0.90,
                    maxlag=10,
                    n=3000,
                    plot_en=False):
    """
    Generate n AWGN vectors lasting L samples.
    Calculate the mean of the ACF of each vector for 0:maxlag
    Return the mean threshold with specified confidence.

    Parameters:
    ----------

    L:          int
                length of vector
    confidence: float
                confidence level for test statistic
    maxlag:     int
                max correlation lag for correlation null hypothesis
                in samples (e.g. indicator decay in samples)
    n:          int
                number of standard normal vectors to generate

    plot_en:    boolean
                plot histogram of pd
    Outputs:
    -------

    mean_th:    float
                value of mean of ACFs of each standard normal vector at CI.
    """

    # th1 = 0
    #print 'confidence is {}'.format(confidence)
    covs_ht = np.zeros(shape=(n,))
    for sample in np.arange(n):
        ht_data = np.random.randn(L)
        covdata = tools_.axcov(ht_data,
                            maxlag)[maxlag:]/ht_data.var()

        covs_ht[sample] = covdata.mean()
        #covs_ht[sample] = np.abs(covdata[1:]).mean()
    #hist, _,_=plt.hist(covs_ht)
    #plt.show()
    mean_th = mean_confidence_interval(covs_ht,confidence)
    #print('th is {}'.format(mean_th))
    return mean_th


def vector_acov(Vt,
                mean_th=0.10,
                maxlag=10,
                iterate=False,
                extra=1,
                min_rank=0,
                verbose=False):
    """
    Calculate auto covariance of row vectors in Vt
    and output indices of vectors which pass correlation null hypothesis.

    Parameters:
    ----------
    Vt:         np.array(k x T)
                row array of compoenents on which to test correlation null hypothesis
    mean_th:    float
                threshold employed to reject components according to correlation null hypothesis
    maxlag:     int
                determined lag until which to average ACF of row-vectors for null hypothesis
    iterate:    boolean
                flag to include components which pass null hypothesis iteratively
                (i.e. if the next row fails, no additional components are added)
    extra:      int
                number of components to add as extra to components which pass null hypothesis
                components are added in ascending order corresponding to order in mean_th
    min_rank:   int
                minimum number of components that should be included
                add additional components given components that (not) passed null hypothesis
    verbose:    boolean
                flag to enable verbose

    Outputs:
    -------
    keep:       list
                includes indices of components which passed the null hypothesis
                and/or additional components added given parameters
    """

    keep = []
    num_components = Vt.shape[0]
    print('mean_th is %s'%mean_th) if verbose else 0
    for vector in range(0, num_components):
        # standarize and normalize
        vi = Vt[vector, :]
        vi =(vi - vi.mean())/vi.std()

        print('vi ~ (mean: %.3f,var:%.3f)'%(vi.mean(),
                                            vi.var())) if verbose else 0

        vi_cov = tools_.axcov(vi,
                            maxlag)[maxlag:]/vi.var()

        if vi_cov.mean() < mean_th:
            if iterate is True:
                break
        else:
            keep.append(vector)
    # Store extra components
    if vector < num_components and extra != 1:
        extra = min(vector*extra,
                    Vt.shape[0])

        for addv in range(1, extra-vector+ 1):
            keep.append(vector + addv)
    # Forcing one components
    if not keep and min_rank>0:
        # vector is empty for once min
        keep.append(0)
        print('Forcing one component') if verbose else 0
    return keep


def compute_svd(M,
                method='vanilla',
                n_components=2,
                n_iter=7,
                random_state=None,
               reconstruct=False):
    """
    Decompose array M given parameters.
    asumme M has been mean_subtracted

    Parameters:
    ----------

    M:          np.array (d xT)
                input array to decompose
    method:     string
                method to decompose M
                ['vanilla','randomized']
    n_components: int
                number of components to extract from M
                if method='randomized'

    Outputs:
    -------

    U:          np.array (d x k)
                left singular vectors of M
                k = n_components if method='randomized'
    Vt:         np.array (k x T)
                right singular vectors of M
                k = n_components if method='randomized'
    s:          np.array (k x 1)
                variance of components
                k = n_components if method='randomized'
    """

    if method == 'vanilla':
        U, s, Vt = numpy.linalg.svd(M, full_matrices=False)

    elif method == 'randomized':
        U, s, Vt = randomized_svd(M,
                                n_components=n_components,
                                n_iter=n_iter,
                                random_state=random_state)

    if reconstruct:
        return U.dot(np.diag(s).dot(Vt))
    return U, s, Vt


def temporal_decimation(data,
                        mb=1):
    """
    Decimate data by mb
    new frame is mean of decimated frames

    Parameters:
    ----------
    data:       np.array (T x d)
                array to be decimated wrt first axis

    mb:         int
                contant by which to decimate data
                (i.e. number of frames to average for decimation)

    Outputs:
    -------
    data0:      np.array (T//mb x d)
                temporally decimated array given data
    """
    data0 = data[:int(len(data)/mb*mb)].reshape(
        (-1, mb) + data.shape[1:]).mean(1).astype('float32')
    return data0


def spatial_decimation(data,
                       dims,
                       ds=1):
    """

    Decimate data by ds
    smaller frame is mean of decimated frames

    Parameters:
    ----------
    data:       np.array (T x d)
                array to be decimated wrt second axis

    ds:         int
                contant by which to decimate data
                (i.e. number of pixels to average for decimation)

    dims:       np.array or tupe (d1,d2,T)
                dimensions of data

    Outputs:
    -------
    data0:      np.array (T x d//ds)
                spatially decimated array given data
    """

    #data0 = data.reshape(len(data0), dims[1]
    # / ds[0], ds[0], dims[2] / ds[1], ds[1]).mean(2).mean(3)
    data0  = data.copy()
    return data0


def denoise_patch(M,
                  maxlag=5,
                  tsub=1,
                  ds=1,
                  noise_norm=False,
                  iterate=False,
                  confidence=0.999,
                  corr=True,
                  kurto=False,
                  tfilt=False,
                  tfide=False,
                  share_th=True,
                  plot_en=False,
                  greedy=True,
                  fudge_factor=0.99,
                  mean_th=None,
                  mean_th_factor=1.15,
                  U_update=False,
                  min_rank=1,
                  verbose=False,
                  pca_method='vanilla',
                  stim_knots=None,
                  stim_delta=200):
    """
    Given single patch, denoise it as outlined by parameters

    Parameters:
    ----------

    M:          np.array (d1xd2xT)
                array to be denoised
    maxlag:     int
                max correlation lag for correlation null hypothesis in samples
                (e.g. indicator decay in samples)
    tsub:       int
                temporal downsample constant
    ds:         int
                spatial downsample constant
    noise_norm: placeholder
    iterate:    boolean
                flag to include correlated components iteratively
    confidence: float
                confidence interval (CI) for correlation null hypothesis
    corr:       boolean
                flag to include components which pass correlation null hypothesis
    kurto:      boolean
                will be removed
                flag to include components which pass kurtosis null hypothesis
    tfilt:      boolean
                flag to temporally filter traces with AR estimate of order p.
    tfide:      boolean
                flag to denoise temporal traces with Trend Filtering
    min_rank:   int
                minimum rank of denoised/compressed matrix
                typically set to 1 to avoid empty output (array of zeros)
                if input array is mostly noise.
    greedy:     boolean
                flag to greedily update spatial and temporal components
                (estimated with PCA) greedily by denoising
                temporal and spatial components
    mean_th_factor: float
                factor to scale mean_th
                typically set to 2 if greedy=True and mean_th=None
                or if mean_th has not been scaled yet.
    share_th:   boolean
                flag to compute a unique thredhold for correlation null hypothesis
                to be used in each tile.
                If false: each tile calculates its own mean_th value.
    fudge_factor: float
                constant to scale estimated noise std st denoising is less
                (lower factor) or more (higher factor) restrictive.
    U_update:   boolean
                flag to (not) update spatial components by imposing L1- constraint.
                True for "round" neurons in 2p.
                False for dendritic data.
    plot_en:    boolean
                flag to enable plots
    verbose:    boolean
                flag to enable verbose
    pca_method: string
                method for matrix decomposition (e.g. PCA, sPCA, rPCA, etc).
                see compute_svd for options

    Outputs:
    -------

    Yd:         np.array (d1 x d2 x T)
                compressed/denoised array given input M
    rlen:       int
                sum of the ranks of all tiles
    """
    dimsM = M.shape
    M = M.reshape((np.prod(dimsM[:2]),dimsM[2]),order='F')

    #print('greedy here 505')
    start = time.time()
    Yd, vtids = denoise_components(M,
                                   maxlag=maxlag,
                                   tsub=tsub,
                                   noise_norm=noise_norm,
                                   iterate=iterate,
                                   confidence=confidence,
                                   corr=corr,
                                   kurto=kurto,
                                   tfilt=tfilt,
                                   tfide=tfide,
                                   mean_th=mean_th,
                                   greedy=greedy,
                                   fudge_factor=fudge_factor,
                                   mean_th_factor=mean_th_factor,
                                   U_update=U_update,
                                   min_rank=min_rank,
                                   plot_en=plot_en,
                                   verbose=verbose,
                                   dims=dimsM,
                                   pca_method=pca_method,
                                   stim_knots=stim_knots,
                                   stim_delta=stim_delta
                                  )

    Yd = Yd.reshape(dimsM, order='F')
    # determine individual rank
    rlen = total_rank(vtids)

    print('\trank:%d\trun_time: %f'%(rlen,time.time()-start))
    return Yd, rlen


def total_rank(vtids,
            verbose=False):
    # determine individual rank
    case1 = ~np.isnan(vtids[0,:])
    if vtids[0,case1].sum()>0:
        ranks = case1
    else:
        ranks = np.nan
    #ranks = np.where(np.logical_or(vtids[0, :] >= 1, vtids[1, :] == 1))[0]
    if np.all(ranks == np.nan):
        rlen = 0
    else:
        rlen = vtids[0,ranks].sum() #len(ranks)
    return rlen


def greedy_spatial_denoiser(Y,
                            V_TF,
                            fudge_factor=1,
                            U_update=False,
                            nus_=None,
                            dims=None,
                            plot_en=False):
    """
    """
    if U_update:
        outs_2 = [c_l1_u_hat(y, V_TF,fudge_factor) for y in Y]
        #outs_2 = update_U_parallel(Y,V_TF,fudge_factor)
        U_hat, nus_ = map(np.asarray,zip(*np.asarray(outs_2)))
    else:
        nus_= np.zeros((V_TF.shape[0]))
        U_hat = np.matmul(Y, np.matmul(V_TF.T, np.linalg.inv(np.matmul(V_TF, V_TF.T))))
    return U_hat, nus_


def greedy_temporal_denoiser(V_hat,
                             solver='SCS',
                             region_thresh_min_pnr=3,
                             region_active_discount=1,
                             lagrange_scaled=True,
                             constraint_segmented=False,
                             fudge_factor=1,
                             verbose=False
                            ):
    """
    Denoise temporal components greedily
    """
    num_components,T =V_hat.shape

    if constraint_segmented:
        outs_ = [trefide.constrained_l1tf(V_hat[idx,:],
                                      solver=solver,
                                      region_thresh_min_pnr=region_thresh_min_pnr,
                                      region_active_discount=region_active_discount,
                                      lagrange_scaled=lagrange_scaled)
             for idx in range(num_components)]

        V_TF, region_indices,lambdas_ = map(np.asarray, zip(*outs_))

    else:
        diff = trefide.difference_operator(num_components)
        #noise_std_ = sp_filters.noise_estimator(V_hat,method='logmexp')
        noise_std_ = noise_estimator.get_noise_fft(V_hat)[0]
        noise_std_ *= fudge_factor

        print('Noise range is %.3e %.3e'%(min(noise_std_),
            max(noise_std_))) if verbose else 0


        outs_ = [c_l1tf_v_hat(V_hat[idx,:],
                              stdv,
                              solver=solver)
                 for idx, stdv in enumerate(noise_std_)]

        V_TF, lambdas_ = map(np.asarray, zip(*np.asarray(outs_)))
        region_indices = None

    #V_TF = preprocessing.normalize(V_TF, norm='l2')
    return V_TF,region_indices,lambdas_


def iterative_temporal_denoiser(Y,
                                U_hat_,
                                V_TF_2,
                                lambdas_=None,
                                region_indices=None,
                                verbose=False,
                                plot_en=False,
                                constraint_segmented=True,
                                fudge_factor=1,
                                solver='SCS',
                                lagrange_scaled=True,
                                region_thresh_min_pnr=3,
                                region_active_discount=1):
    # iterative_update_V
    """
    Update temporal components V_TF_ iteratively
    V_i = argmin ||Y-UV||_2^2 + sum_i lambda_i ||D^2V_i||_1
    (i.e. subtract off other components from Y, project
    normalized U_i onto residual and denoise V with TF)

    Parameters:
    ----------
    Y:              np.array (d x T)
                    2D video array (pixels x Time)
    U_hat:          np.array (d x k)
                    spatial components of Y
                    k is the estimated rank of Y
    V_TF_:          np.array (k x T)
                    temporal components of Y
                    k is the estimated rank of Y
    lambdas_:       np.array (k x 1)
                    lagrange multipliers to enforce ||D^2V_i||_1
                    where i corresponds to a single pixel
                    and D^2 is the second difference operator
                    if None: lambdas_s in initialized by recalculating
                    the noise if the temporal component
    verbose:        boolean
                    flag to enable verbose
    plot_en:        string
                    flag to enable plots

    Outputs:
    -------
    V_TF_2:         np.array (d x T)
                    updated V_TF_

    """
    num_components, T = V_TF_2.shape
    norm_U = (U_hat_**2).sum(axis=0)

    if plot_en:
        V_TF_=V_TF_2.copy()

    for ii in range(num_components):
        idx_ = np.setdiff1d(np.arange(num_components),ii)
        R_ = Y - U_hat_[:,idx_].dot(V_TF_2[idx_,:])
        V_ = U_hat_[:,ii].T.dot(R_)/norm_U[ii]

        if lambdas_ is None:
            if constraint_segmented:
                V_2,_, _ = trefide.constrained_l1tf(V_,
                                                    solver=solver,
                                                    region_thresh_min_pnr=region_thresh_min_pnr,
                                                    region_active_discount=region_active_discount,
                                                    verbose=verbose,
                                                    lagrange_scaled=lagrange_scaled)
            else:
                # deprecated
                noise_std_ = noise_estimator.get_noise_fft(V_[np.newaxis,:])[0]
                if fudge_factor is not None:
                    noise_std_ *=fudge_factor

                print('Solve V_[%d])'%(ii)) if verbose else 0
                print('Noise range is %.3e %.3e'%
                        (min(noise_std_),max(noise_std_))) if verbose else 0

                V_2 = c_l1tf_v_hat(V_, noise_std_)[0]

        else:
            if isinstance(lambdas_,int) or isinstance(lambdas_,float):
                V_2 = c_update_V(V_, lambdas_/norm_U[ii])
            else:
                if region_indices is None:
                    V_2 = c_update_V(V_, lambdas_[ii]/norm_U[ii])
                else:
                    V_2 = c_update_V(V_, lambdas_[ii]/norm_U[ii], region_indices[ii])

        V_TF_2[ii,:] = V_2.copy()

    uplot.plot_temporal_traces(V_TF_2,V_TF_) if plot_en else 0
    return V_TF_2


def iteration_error(Y,
                    U_hat,
                    V_TF,
                    region_indices=None,
                    lambdas_=None,
                    nus_=None,
                    U_update=False):
    """
    F(U,V)=||Y-UV||^2_2 + sum_i lambda_i ||D^2 V_i||_1 + sum_j nu_j ||U_j||_1
    # due to normalization F(U,V) may not decrease monotonically. problem?
    """

    T = Y.shape[1]
    F_uv1 = np.linalg.norm(Y - U_hat.dot(V_TF))**2

    if region_indices is None:
        diff = trefide.difference_operator(T)
        F_uv2  = np.sum(lambdas_*np.sum(np.abs(diff.dot(V_TF.T)),0),axis=0)
    else:
        F_uv2 = 0
        for ii_, region_ in enumerate(region_indices):
            for kk, reg_ in enumerate(region_):
                v_ = V_TF[ii_,reg_.flatten()]
                len_ = len(reg_.flatten())
                c_diff = trefide.difference_operator(len_)
                cte2 = (lambdas_[ii_][kk])*np.sum(np.abs(c_diff.dot(v_)),0)
                F_uv2  = F_uv2 + cte2

    F_uv3  = np.sum(nus_*np.sum(np.abs(U_hat),1)) if U_update else 0

    return F_uv1 , F_uv2, F_uv3


def greedy_component_denoiser(Y,
                              U_hat,
                              V_hat,
                              dims=None,
                              fudge_factor=1,
                              maxlag=5,
                              confidence=0.999,
                              corr=True,
                              mean_th=None,
                              kurto=False,
                              verbose=False,
                              plot_en=False,
                              U_update=False,
                              pca_method='vanilla',
                              final_regression=True,
                              max_num_iters=20,
                              max_num_components=20,
                              solver='SCS',
                              constraint_segmented=False, # call trefide
                              ignore_frames=None,
                              stim_knots=None,
                              stim_delta=10):
    """
    Denoise spatial and temporal components greedily
    F(U,V) = ||Y-UV||^2_2 + sum_i lambda_i ||D^2V_i||_1 + sum_j nu_j ||U_j||_1
    applying a soft constraint that the L2 norms of U and V are both 1
    i,j index components and pixels respectively
    lambda_i and nu_j are the lagrange multipliers

    Parameters:
    ----------

    Y:          np.array (d x T)
                2D video array (pixels x Time)
    U_hat:      np.array (d x k)
                spatial components of Y
                k is the estimated rank of Y
    V_hat:      np.array (k x T)
                temporal components of Y
                k is the estimated rank of Y
    dims:       tuple (d1 x d2 x T)
                dimensions of video array used for plotting
    fudge_factor: float
                constant to scale estimated noise std st denoising st denoising is less
                (lower factor) or more (higher factor) restrictive.
    maxlag:     int
                max correlation lag for correlation null hypothesis in samples
                (e.g. indicator decay in samples)
                employed to look for additional components in residual
    confidence: float
                confidence interval (CI) for correlation null hypothesis
                employed to look for additional components in residual
    corr:       boolean
                flag to include components which pass correlation null hypothesis
    kurto:      boolean
                flag to include components which pass kurtosis null hypothesis
    mean_th:    float
                threshold employed to reject components according to correlation null hypothesis
    U_update:   boolean
                flag to (not) update spatial components by imposing L1-constraint.
                True for "round" neurons in 2p.
                False for dendritic data.
    plot_en:    boolean
                flag to enable plots
    verbose:    boolean
                flag to enable verbose
    pca_method: string
                method for matrix decomposition (e.g. PCA, sPCA, rPCA, etc).
                see compute_svd for options

    Outputs:
    -------

    U_hat:      np.array (d x k1)
                spatial components of Y
                k1 is the updated estimated rank of Y
                k1 >= k depending on additional structured components added from residual.
    V_TF:       np.array (k1 x T)
                temporal components of Y
                k1 is the updated estimated rank of Y
                k1 >= k depending on additional structured components added from residual.
    """

    #mean_th = mean_th+0.05
    num_components, T = V_hat.shape

    uplot.plot_temporal_traces(V_hat) if plot_en else 0
    uplot.plot_spatial_component(U_hat,dims) if plot_en and (not dims==None) else 0

    #####################################################
    # Running iterations
    print('Initialization with %d components'%(num_components)) if verbose else 0
    print('solve V(i) = argmin_W ||D^2 W||_1 \n'
          +'\t st ||V_i-W||_2<'
          +'fudge_factor*sigma_i*sqrt(T)')if verbose else 0
    print('Max # of greedy loops: %d (relative convergence)'%max_num_iters) if verbose else 0

    rerun_1 = 1 # flag to run part (1)
    iteration = 0 # iteration # for part(1)

    if U_update:
        print('solve U(j) = argmin_W ||W||_1 '
              + 'st ||Y_j-W\'V_TF(j)||_2^2<T*fudge^2') if verbose else 0
    else:
        print('U = Y*pinv(V)') if verbose else 0

    while rerun_1:
        num_components = V_hat.shape[0]
        print('*Iteration %d: initialization %d components'
                %(iteration, num_components)) if verbose else 0

        print('solve V(i)') if verbose else 0

        V_TF,region_indices,lambdas_= greedy_temporal_denoiser(V_hat,
                                                               solver=solver,
                                                               constraint_segmented=constraint_segmented,
                                                               fudge_factor=fudge_factor,
                                                               verbose=verbose
                                                              )
        uplot.plot_temporal_traces(V_TF,V_hat) if plot_en else 0

        U_hat, nus_ = greedy_spatial_denoiser(Y,
                                              V_TF,
                                              fudge_factor=fudge_factor,
                                              U_update=U_update,
                                              dims=dims,
                                              plot_en=plot_en)

        lambdas_ = lambdas_ * (U_hat**2).sum(axis=0)

        uplot.plot_spatial_component(U_hat,dims) if plot_en and (not dims==None) else 0

        print('Iteration %d: begin greedy loops'%(iteration)) if verbose else 0
        F_UVs = []

        for loop_ in range(max_num_iters):
            print('\tIteration %d loop %d with %d components'%(
                            iteration,loop_,num_components)) if verbose else 0

            print('\tupdate V_i in closed form') if verbose else 0
            V_TF = iterative_temporal_denoiser(Y,
                                               U_hat,
                                               V_TF,
                                               solver=solver,
                                               lambdas_=lambdas_,
                                               region_indices=region_indices,
                                               plot_en=plot_en,
                                               verbose=verbose,
                                               constraint_segmented=constraint_segmented,
                                               fudge_factor=fudge_factor)

            print('\tupdate U_j in closed form/regression') if verbose else 0

            U_hat, _ = greedy_spatial_denoiser(Y,
                                               V_TF,
                                               nus_=nus_,
                                               fudge_factor=fudge_factor,
                                               U_update=U_update,
                                               dims=dims,
                                               plot_en=plot_en)

            F_uv1, F_uv2, F_uv3 = iteration_error(Y,
                                                  U_hat,
                                                  V_TF,
                                                  region_indices=region_indices,
                                                  lambdas_=lambdas_,
                                                  nus_=nus_,
                                                  U_update=U_update)

            F_uv = F_uv1 + F_uv2 + F_uv3
            F_UVs.append(F_uv)
            print('\tIteration %d loop %d error (%d+%d+%d)=%d'%(
                                    iteration,loop_,
                                    F_uv1,F_uv2,F_uv3,F_uv)) if verbose else 0

            if loop_ >=1:
                no_change = np.isclose(F_uv,F_UVs[loop_-1],rtol=1e-03, atol=1e-04)
                bad_iter = F_uv > F_UVs[loop_-1]
                if no_change or bad_iter:
                    print('\tIteration %d loop %d end - no significant updates'%(
                                            iteration,loop_)) if verbose else 0
                    break

        if plot_en:
            plt.title('Error F(u,v)')
            plt.plot(F_UVs,'o-')
            plt.show()

        print('*Iteration %d residual search with %d components'%(
                                iteration, V_TF.shape[0])) if verbose else 0

        ### (2) Compute PCA on residual R  and check for correlated components
        U_r, s_r, Vt_r = compute_svd((Y-U_hat.dot(V_TF)).astype('float32'),
                                     method=pca_method)

        ignore_frames =stimulus_segmentation(Vt_r.shape[1],
                                           stim_knots=stim_knots,
                                           stim_delta=stim_delta)

        ctid = choose_rank(Vt_r[:,~ignore_frames],
                           maxlag=maxlag,
                           confidence=confidence,
                           corr=corr,
                           kurto=kurto,
                           mean_th=mean_th)

        keep1_r = np.where(np.logical_or(ctid[0, :] == 1, ctid[1, :] == 1))[0]
        uplot.plot_vt_cov(Vt_r,keep1_r,maxlag) if plot_en else 0

        if len(keep1_r)==0:
            print('Final number of components %d'%V_TF.shape[0]) if verbose else 0
            rerun_1 = 0
        else:
            print('Rerun Iterations - adding %d components'%(len(keep1_r))) if verbose else 0
            num_components = num_components + len(keep1_r)
            if max_num_components <= num_components:
                print('Number of components %d > max allowed %d'%(num_components,max_num_components))
                rerun_1 = 0
            else:
                rerun_1 = 1
                iteration +=1
                V_hat = np.vstack((V_TF, Vt_r[keep1_r,:]))
                U_hat = np.hstack((U_hat,U_r[:,keep1_r].dot(np.diag(s_r[keep1_r]))))


    ##################
    ### Final update
    ##################

    print('Final update after %d iterations %d'%iteration) if verbose else 0
    print('\tFinal update of temporal components') if verbose else 0
    V_TF = iterative_temporal_denoiser(Y,
                                       U_hat,
                                       V_TF,
                                       solver=solver,
                                       plot_en=plot_en,
                                       constraint_segmented=constraint_segmented,
                                       fudge_factor=fudge_factor)

    print('\tFinal update of spatial components') if verbose else 0

    U_hat, _ = greedy_spatial_denoiser(Y,
                                       V_TF,
                                       fudge_factor=fudge_factor,
                                       U_update=U_update,
                                       dims=dims,
                                       plot_en=plot_en)

    if final_regression:
        print('\tFinal regression for V(j)') if verbose else 0
        V_TF_i = np.matmul(np.matmul(np.linalg.inv(np.matmul(U_hat.T, U_hat)), U_hat.T), Y)
    else:
        V_TF_i = V_TF.copy()
    uplot.plot_spatial_component(U_hat,dims) if plot_en and (not dims==None) else 0
    uplot.plot_temporal_traces(V_TF_i,V_TF) if plot_en else 0
    # this needs to be updated to reflect any new rank due to new numb of iterations
    return U_hat , V_TF_i


def denoise_components(data_all,
                       dims=None,
                       maxlag=5,
                       tsub=1,
                       ds=1,
                       noise_norm=False,
                       iterate=False,
                       confidence=0.99,
                       corr=True,
                       kurto=False,
                       tfilt=False,
                       tfide=False,
                       mean_th=None,
                       greedy=True,
                       mean_th_factor=1.,
                       p=1.,
                       fudge_factor=1.,
                       plot_en=False,
                       verbose=False,
                       U_update=False,
                       min_rank=0,
                       pca_method='vanilla',
                       detrend=False,
                       stim_knots=None,
                       stim_delta=200):
    """
    Compress array data_all as determined by parameters.

    Parameters:
    ----------

    data_all:   np.array (d x T) or (d1 x d2 xT)
                2D or 3D video array (pixels x Time) or (pixel x pixel x Time)
    dims:       tuple (d1 x d2 x T)
                dimensions of video array used for plotting
    maxlag:     int
                max correlation lag for correlation null hypothesis in samples
                (e.g. indicator decay in samples)
    tsub:       int
                temporal downsample constant
    ds:         int
                spatial downsample constant
    noise_norm: placeholder
    iterate:    boolean
                flag to include correlated components iteratively
    confidence: float
                confidence interval (CI) for correlation null hypothesis
    corr:       boolean
                flag to include components which pass correlation null hypothesis
    kurto:      boolean
                flag to include components which pass kurtosis null hypothesis
    tfilt:      boolean
                flag to temporally filter traces with AR estimate of order p.
    p:          int
                order of AR estimate, used only if tfilt is True
    tfide:      boolean
                flag to denoise temporal traces with Trend Filtering
    mean_th:    float
                threshold employed to reject components according to correlation null hypothesis
    min_rank:   int
                minimum rank of denoised/compressed matrix
                typically set to 1 to avoid empty output (array of zeros)
                if input array is mostly noise.
    greedy:     boolean
                flag to greedily update spatial and temporal components (estimated with PCA)
                greedyly by denoising temporal and spatial components
    mean_th_factor: float
                factor to scale mean_th
                typically set to 2 if greedy=True and mean_th=None or if mean_th has not been scaled yet.
    fudge_factor: float
                constant to scale estimated noise std st denoising st denoising is less
                (lower factor) or more (higher factor) restrictive.
    U_update:   boolean
                flag to (not) update spatial components by imposing L1- constraint.
                True for "round" neurons in 2p.
                False for dendritic data.
    plot_en:    boolean
                flag to enable plots
    verbose:    boolean
                flag to enable verbose
    pca_method: string
                method for matrix decomposition (e.g. PCA, sPCA, rPCA, etc).
                see compute_svd for options
    Outputs:
    -------

    Yd_out:     np.array (d x T)
                compressed/denoised array (dxT)
    ctids:      np.array (3,d)
                indicator 3D matrix (corr-kurto-reject) which points which statistic
                a given component passed and thus it is included.
                If greedy=True, all components added are included as corr components.
    """
    if data_all.ndim == 3:
        dims = data_all.shape
        data_all = data_all.reshape((np.prod(dims[:2]),dims[2]), order='F')
    data_all = data_all.T.astype('float32')
    # In a 2d matrix, we get rid of any broke (inf) pixels
    # we assume fixed broken across time
    broken_idx = np.isinf(data_all[0,:])
    # Work only on good pixels
    if np.any(broken_idx):
        print('broken pixels') if verbose else 0
        data = data_all[:, ~broken_idx]
    else:
        data = data_all.copy()

    # Remove the mean
    mu = data.mean(0, keepdims=True)
    data = data - mu

    # temporally filter the data
    if tfilt:
        print('Apply exponential filter') if verbose else 0
        data0 = np.zeros(data.shape)
        #T, num_pxls = data.shape
        for ii, trace in enumerate(data.T):
            # Estimate tau for exponential
            tau = cnmf.deconvolution.estimate_time_constant(
                    trace,p=p,sn=None,lags=5,fudge_factor=1.)
            window = tau **range(0,100)
            data0[:,ii] = np.convolve(fluor,window,mode='full')[:T]/np.sum(window)
    else:
        data0 = data.copy()

    # temporally decimate the data
    if tsub > 1:
        print('Temporal decimation %d'%tsub) if verbose else 0
        data0 = temporal_decimation(data0, tsub)

    # spatially decimate the data
    if ds > 1:
        print('Spatial decimation %d'%ds) if verbose else 0
        #D = len(dims)
        #ds = np.ones(D-1).astype('uint8')
        #data0 = spatial_decimation(data0, ds, dims)
        data0 = data0.copy()
    # Run svd

    start = time.time()
    #print(data0.shape)
    #print('Begin svd')
    U, s, Vt = compute_svd(data0.T,
                            method=pca_method)
    #print('svd run for %f'%(time.time()-start))

    # Project back if temporally filtered or downsampled
    if tfilt or tsub > 1:
        Vt = U.T.dot(data.T)

    # if greedy Force x2 mean_th (store only big components)
    if greedy and (mean_th_factor <= 1.):
        mean_th_factor = 2.

    # Select components
    #start = time.time()

    if mean_th is None:
        mean_th = wnoise_acov_CI(Vt.shape[1],
                                 confidence=confidence,
                                 maxlag=maxlag)
        mean_th *= mean_th_factor

    ignore_segments =stimulus_segmentation(Vt.shape[1],
                                           stim_knots=stim_knots,
                                           stim_delta=stim_delta
                                          )
    ctid = choose_rank(Vt[:,~ignore_segments],
                       maxlag=maxlag,
                       iterate=iterate,
                       confidence=confidence,
                       corr=corr,
                       kurto=kurto,
                       mean_th=mean_th)

    keep1 = np.where(np.logical_or(ctid[0, :] == 1, ctid[1, :] == 1))[0]

    #print('select run for %f'%(time.time()-start))

    # Plot temporal correlations
    uplot.plot_vt_cov(Vt,keep1,maxlag) if plot_en else 0

    # If no components to store, return min rank
    if np.all(keep1 == np.nan):
        if min_rank == 0:
            Yd = np.zeros(data.T.shape)
        else:
            min_rank=min_rank+1
            print('Forcing %d component(s)'%min_rank) if verbose else 0
            ctid[0,:min_rank]=1
            U = U[:,:min_rank].dot(np.eye(min_rank)*s[:min_rank])
            Yd = (U.dot(Vt[:min_rank,:]))
        Yd += mu.T
        return Yd, ctid

    Vt = Vt[keep1,:]
    #print('here 1191')

    # Denoise each temporal component
    if tfide:
        noise_levels = noise_estimator.get_noise_fft(Vt)[0]
        Vt = trefide.denoise(Vt, stdvs = noise_levels)

    if tfide and (tfilt or tsub > 1):
        U = data.T.dot(np.linalg.pinv(Vt).T)
    else:
        U = U[:,keep1].dot(np.eye(len(keep1))*s[keep1.astype(int)])

    # call greedy
    if greedy:
        try:
            U, Vt = greedy_component_denoiser(data.T,
                                              U,
                                              Vt,
                                              dims=dims,
                                              fudge_factor=fudge_factor,
                                              maxlag=maxlag,
                                              confidence=confidence,
                                              corr=corr,
                                              kurto=kurto,
                                              mean_th=mean_th/mean_th_factor,
                                              U_update=U_update,
                                              plot_en=plot_en,
                                              verbose=verbose,
                                              pca_method=pca_method,
                                              stim_knots=stim_knots,
                                              stim_delta=stim_delta)
            ctid[0,np.arange(Vt.shape[0])]=1
        except:
            print('\tERROR: Greedy solving failed, keeping %d parameters'%
                    (len(keep1)))
            ctid[0,0] = 100

    Yd = U.dot(Vt) + mu.T

    if broken_idx.sum() > 0:
        print('ERROR: There are %d broken pixels.'%(broken_idx.sum()))
        Yd_out = np.ones(shape=data_all.shape).T*np.inf
        Yd_out[~broken_idx,:] = Yd
    else:
        Yd_out =  Yd

    return Yd_out, ctid


#########################################################################################
## DURING MERGING PHASE
#########################################################################################


def c_l1tf_v_hat(v,
                 sigma,
                 verbose=False,
                 solver='SCS',
                 max_iters=1000):
    """
    Update vector v according to difference fctn diff
    with noise_std(v) = sigma

    V(i) = argmin_W ||D^2 W||_1
    st ||V_i-W||_2<sigma_i*sqrt(T)
    Include optimal lagrande multiplier for constraint

    """
    print('c_l1tf_v_hat') if verbose else 0
    if np.abs(sigma)<=1e-3:
        print('Do not denoise (high SNR: noise_level=%.3e)'%
                sigma) if verbose else 0
        return v , 0
    T = len(v)
    v_hat = cp.Variable(T)
    print(sigma*np.sqrt(T)) if verbose else 0
    diff = trefide.difference_operator(T)

    objective = cp.Minimize(cp.norm(diff*v_hat,1))
    constraints = [cp.norm(v-v_hat,2)<=sigma*np.sqrt(T)]

    cp.Problem(objective, constraints).solve(solver=solver,
                                             max_iters=max_iters,
                                             verbose=False)
    lambda_= constraints[0].dual_value
    if lambda_ !=0:
        lambda_=1./lambda_

    return np.asarray(v_hat.value).flatten(), lambda_


def c_l1_u_hat(y,
               V_TF,
               fudge_factor=1,
               verbose=False,
              solver='SCS',
              max_iters=1000):
    """
    update array U given Y and V_TF

    U(j) = argmin_W ||W||_1
    st ||Y_j-W'V_TF(j)||_2^2 < T
    if problem infeasible:
        set U = regression Vt onto Y and \nu = 0
    """
    print('c_l1_u_hat') if verbose else 0
    num_components = V_TF.shape[0]
    u_hat = cp.Variable(num_components)
    objective = cp.Minimize(cp.norm(u_hat,1))
    constraints = [cp.norm(y[np.newaxis,:] - u_hat.T*V_TF,2) < np.sqrt(len(y))*fudge_factor]
    problem = cp.Problem(objective,constraints)
    problem.solve(solver=solver,
                  max_iters=max_iters,
                  verbose=False)
    if problem.status in ["infeasible", "unbounded"]:
        return y[np.newaxis,:].dot(np.linalg.pinv(V_TF)).flatten(), 0
    else:
        return np.asarray(u_hat.value).flatten(), 1/constraints[0].dual_value


def l1tf_dual(v,
              lambda_=1.0,
               solver='cython',
               cvxpy_solver='SCS',
               max_iters=400,
               verbose=False,
               max_lambda=False,
               rescale=1
                  ):
    """
    Solve with either of these solvers

    """
    if solver == 'cython':
        x_hat = l1_tf(y,
                    lambda_,
                    max_lambda,
                    max_iters,
                    1)

    elif solver =='cvxpy':
        T = len(v)
        diff = trefide.difference_operator(T)
        print(lambda_) if verbose else 0
        v_hat = cp.Variable(T)
        cte2 = lambda_*cp.norm(diff*v_hat,1)
        objective = cp.Minimize(cp.norm(v-v_hat,2)**2 + cte2)

        cp.Problem(objective).solve(solver=cvxpy_solver,
                                    max_iters=max_iters,
                                    verbose=False)

        v_hat = np.asarray(v_hat.value).flatten()

    elif solver == 'g_or':
        pass

    return v_hat


def c_update_V(v,
                lambda_,
                region_indices=None,
                verbose=False,
                plot_en=False,
                solver='SCS',
                max_iters=1000,
               max_lambda=False):
    """
    merging in progress
    min ||Y-UV||_2^2 + sum_i lambda_i||D^2V_i||_1
    # Fixing U we have
    min ||v-v_hat||_2^2 + lambda_i||D^2V_i||_1
    """
    T = len(v)
    diff = trefide.difference_operator(T)
    print(lambda_) if verbose else 0

    if isinstance(lambda_,int) or isinstance(lambda_,float):
        print('1255') if verbose else 0
        v_hat = l1tf_dual(v,
                      lambda_=lambda_,
                      cvxpy_solver=solver,
                      max_iters=max_iters,
                      max_lambda=max_lambda,
                      verbose=verbose)
        return v_hat

    if region_indices is None:
        # deprecated- editing merge
        print('1267') if verbose else 0
        v_hat = l1tf_dual(v,
                      lambda_=lambda_,
                      cvxpy_solver=solver,
                      max_iters=max_iters,
                      max_lambda=max_lambda,
                      verbose=verbose)
        return v_final

    print('1092') if verbose else 0
    v_final = np.zeros(len(v))
    for ii, region_ in enumerate(region_indices):
        print('Denoising component %d'%ii) if verbose else 0
        v_ = v[region_.flatten()]
        v_hat = l1tf_dual(v_,
              lambda_=lambda_[ii],
              cvxpy_solver=solver,
              max_iters=max_iters,
              max_lambda=max_lambda,
              verbose=verbose)

        v_final[region_.flatten()]=v_hat
        uplot.plot_temporal_traces(v_hat[:,np.newaxis].T,
                                   v_[:,np.newaxis].T) if plot_en else 0
    return v_final


def c_update_V_deprecated(v,
                lambda_,
                region_indices=None,
                verbose=False,
                plot_en=False,
                solver='SCS',
                max_iters=1000):
    """
    DEPRECATED ---
    min ||Y-UV||_2^2 + sum_i lambda_i||D^2V_i||_1
    # Fixing U we have
    min ||v-v_hat||_2^2 + lambda_i||D^2V_i||_1
    """
    T = len(v)
    diff = trefide.difference_operator(T)
    print(lambda_) if verbose else 0

    if isinstance(lambda_,int) or isinstance(lambda_,float):
        print('1255') if verbose else 0
        v_hat = cp.Variable(T)
        cte2 = lambda_*cp.norm(diff*v_hat,1)

        objective = cp.Minimize(cp.norm(v-v_hat,2)**2 + cte2)

        cp.Problem(objective).solve(solver=solver,
                                    max_iters=max_iters,
                                    verbose=False)
        return np.asarray(v_hat.value).flatten()

    if region_indices is None:
        print('1267') if verbose else 0
        v_hat = cp.Variable(len(v))
        cte2 = lambda_*cp.norm(diff*v_hat,1)
        objective = cp.Minimize(
                                cp.norm(v-v_hat,2)**2
                                + cte2)

        cp.Problem(objective).solve(solver=solver,
                                    max_iters=max_iters,
                                    verbose=False)
        v_final = np.asarray(v_hat.value).flatten()
        return v_final

    print('1092') if verbose else 0
    v_final = np.zeros(len(v))
    for ii, region_ in enumerate(region_indices):
        print('Denoising component %d'%ii) if verbose else 0
        v_ = v[region_.flatten()]
        len_ = len(v_)
        v_hat = cp.Variable(len_)
        c_diff = (np.diag(2*np.ones(len_),0)+
                np.diag(-1*np.ones(len_-1),1)+
                np.diag(-1*np.ones(len_-1),-1)
                )[1:len_]
        cte2 = lambda_[ii]*cp.norm(c_diff*v_hat,1)
        objective = cp.Minimize(cp.norm(v_-v_hat,2)**2 + cte2)
        cp.Problem(objective).solve(solver=solver,
                                    max_iters=max_iters,
                                    verbose=False)
        v_final[region_.flatten()]=np.asarray(v_hat.value).flatten()
        uplot.plot_temporal_traces(np.asarray(v_hat.value).flatten()
                [:,np.newaxis].T,v_[:,np.newaxis].T) if plot_en else 0
    return v_final


def c_update_U(y,
               V_TF,
               nu_,
              solver='ECOS',
              max_iters=1000):
    """
    min ||Y-UV||^2_2 + sum_j nu_j ||U_j||_1.
    for each pixel
    min  ||y_j-u_j*v||_2^2 + nu_j ||u_j||_1.
    """
    num_components = V_TF.shape[0]
    u_hat = cp.Variable(num_components)
    if nu_ == 0:
        objective = cp.Minimize(
                cp.norm(y[np.newaxis,:]-u_hat.T*V_TF,2)**2)
    else:
        objective = cp.Minimize(
                            cp.norm(
                                    y[np.newaxis,:]-u_hat.T*V_TF
                                ,2)**2
                                + nu_*cp.norm(u_hat,1)
                                )
    problem = cp.Problem(objective)
    problem.solve(solver=solver,
                    verbose=False,
                    max_iters=max_iters)
    return np.asarray(u_hat.value).flatten()


#### DEPRECATED parallel implementation

def update_u_parallel(Y,V_TF,fudge_factor):
    pool = multiprocessing.Pool(processes=20)
    c_outs = pool.starmap(c_l1_u_hat, itertools.product(y, V_TF, fudge_factor))
    pool.close()
    pool.join()
    return c_outs


def c_update_U_parallel(Y,V_TF,nus_):
    """
    call c_update_U as queue
    """
    pool = multiprocessing.Pool()#processes=20)
    args = [(y,V_TF,nus_[idx]) for idx, y in enumerate(Y)]
    c_outs = pool.starmap(c_update_U, args)

    pool.close()
    pool.join()
    return c_outs
