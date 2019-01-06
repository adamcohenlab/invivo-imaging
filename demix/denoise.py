import numpy as np


import noise_estimator
import spatial_filtering
import tool_grid
import time

#from trefide.utils.noise import estimate_noise


# ____________________________
# Wrapper to call denoisers
# see individual functions for more information
# ____________________________


def spatial(Y_new,
            gHalf=[2,2],
            sn=None):
    """
    Calls spatial wiener filter in pixel neighborhood
    """
    mov_wf,_ = spatial_filtering.spatial_filter_image(Y_new,
                                                gHalf=gHalf,
                                                sn=sn)
    return mov_wf


def temporal(W,
             confidence=0.999,
             dx=1,
             fudge_factor=1,
             greedy=False,
             maxlag=3,
             mean_th_factor=1.15,
             min_rank=1,
             nblocks=[10,10],
             snr_threshold=2,
             U_update=False,
             verbose=False):
    """
    Calls greedy temporal denoiser in pixel neighborhood
    """
    #start = time.time()
    mov_d, ranks = tool_grid.denoise_dx_tiles(W,
                                              confidence=confidence,
                                              dx=dx,
                                              fudge_factor=fudge_factor,
                                              greedy=greedy,
                                              maxlag=maxlag,
                                              mean_th_factor=mean_th_factor,
                                              min_rank=min_rank,
                                              nblocks=nblocks,
                                              snr_threshold=snr_threshold,
                                              U_update=U_update,
                                              verbose=verbose)
    #print('Temporal denoiser run for %.3f sec'%(time.time()-start))
    return mov_d, ranks


def noise_level(mov_wf,
                range_ff =[0.25,0.5],
               noise_method='mean'):
    """
    Calculate noise level in movie pixels
    use mean instead of logmexp
    """
    ndim_ = np.ndim(mov_wf)
    if ndim_==3:
        dims_ = mov_wf.shape
        mov_wf = mov_wf.reshape((np.prod(dims_[:2]), dims_[2]),order='F')
    noise_level = noise_estimator.get_noise_fft(mov_wf,
                                                      noise_range=range_ff)[0]
    #noise_level = noise_estimator.estimate_noise(mov_wf,
    #                                             summarize=noise_method)#[0]
    if ndim_ ==3:
        noise_level = noise_level.reshape(dims_[:2], order='F')

    return noise_level