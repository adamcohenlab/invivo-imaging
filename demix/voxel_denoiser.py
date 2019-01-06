from tool_grid import run_single,combine_blocks,combine_4xd,offset_tiling,combine_4xd
import numpy as np


def offset_tiling_voxel(W1,
                        nblocks=[10,10],
                        offset_case='rc'):

    arrays = []
    for W in W1:
        W_rs1, drs = offset_tiling(W,
                     nblocks=nblocks,
                     offset_case=offset_case)
        arrays.append(W_rs1)
    return arrays,drs



def interleave_3D(ps):
    cc1 = np.concatenate(ps,axis=2)
    cc = np.zeros(cc1.shape)
    cc[:,:,0::3]=ps[0]
    cc[:,:,1::3]=ps[1]
    cc[:,:,2::3]=ps[2]
    return cc

def deinterleave_3D(cc):
    ps=[]
    ps.append(cc[:,:,0::3])
    ps.append(cc[:,:,1::3])
    ps.append(cc[:,:,2::3])
    return ps


def denoise_dx_voxel_tiling(W1,W2,W3,
                             nblocks=[10,10],
                             dx=1,
                             maxlag=5,
                             confidence=0.99,
                             greedy=False,
                             fudge_factor=0.99,
                             mean_th_factor=1.15,
                             U_update=False,
                             min_rank=1,
                            interleave=False,
                             stim_knots=None,
                             stim_delta=200):


    """
    Given matrix W, denoise it according
    Input:
    ------

    Output:
    ------
    """
    dims = W1.shape

    
    dW_1, dW_2, dW_3, rank_W_, dims_ = denoise_voxel_tiling(W1,W2,W3,
                                                            nblocks=nblocks,
                                                            offset_case=None,
                                                             interleave=interleave,
                                                             maxlag=maxlag,
                                                             confidence=confidence,
                                                             greedy = greedy,
                                                             fudge_factor=fudge_factor,
                                                             mean_th_factor=mean_th_factor,
                                                             U_update=U_update,
                                                             min_rank=min_rank,
                                                             stim_knots=stim_knots,
                                                             stim_delta = stim_delta)

    if dx == 1:
        return dW_1, dW_2, dW_3, rank_W_

    dW_rs1, dW_rs2, dW_rs3, rank_W_rs, dims_rs = denoise_voxel_tiling(W1,W2,W3,
                                                                     nblocks=nblocks,
                                                                      offset_case='r',
                                                                     interleave=interleave,
                                                                     maxlag=maxlag,
                                                                     confidence=confidence,
                                                                     greedy = greedy,
                                                                     fudge_factor=fudge_factor,
                                                                     mean_th_factor=mean_th_factor,
                                                                     U_update=U_update,
                                                                     min_rank=min_rank,
                                                                     stim_knots=stim_knots,
                                                                     stim_delta = stim_delta)


    dW_cs1, dW_cs2, dW_cs3, rank_W_cs, dims_cs = denoise_voxel_tiling(W1,W2,W3,
                                                                     nblocks=nblocks,
                                                                      offset_case='c',
                                                                     interleave=interleave,
                                                                     maxlag=maxlag,
                                                                     confidence=confidence,
                                                                     greedy = greedy,
                                                                     fudge_factor=fudge_factor,
                                                                     mean_th_factor=mean_th_factor,
                                                                     U_update=U_update,
                                                                     min_rank=min_rank,
                                                                     stim_knots=stim_knots,
                                                                     stim_delta = stim_delta)


    dW_rcs1, dW_rcs2, dW_rcs3, rank_W_rcs, dims_rcs = denoise_voxel_tiling(W1,W2,W3,
                                                                     nblocks=nblocks,
                                                                        offset_case='rc',
                                                                     interleave=interleave,
                                                                     maxlag=maxlag,
                                                                     confidence=confidence,
                                                                     greedy = greedy,
                                                                     fudge_factor=fudge_factor,
                                                                     mean_th_factor=mean_th_factor,
                                                                     U_update=U_update,
                                                                     min_rank=min_rank,
                                                                     stim_knots=stim_knots,
                                                                     stim_delta = stim_delta)

    W_four_1 = combine_4xd(nblocks,
                         dW_1,
                         dW_rs1,
                         dW_cs1,
                         dW_rcs1,
                         dims_,
                         dims_rs,
                         dims_cs,
                         dims_rcs)
    
    W_four_2 = combine_4xd(nblocks,
                     dW_2,
                     dW_rs2,
                     dW_cs2,
                     dW_rcs2,
                     dims_,
                     dims_rs,
                     dims_cs,
                     dims_rcs)
    
    W_four_3 = combine_4xd(nblocks,
                     dW_3,
                     dW_rs3,
                     dW_cs3,
                     dW_rcs3,
                     dims_,
                     dims_rs,
                     dims_cs,
                     dims_rcs)

    return W_four_1, W_four_2, W_four_3, [rank_W_,rank_W_rs,rank_W_cs,rank_W_rcs]


def denoise_voxel_tiling(W1,W2,W3,
                         dims=None,
                         nblocks=[10,10],
                         offset_case=None,
                         interleave=False,
                         maxlag=5,
                         confidence=0.9,
                         greedy = False,
                         fudge_factor=1,
                         mean_th_factor=1.25,
                         U_update=False,
                         min_rank=1,
                         stim_knots=None,
                         stim_delta = 200,
                         reconstruct=True,
                        ):
    """

    """

    if not isinstance(W1,list):
        W1,dims = offset_tiling(W1,
                         nblocks=nblocks,
                         offset_case=offset_case)
        W2,_ = offset_tiling(W2,
                         nblocks=nblocks,
                         offset_case=offset_case)

        W3,_ = offset_tiling(W3,
                         nblocks=nblocks,
                         offset_case=offset_case)
        
    dims_ = np.asarray(list(map(np.shape,W1)))

 
    if interleave:
        args = [interleave_3D(ps) for ii, ps in enumerate(zip(W1,W2,W3))]
    else:
        args = [np.concatenate(ps,axis=0) for ii,ps in enumerate(zip(W1,W2,W3))]

    Yds, vtids = run_single(args,
                   maxlag=maxlag,
                   confidence=confidence,
                   greedy=greedy,
                   fudge_factor=fudge_factor,
                   mean_th_factor=mean_th_factor,
                   U_update=U_update,
                   min_rank=min_rank,
                   stim_knots=stim_knots,
                   stim_delta=stim_delta)


    if interleave:
        a = [deinterleave_3D(x) for x in Yds]
    else:
        #dims_ = np.asarray(list(map(np.shape,W1)))
        idx_ = np.arange(1,3)[:,np.newaxis].dot(dims_.T[0][np.newaxis,:]).T
        a = [np.split(x,c_idx,axis=0) for x ,c_idx in zip(Yds,idx_)]

    dW1_cs = [a_[0] for a_ in a]
    dW2_cs = [a_[1] for a_ in a]
    dW3_cs = [a_[2] for a_ in a]

    if reconstruct:
        dW1_cs = combine_blocks(dims,dW1_cs,list_order='C')
        dW2_cs = combine_blocks(dims,dW2_cs,list_order='C')
        dW3_cs = combine_blocks(dims,dW3_cs,list_order='C')

    return dW1_cs,dW2_cs,dW3_cs,vtids,dims_