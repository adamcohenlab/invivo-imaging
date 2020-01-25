# invivo-imaging

Spike-Guided Penalized Matrix Decomposition-Non-negative Matrix Factorization (SGPMD-NMF) pipeline code for in vivo voltage imaging

1.	NoRMCorre correction of x-y motion
2.	Photobleach correction with b-spline fit
3.	PMD denoiser with parameters optimized on simulated data
4.	From each pixel, regress out residual motion artifacts using NoRMCorre-estimated motion trajectories   
5.	Manually crop blood vessels
6.	localNMF demixing on temporally high-pass filtered movie to obtain spatial support of cells
7.  fastHALS iterations on full unfiltered movie to calculate full spatial footprints of cells and background
8.	Update spatial footprint of cells and background by ensuring smoothness of background spatial components around edges of cells
9.	Apply the updated spatial and background footprints to calculate the temporal traces from the full denoised movie

## Instructions

Instructions to run the pipeline here: http://bit.ly/sgpmdnmf-steps

Detailed instructions for setup on Harvard's Cannon cluster: http://bit.ly/harvard-sgpmdnmf-steps

## Dependencies

* [TREFIDE](http://github.com/ikinsella/trefide)

## References

[1] Xie, M., Adam, Y., Fan, L., Böhm, U., Kinsella, I., Zhou, D., Paninski, L., Cohen, A. High fidelity estimates of spikes and subthreshold waveforms from 1-photon voltage imaging in vivo. Submitted (2020)
