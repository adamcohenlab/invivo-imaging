# invivo-imaging

Spike-Guided Penalized Matrix Decomposition-Non-negative Matrix Factorization (SGPMD-NMF) pipeline code for in vivo voltage imaging and used for data in [1]

1.	NoRMCorre correction of x-y motion [2]
2.	Trimming edges of movie to remove pixels that entered or exited the field of view
3.	Photobleach correction with b-spline fit
4.	PMD denoiser with parameters optimized on simulated data [3]
5.	From each pixel, regress out residual motion artifacts using NoRMCorre-estimated motion trajectories   
6.	Manually crop blood vessels
7.	localNMF demixing on temporally high-pass filtered movie to obtain spatial support of cells
8.  fastHALS iterations on full unfiltered movie to calculate full spatial footprints of cells and background
9.	Update spatial footprint of cells and background by ensuring smoothness of background spatial components around edges of cells
10.	Apply the updated spatial and background footprints to calculate the temporal traces from the full denoised movie

## Instructions

Detailed instructions for use on Harvard's Cannon cluster: https://docs.google.com/document/d/1dV4yBQQuKRmscEwlC_2DaNGKz9M_jCuN6qnOgfYe2ew/edit?usp=sharing

## Dependencies

* TREFIDE (http://github.com/ikinsella/trefide)

## References

[1] Adam, Y. et al. Voltage imaging and optogenetics reveal behavior dependent changes in hippocampal dynamics. Nature (2019)

[2] Pnevmatikakis, E. A. & Giovannucci, A. NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. J Neurosci Methods (2017). 

[3] Buchanan, E. K. et al. Penalized matrix decomposition for denoising, compression, and improved demixing of functional imaging data. bioRxiv, 334706 (2019). 
