# invivo-imaging

analysis pipeline code for in vivo voltage imaging, modified from [1]

1.	NoRMCorre correction of x-y motion [2]
2.	Trimming edges of movie to remove pixels that entered or exited the field of view
3.	Photobleach correction with b-spline fit
4.	PMD denoiser with parameters optimized on simulated data
5.	From each pixel, regress out residual motion artifacts using NoRMCorre-estimated motion trajectories   
6.	Manually crop blood vessels
7.	Select areas of background to initialize background estimator
8.	Superpixel-based NMF demixing to calculate spatial footprints of cells and background
9.	Regress the spatial footprint of cells out of the background (using a constant intercept term) to update spatial footprint of cells and background
10.	Apply the updated spatial and background footprints to calculate the temporal traces from the full denoised movie

References:

[1] Buchanan, E. K. et al. Penalized matrix decomposition for denoising, compression, and improved demixing of functional imaging data. bioRxiv, 334706 (2018). 

[2] Pnevmatikakis, E. A. & Giovannucci, A. NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. bioRxiv, 108514 (2017). 
