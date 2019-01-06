# Cocaim
Compression and denoising of Calcium Imaging recordings

Example Demo_denoiser.ipynb <br />
Dataset from https://github.com/simonsfoundation/CaImAn

## Comments:
1. Runs in python3
2. We assume the data is motion corrected and detrended by your method of choice.
3. Use four greedy denoisers (M1 applied in 4 grids) if M1 alone results in several block artifacts.
4. Ongoing: NMF step for different datasets.

## Dependencies
Dependencies:
1. [CaImAn](https://github.com/simonsfoundation/CaImAn) - follow their setup instructions
2. [CVXPY](https://cvxgrp.github.io/cvxpy/install/index.html) - follow their setup instructions, however on Ubuntu and Py3.6 I found that I first needed to use conda to install blas and lapack as described in the last comment [here](https://github.com/cvxgrp/cvxpy/issues/357). 
3. [CVXOPT] (http://cvxopt.org/).
