# General Dependencies
import timeit, os
import numpy as np
import scipy.io as io

# Preprocessing Dependencies
from trefide.utils import psd_noise_estimate

# PMD Model Dependencies
from trefide.pmd import batch_decompose, batch_recompose, overlapping_batch_decompose, overlapping_batch_recompose

# Plot & Video Dependencies
import matplotlib.pyplot as plt
from trefide.plot import pixelwise_ranks
from trefide.extras.util_plot import comparison_plot, correlation_traces, snr_per_frame, nearest_frame_corr
from trefide.video import write_mpl, play_cv2
from trefide.preprocess import flag_outliers, detrend

#from trefide.decimation import decimated_decompose, decimated_batch_decompose, downsample_video, downsample_image, downsample_signal, overlapping_dbd

from skimage import io as imio

import time

if not os.path.isfile('./detr.tif'):
	# load movie from tif
	raw_mov = imio.imread('trimmed.tif').transpose(1,2,0)
	raw_stim = 10 * np.ones(raw_mov.shape[2]) # simulate stimulation values for in vivo data

	outliers = False # outliers already removed
	print('Trimmed movie loaded\n')

	if outliers:
		outlier_idx, disc_idx = flag_outliers(raw_mov[20,140], thresh_stdv=4, buffer=10, visualize=None)
		mov = np.delete(raw_mov, outlier_idx, axis=-1)
		stim = np.delete(raw_stim, outlier_idx)
	else:
		mov = raw_mov
		stim = raw_stim
		disc_idx = np.array([])    

	start = time.time()

	mov_detr, trend, stim, disc_idx = detrend(mov, stim, disc_idx.squeeze(), visualize=None, spacing=5000)

	print("Detrending took: " + str(time.time()-start) + ' sec\n')

	fov_height, fov_width, num_frames = mov_detr.shape

	mov_detr_nnorm = np.ascontiguousarray(mov_detr)
	Sn_image = np.sqrt(np.reshape(psd_noise_estimate(np.reshape(mov_detr_nnorm, 
													 (fov_height*fov_width, num_frames))),
						  (fov_height, fov_width, 1)))
	mov_detr_nnorm = mov_detr_nnorm / Sn_image

	imio.imsave("detr.tif", mov_detr)
	imio.imsave("detr_nnorm.tif", mov_detr_nnorm)
	imio.imsave("Sn_image.tif", Sn_image)
	imio.imsave("trend.tif", mov - mov_detr)
	print("Detrended movies saved\n")
else:
	mov_detr_nnorm = imio.imread('detr_nnorm.tif')
	Sn_image = imio.imread("Sn_image.tif")
	print("Detrended movies loaded\n")

if not os.path.isfile('./denoised.tif'):
	def compute_thresh(samples, conf=5):
		return np.percentile(samples, conf)
	def tv_norm(image):
		return np.sum (np.abs(image[:,:-1] - image[:,1:])) + np.sum(np.abs(image[:-1,:] - image[1:,:]))

	def spatial_test_statistic(component):
		d1, d2 = component.shape
		return (tv_norm(component) *d1*d2)/ (np.sum(np.abs(component)) * (d1*(d2-1) + d2 * (d1-1)))

	def temporal_test_statistic(signal):
		return np.sum(np.abs(signal[2:] + signal[:-2] - 2*signal[1:-1])) / np.sum(np.abs(signal))
	def determine_thresholds(mov_dims, block_dims, num_components, conf=5, plot=False):
		
		# Simulate Noise Movie
		noise_mov = np.ascontiguousarray(np.reshape(np.random.randn(np.prod(mov_dims)), mov_dims))
		
		# Perform Blockwise PMD Of Noise Matrix In Parallel
		spatial_components,\
		temporal_components,\
		block_ranks,\
		block_indices = batch_decompose(mov_dims[0], mov_dims[1], mov_dims[2],
										noise_mov, block_dims[0], block_dims[1],
										1e3, 1e3,
										num_components, consec_failures,
										max_iters_main, max_iters_init, tol,
										d_sub=d_sub, t_sub=t_sub)
		
		# Gather Test Statistics
		spatial_stat = []
		temporal_stat = []
		num_blocks = int((mov_dims[0] / block_dims[0]) * (mov_dims[1] / block_dims[1]))
		for block_idx in range(num_blocks): 
			for k in range(int(block_ranks[block_idx])):
				spatial_stat.append(spatial_test_statistic(spatial_components[block_idx,:,:,k]))
				temporal_stat.append(temporal_test_statistic(temporal_components[block_idx,k,:]))

		# Compute Thresholds
		spatial_thresh = compute_thresh(spatial_stat, conf=conf)
		temporal_thresh = compute_thresh(temporal_stat, conf=conf)
		
		if plot:
			fig, ax = plt.subplots(2,2,figsize=(8,8))
			ax[0,0].scatter(spatial_stat, temporal_stat, marker='x', c='r', alpha = .2)
			ax[0,0].axvline(spatial_thresh)
			ax[0,0].axhline(temporal_thresh)
			ax[0,1].hist(temporal_stat, bins=20, color='r')
			ax[0,1].axvline(temporal_thresh)
			ax[0,1].set_title("Temporal Threshold: {}".format(temporal_thresh))
			ax[1,0].hist(spatial_stat, bins=20, color='r')
			ax[1,0].axvline(spatial_thresh)
			ax[1,0].set_title("Spatial Threshold: {}".format(spatial_thresh))
			plt.show()
		
		return spatial_thresh, temporal_thresh
	tol = 5e-3

	mov = mov_detr_nnorm

	fov_height, fov_width, num_frames = mov.shape
	num_pixels = fov_height * fov_width

	blocks = io.loadmat('num_blocks.mat')
	blocks = blocks['var'][0]
	block_height = int(fov_height / blocks[0])
	block_width = int(fov_width / blocks[1])
	max_components = 50
	max_iters_main = 10
	max_iters_init = 40
	consec_failures = 3
	d_sub=1
	t_sub=1

	start = time.time()

	spatial_thresh, temporal_thresh = determine_thresholds((fov_height, fov_width, num_frames),
														   (block_height, block_width),
														   consec_failures, conf=5, plot=False)

	print("Simulation took: " + str(time.time()-start) + ' sec')

	start = time.time()

	# Perform 4x Overlapping Blockwise PMD In Parallel
	spatial_components,\
	temporal_components,\
	block_ranks,\
	block_indices,\
	block_weights = overlapping_batch_decompose(fov_height, fov_width, num_frames,
												mov, block_height, block_width,
												spatial_thresh, temporal_thresh,
												max_components, consec_failures,
												max_iters_main, max_iters_init, tol,
												d_sub=d_sub, t_sub=t_sub)

	# Use Compressed Components To Reconstruct Denoise Video
	mov_denoised = np.asarray(overlapping_batch_recompose(fov_height, fov_width, num_frames,
														  block_height, block_width,
														  spatial_components,
														  temporal_components,
														  block_ranks,
														  block_indices,
														  block_weights)) 

	print("Denoising took: " + str(time.time()-start) + ' sec')

	imio.imsave('denoised.tif',mov_denoised)
	imio.imsave('PMD_residual.tif',mov - mov_denoised)
	np.save('block_ranks.npy', block_ranks)
else:
	mov_denoised = imio.imread('denoised.tif')
	print('Denoised movie loaded\n')
if not os.path.isfile('./denoised_15s.tif'):
	imio.imsave('denoised_15s.tif',mov_denoised[:,:,:15000] * np.squeeze(np.repeat(np.expand_dims(Sn_image,2),15000,axis=2)))
	print('Denoised snippet saved\n')

print('Finished\n')
