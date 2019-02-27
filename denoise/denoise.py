# General Dependencies
import timeit, os, sys
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

import time, math

data_dir = sys.argv[1]
out_dir = sys.argv[3]
mov_in = sys.argv[2]
detr_spacing = int(sys.argv[4])
rblocks = int(sys.argv[5])
cblocks = int(sys.argv[6])

if not os.path.isfile(out_dir + '/detr.tif'):
	if mov_in[-4:] == '.tif':
	# load movie from tif
		raw_mov = imio.imread(data_dir + '/' + mov_in).transpose(1,2,0)
		nrows = raw_mov.shape[0]
		ncols = raw_mov.shape[1]
	elif mov_in[-4:] == '.bin':
		raw_mov = np.fromfile(data_dir + '/' + mov_in,dtype=np.int16)
		with open(data_dir + '/experimental_parameters.txt') as file:
			x = file.read()
			params = [int(i) for i in re.findall(r'\d+',x)]
		nrows = params[0]
		ncols = params[1]
		raw_mov = np.reshape(raw_mov,(-1, nrows, ncols)).transpose(1,2,0)
	else
		raise ValueError('Only .tif and .bin files supported for input.')

	raw_mov = raw_mov[math.floor((nrows % (2 * rblocks))/2):-math.ceil((nrows % (2 * rblocks))/2),math.floor((ncols % (2 * cblocks))/2):-math.ceil((ncols % (2 * cblocks))/2),:]

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

	mov_detr, trend, stim, disc_idx = detrend(mov, stim, disc_idx.squeeze(), visualize=None, spacing=detr_spacing)

	print("Detrending took: " + str(time.time()-start) + ' sec\n')

	fov_height, fov_width, num_frames = mov_detr.shape

	mov_detr_nnorm = np.ascontiguousarray(mov_detr)
	Sn_image = np.sqrt(np.reshape(psd_noise_estimate(np.reshape(mov_detr_nnorm, 
													 (fov_height*fov_width, num_frames))),
						  (fov_height, fov_width, 1)))
	mov_detr_nnorm = mov_detr_nnorm / Sn_image

	imio.imsave(out_dir + "/detr.tif", mov_detr)
	imio.imsave(out_dir + "/detr_nnorm.tif", mov_detr_nnorm)
	imio.imsave(out_dir + "/Sn_image.tif", Sn_image)
	imio.imsave(out_dir + "/trend.tif", mov - mov_detr)
	print("Detrended movies saved\n")
else:
	mov_detr_nnorm = imio.imread(out_dir + '/detr_nnorm.tif')
	Sn_image = imio.imread(out_dir + "/Sn_image.tif")
	print("Detrended movies loaded\n")

if not os.path.isfile(out_dir + '/denoised.tif'):
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

	# blocks = io.loadmat(data_dir + '/num_blocks.mat')
	# blocks = blocks['var'][0]
	block_height = int(fov_height / rblocks)
	block_width = int(fov_width / cblocks)
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

	imio.imsave(out_dir + '/denoised.tif',mov_denoised)
	imio.imsave(out_dir + '/PMD_residual.tif',mov - mov_denoised)
	np.save(out_dir + '/block_ranks.npy', block_ranks)

print('Finished\n')
