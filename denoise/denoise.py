# General Dependencies
import timeit, os, sys, re
import numpy as np
import scipy
import trefide.pmd
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
trunc_start = int(sys.argv[7])
trunc_length = int(sys.argv[8])
trunc_end = trunc_start + trunc_length

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
	else:
		raise ValueError('File: ' + mov_in + ' invalid. Only .tif and .bin files supported for input.')

	######## PUT MOVIE TRIMMING HERE #######

	row_cut_lower = math.floor((nrows % (2 * rblocks))/2)
	row_cut_upper = raw_mov.shape[0]-math.ceil((nrows % (2 * rblocks))/2)
	col_cut_lower = math.floor((ncols % (2 * cblocks))/2)
	col_cut_upper = raw_mov.shape[1]-math.ceil((ncols % (2 * cblocks))/2)
	raw_mov = raw_mov[row_cut_lower:row_cut_upper,col_cut_lower:col_cut_upper,:]

	intens = np.squeeze(np.mean(np.mean(raw_mov,axis=1),axis=0))
	intensS = np.cumsum(intens,dtype=float)
	intensS[20:] = intensS[20:] - intensS[:-20]
	intensS = intensS[100:] / 20
	intensS = intens[100:] / intensS
	badFrames = np.where(intensS <= 0.991)
	
	raw_mov = raw_mov[:,:,100:]
	for i in badFrames:
		raw_mov[:,:,i] = raw_mov[:,:,i-1]

	print('Movie size: {0}\n'.format(raw_mov.shape))

	if len(sys.argv) <= 9 or sys.argv[9] == '':
		raw_stim = 10 * np.ones(raw_mov.shape[2]) # simulate stimulation values for in vivo data
	else:
		wf = np.fromfile(data_dir + sys.argv[9],dtype="float64")
		raw_stim = (wf.newbyteorder().reshape(2,-1))[0,::10]
		raw_stim = raw_stim[100:]

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
	Y = mov_detr_nnorm
	d1, d2, T = Y.shape
	t = trunc_length;

	# Specify Decomp Parameters
	bheight = 	int(d1 / rblocks)
	hbheight = int(bheight / 2)
	bwidth = int(d2 / cblocks)
	hbwidth = int(bwidth / 2)
	max_comp = 30

	def eval_spatial_stat(u):
		tmp1 = np.abs(u[1:,:] - u[:-1,:])
		tmp2 = np.abs(u[:,1:] - u[:,:-1])
		return (((np.sum(tmp1) + np.sum(tmp2)) * np.prod(u.shape)) / 
				((np.prod(tmp1.shape) + np.prod(tmp2.shape)) * np.sum(np.abs(u))))

	def eval_temporal_stat(v):
		return np.sum(np.abs(v[:-2] + v[2:] - 2*v[1:-1])) / np.sum(np.abs(v))

	num_blocks = rblocks * cblocks
	num_comps = 3
	num_samples = 240
	num_repeats = int(num_samples / (num_blocks * num_comps))

	start = time.time()

	# Iteratively Simulate & Fit Noise To Collect Samples
	spatial_stats = []
	temporal_stats = []

	for rep in range(num_repeats):
		
		# Generate Noise Movie Of NumSimul Blocks
		Y_sim = np.reshape(np.random.randn(num_blocks*bwidth*bheight*t),
							(bheight, num_blocks*bwidth, t))

		# Run Denoiser W/ Max Comp 3 and absurdly large thresholds
		out = trefide.pmd.batch_decompose(bheight, num_blocks*bwidth, t,
											Y_sim, bheight, bwidth,
											1e5, 1e5, 
											num_comps, 3, 10, 40, 5e-3)
	
		# Collect Test Statistics On All Samples
		for bdx in range(num_blocks):
			for cdx in range(num_comps):
				spatial_stats.append(eval_spatial_stat(out[0][bdx,cdx]))
				temporal_stats.append(eval_temporal_stat(out[1][bdx,cdx]))

	spatial_stats = np.array(spatial_stats)
	temporal_stats = np.array(temporal_stats)

	conf = .95

	# Compute Thresholds
	spatial_thresh =  np.percentile(spatial_stats, conf)
	temporal_thresh = np.percentile(temporal_stats, conf)

	print("Simulation took: " + str(time.time()-start) + ' sec')

	start = time.time()

	outs = []

	#Run on standard tiling
	Y_tmp = np.ascontiguousarray(Y[:,:,trunc_start:trunc_end])
	outs.append(trefide.pmd.batch_decompose(d1, d2, t,
											Y_tmp, bheight, bwidth,
											spatial_thresh, temporal_thresh,
											max_comp, 3, 40, 40, 5e-3))

	#Run again on vertical offset
	Y_tmp = np.ascontiguousarray(Y[hbheight:-hbheight,:,trunc_start:trunc_end])
	outs.append(trefide.pmd.batch_decompose(d1-bheight, d2, t,
											Y_tmp, bheight, bwidth,
											spatial_thresh, temporal_thresh, 
											max_comp, 3, 40, 40, 5e-3))

	#Run again on horizontal offset
	Y_tmp = np.ascontiguousarray(Y[:,hbwidth:-hbwidth,trunc_start:trunc_end])
	outs.append(trefide.pmd.batch_decompose(d1, d2-bwidth, t,
											Y_tmp, bheight, bwidth,
											spatial_thresh, temporal_thresh, 
											max_comp, 3, 40, 40, 5e-3))

	# Run again on diagonal offset
	Y_tmp = np.ascontiguousarray(Y[hbheight:-hbheight,hbwidth:-hbwidth,trunc_start:trunc_end])
	outs.append(trefide.pmd.batch_decompose(d1-bheight, d2-bwidth, t,
											Y_tmp, bheight, bwidth,
											spatial_thresh, temporal_thresh, 
											max_comp, 3, 40, 40, 5e-3))

	outs = [list(out) for out in outs]

	# Iterate Over Different "Offsets" Of Overlapping Decomp
	for odx, out in enumerate(outs):
		
		# Determine Index Offsets Induced By Tiling
		vertical_offset = 0 if odx % 2 == 0 else hbheight
		horizontal_offset = 0 if odx < 2 else hbwidth
		
		# IMPORTANT: Appending "Full" V to original outputs, idenically formatted
		out.append(np.zeros(out[1].shape[:-1]+(T,)))  

		# Iterate Over Each Decomposed Block
		for bdx, (rank, block_inds) in enumerate(zip(out[2], out[3])):
			if rank > 0:
				# Extract Rank & Compute Index w.r.t full movie
				rank = int(rank)
				ydx = vertical_offset + int(bheight * block_inds[0])
				xdx = horizontal_offset + int(bwidth * block_inds[1])

				# Extact Block Out From Full FOV Matching Subset Of U
				Y_tmp = np.ascontiguousarray(
					np.reshape(Y[ydx:ydx+bheight,xdx:xdx+bwidth,:], (-1, T))
				)
				U_tmp = np.ascontiguousarray(
					np.reshape(out[0][bdx,:,:,:rank], (-1,rank)).T
				)

				# TODO: Use More Efficient Regression Method
				out[-1][bdx,:rank,:] = np.dot(np.linalg.inv(np.dot(U_tmp, U_tmp.T)),
												np.dot(U_tmp, Y_tmp))
	
	def recompose(denoiser_outputs, weights=None, temporal_index=1):
		#Preallocate
		bheight, bwidth = denoiser_outputs[0].shape[1:3]
		T = denoiser_outputs[temporal_index].shape[-1]
		nbh, nbw = np.max(denoiser_outputs[3], axis=0)
		Y_den = np.zeros((int(nbh+1)*bheight, int(nbw+1)*bwidth, T))
		if weights is None:
			weights = np.ones((bheight,bwidth))
			
		# Accumulate Scaled Blocks
		for bdx, (rank, block_inds) in enumerate(zip(denoiser_outputs[2],
													 denoiser_outputs[3])):
			if rank > 0:
				rank = int(rank)
				ydx = int(bheight * block_inds[0])
				xdx = int(bwidth * block_inds[temporal_index])
				Y_den[ydx:ydx+bheight,xdx:xdx+bwidth,:] = np.dot(
					denoiser_outputs[0][bdx,:,:,:rank] * weights[:,:,None],
					denoiser_outputs[temporal_index][bdx,:rank,:]
				)
		return Y_den

	# Generate Single Quadrant Weights
	ul_weights = np.empty((hbheight, hbwidth), dtype=np.float64)
	for i in range(hbheight):
		for j in range(hbwidth):
			ul_weights[i,j] = min(i, j)+1


	# Construct Full Tile Weights From Quadrant
	tile_weights = np.hstack([np.vstack([ul_weights, 
										 np.flipud(ul_weights)]),
							  np.vstack([np.fliplr(ul_weights), 
										 np.fliplr(np.flipud(ul_weights))])]) 

	# Construct Full FOV Weights By Repeating
	weights = np.tile(tile_weights, (int(d1/bheight), int(d2/bwidth)))

	# Sum All Weights At Get FOV Pixelwise-Normalization
	cumulative_weights = np.zeros((d1,d2))
	cumulative_weights += weights
	cumulative_weights[hbheight:-hbheight,:] += weights[:-bheight, :]
	cumulative_weights[:,hbwidth:-hbwidth] += weights[:, :-bwidth]
	cumulative_weights[hbheight:-hbheight,hbwidth:-hbwidth] += weights[:-bheight, :-bwidth]

	# Compose Original Tiling
	Y_den = recompose(outs[0], weights=tile_weights, temporal_index=-1)
	# Add Horizontal Offset
	Y_tmp = recompose(outs[1], weights=tile_weights, temporal_index=-1)
	Y_den[hbheight:-hbheight,:,:] += Y_tmp
	# Add Vertical Offset
	Y_tmp = recompose(outs[2], weights=tile_weights, temporal_index=-1)
	Y_den[:,hbwidth:-hbwidth,:] += Y_tmp
	# Add Diagonal Offset
	Y_tmp = recompose(outs[3], weights=tile_weights, temporal_index=-1)
	Y_den[hbheight:-hbheight,hbwidth:-hbwidth,:] += Y_tmp

	# Normalize Movie With recombination weights
	Y_den /= cumulative_weights[:,:,None]
	
	print("Denoising took: " + str(time.time()-start) + ' sec')

	start = time.time()

	imio.imsave(out_dir + '/denoised.tif',Y_den)
	#io.savemat(out_dir + '/denoised.mat',{'denoised':mov_denoised})
	imio.imsave(out_dir + '/PMD_residual.tif',Y - Y_den)
	# np.save(out_dir + '/block_ranks.npy', block_ranks)
	
	print("Denoising Saveout took: " + str(time.time()-start) + ' sec')
else:
	Y_den = imio.imread(out_dir + '/denoised.tif')
	d1,d2,T = Y_den.shape
	print('Denoised movie loaded\n')
if T > 15000 and (not os.path.isfile(out_dir + '/denoised_15s.tif')):
	imio.imsave(out_dir + '/denoised_15s.tif',Y_den[:,:,:15000])# * np.squeeze(np.repeat(np.expand_dims(Sn_image,2),15000,axis=2)))
	print('Denoised snippet saved\n')


print('Finished\n')