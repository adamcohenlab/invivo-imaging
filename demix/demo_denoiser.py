import os

import sys
#sys.path.append('../CaImAn/')

import time
import numpy as np

import skimage.io
import skvideo.io

import denoise
import util_plot
import util_movie

def tic():
	global startTime_for_tictoc
	startTime_for_tictoc = time.time()
	return

def toc():
	if 'startTime_for_tictoc' in globals():
	    print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
	else:
	    print ("Toc: start time not set")
	return

def store_output_movie(movie,
						file_name):
	np.savez(file_name,
			data = movie)
	return

def main(mov,
	den_spatial=False,
	den_temporal=False,
	store_outputs=False,
	make_movie=False,
	make_plots=False):

	gHalf=[2,2]
	range_ff=[0.25,0.5]
	nblocks = [10,10]
	greedy = False
	dx = 1
	movie_length = 1000

	# File name outs
	fname_out_wf = data_file_out +'_wf.npz'
	fname_out_pca = data_file_out +'_pca.npz'
	fname_out_movie = data_file_out + '.mp4'

	if den_spatial:
		tic()
		print('Run spatial denoiser')
		mov_wf = denoise.spatial(mov,
								gHalf=gHalf)
		toc()

		if store_outputs:
			store_output_movie(mov_wf,
								fname_out_wf)
		if make_plots:
			util_plot.comparison_plot([mov, mov_wf],
			                          option='corr',
			                          titles_=['Movie', 'WF Movie'],
			                          plot_orientation='vertical',
			                          cbar_orientation='horizontal'
			                         )
	else:
		mov_wf = mov

	if den_temporal:
		print('Calculate/scale wrt noise level')
		tic()
		noise_level = denoise.noise_level(mov_wf,
		                                  range_ff=range_ff)
		mov_nn = mov_wf/noise_level[:,:,np.newaxis]
		toc()
		print('Run temporal denoiser')
		tic()
		mov_d,ranks = denoise.temporal(mov_nn,
		                               nblocks=nblocks,
		                               greedy= greedy,
		                               dx=dx)

		mov_den = mov_d*noise_level[:,:,np.newaxis]
		toc()

		if store_outputs:
			store_output_movie(mov_den,
								fname_out_pca)
	else:
		mov_den = mov_wf

	if make_movie:
		print('Write movie')
		tic()
		util_movie.movie_writer(mov,
					mov_den,
					fname_out_movie,
					movie_length)
		toc()
	return

if __name__ == "__main__":
	## Input movie
	data_file_path ='example_movies/demoMovie.tif'
	data_dir_out = 'output_run/'
	data_file_name = os.path.split(data_file_path)[-1]
	data_base_name = os.path.splitext(data_file_name)[0]
	data_file_out = os.path.join(data_dir_out,
								data_base_name)

	mov = skimage.io.imread(data_file_path).transpose([1,2,0])
	print(mov.shape)

	store_outputs = True
	den_spatial = True
	den_temporal = True
	make_movie = True
	make_plots = False

	main(mov,
		den_spatial=den_spatial,
		den_temporal=den_temporal,
		store_outputs=store_outputs,
		make_movie = make_movie,
		make_plots=make_plots)
	#main(mov,data_file_out)




