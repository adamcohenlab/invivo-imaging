import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import animation

from IPython.display import HTML
import pylab as pl

from mpl_toolkits.axes_grid1 import make_axes_locatable

import skvideo.io
import util_plot
import cv2

def movie_writer(mov,
                mov_den,
                fname_out_movie,
                video_length,
                min_unit=255,
                 ax_concat=1,
                min_single_frame=False):

    # Truncate movie
    video_length = min(mov.shape[2],
                        video_length)

    mov = mov[:,:,:video_length].transpose([2,0,1])
    mov_den = mov_den[:,:,:video_length].transpose([2,0,1])
    mov_res = mov-mov_den

    ax_ = (0)
    if min_single_frame:
        mov = mov - mov.min(axis=ax_,keepdims=True)
        mov_den = mov_den - mov_den.min(axis=ax_,keepdims=True)
        mov_res = mov_res - mov_res.min(axis=ax_,keepdims=True)

    mov = (mov - mov.min())/(mov.max() - mov.min())*min_unit;
    mov_den = (mov_den - mov_den.min())/(mov_den.max() - mov_den.min())*min_unit
    mov_res = (mov_res - mov_res.min())/(mov_res.max() - mov_res.min())*min_unit

    #mov = video_normalization(mov)
    #mov_den = video_normalization(mov_den)
    #mov_res = video_normalization(mov_res)

    movie_frames = np.concatenate([mov,
                                mov_den,
                                mov-mov_den
                                ],
                                axis=ax_concat)

    writer = skvideo.io.FFmpegWriter(fname_out_movie)

    for frame in range(video_length):
        writer.writeFrame(movie_frames[frame,:,:])
    writer.close()
    return


def play(movie, gain=3, fr=120, offset=0, magnification=3,
    frame_range=[350,1000]):
    maxmov = np.max(movie)
    looping=True
    terminated=False
    while looping:
        for t in range(frame_range[0], frame_range[1]):
            if magnification != 1:
                frame = cv2.resize(movie[:,:,t],
                                   None,
                                   fx=magnification,
                                   fy=magnification,
                                   interpolation=cv2.INTER_LINEAR)
            cv2.imshow('frame', (frame - offset) / maxmov*gain)
            if cv2.waitKey(int(1. / fr * 1000)) & 0xFF == ord('q'):
                looping = False
                terminated = True
                break
        if terminated:
            break

    cv2.waitKey(100)
    cv2.destroyAllWindows()
    for i in range(10):
        cv2.waitKey(100)
    return


def video_normalization(Y,
                        smin=True,
                        dmax=True):
    """
    Normalize video
    min is wrt all pixels in all frames
    max is wrt all pixels in all frames
    """
    if len(Y) ==0:
        return Y

    Yout = np.copy(Y)
    if smin is True:
        Yout = Yout-Yout.min()
    if dmax is True:
        Yout = np.true_divide(Yout,Yout.max())
    return Yout


def draw_video_init(video,
                    max_time=None,
                    xtick=None,
                    ytick=None,
                    share_scale=False,
                    gains=1):
    wsize , hsize , num_tframe = video[0][0].shape
    # define num_plots (num rows and num cols)
    num_rows = len(video)
    if num_rows > 1:
        num_cols = max(map(len, video))
    else:
        num_cols = len(video[0])

    num_plots = num_rows*num_cols

    wratio =[wsize/wsize]*num_cols
    hratio =[hsize/wsize]*num_rows

    fig = plt.figure(figsize=(6*num_cols,
                            6*num_rows))
    gs1 = gridspec.GridSpec(num_rows, num_cols,
                           width_ratios = wratio,
                            height_ratios = hratio)

    if isinstance(gains,int):
        gains = np.ones(shape=(num_rows,num_cols))
    else:
        print('gains were given')

   # Initialize subplot list
    ims = [None]*num_plots

    # for each
    counter = 0

    for row in range(num_rows):
        for col in range(num_cols):
            cvideo = video[row][col]
            if np.any(cvideo):
                ax = plt.subplot(gs1[counter])
                if xtick is None and ytick is None:
                    ax.axis('off')
                else:
                    if ~(xtick is None):
                        ax.set_xticks(xtick[1:])
                        ax.tick_params(bottom='on',top='on')
                    if ~(ytick is None):
                        ax.set_yticks(ytick[1:])
                        ax.tick_params(right='on',left='on')
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])

                if share_scale is False:
                    vmin = cvideo.min()
                    vmax = cvideo.max()
                    cmin = cvideo[:,:,0].min()
                cvideo = cvideo*gains[row][col]/cvideo.max()

                # Initialize the video
                ims[counter] = ax.imshow(cvideo[:,:,0],cmap=plt.cm.gist_gray,
                                        interpolation='nearest',
                                        vmin=vmin,vmax=vmax)
            counter += 1

    gs1.update(wspace =0.1, hspace =0)
    # remove if empty

    def remove_empty(x):
        return np.any(x)

    ims = list(filter(remove_empty,ims))

    # Write animator function
    def animate(frame):
        counter = 0
        for row in range(num_rows):
            for col in range(num_cols):
                cvideo = video[row][col]
                if np.any(cvideo):
                    data_in = cvideo[:, :, frame]
                    ims[counter].set_data(data_in*gains[row][col]/data_in.max())
                    counter += 1
        return ims

    if max_time is None:
        num_frames = num_tframe
    else:
        num_frames = min(max_time,num_tframe)
    print('Making {} frames'.format(num_frames))

    return animation.FuncAnimation(fig, animate, frames=range(num_frames), interval=30, blit=True)


def draw_video(anim, fps=10, store=False, name='tmp.mp4'):
    if store:
        anim.save(name, fps=fps, extra_args=['-vcodec', 'libx264'])
        #video = open(name, "rb").read()
        #anim._encoded_video = video.encode("base64")
    return HTML(anim.to_html5_video())


def Vt_denoised_video(Vt,Vts):
    """
    Make video of Vt and Vts
    Vt = num_componentx x num_timeframes
    """
    num_components, num_timeframes= Vt.shape

    Vdiff = Vt-Vts

    tts = np.arange(num_timeframes)

    fig, ax = plt.subplots(2,1,figsize=(15,10))

    ymax, ymin = Vt.max(), Vt.min()
    ydmax, ydmin = Vdiff.max(), Vdiff.min()

    ax[0].set_ylim(ymin,ymax)
    ax[1].set_ylim(ydmin,ydmax)
    l1,l2,l3=[],[],[]

    l1, = ax[0].plot(tts,Vt[0,:],'b')
    l2, = ax[0].plot(tts,Vts[0,:],'r')
    l3, = ax[1].plot(tts,Vdiff[0,:],'r')

    ax[1].legend('line')

    def animate(ii):
        l1.set_ydata(Vt[ii,:])
        l2.set_ydata(Vts[ii,:])
        l3.set_ydata(Vdiff[ii,:])
        l3.set_label('Iter %d '%(ii))
        #l3.set_label('Tile %d sv: %d'%(num[ii], num2[ii]))
        ax[1].legend()
        return l1,l2,l3


    # Init only required for blitting to give a clean slate.
    def init():
        l1.set_ydata(Vt[0,:])
        l2.set_ydata(Vts[0,:])
        l3.set_ydata(Vdiff[0,:])
        l1.set_label('line init')
        ax[1].legend()
        return l1,l2,l3


    anim = animation.FuncAnimation(fig, animate,
            frames=np.arange(num_components-1), init_func=init,blit=False,
            interval=1000)
    # HTML(anim.to_html5_video())
    # anim.save('traces_Vt_Vtd.mp4',fps=5,extra_args=['-vcodec','libx264'])
    return anim