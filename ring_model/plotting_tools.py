import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt

import sys
sys.path.append('../analysis/')
from phase_latency import find_latencies_over_space

def space_time_vsd_style_plot(t, array, zlabel='rate (Hz)',\
                              xlabel='time (ms)', ylabel='cortical distance (mm)', title='',
                              zlim=None, with_latency_analysis=False,\
                              phase_criteria=-np.pi/2.+np.pi/6.,
                              params=None, xzoom=None, yzoom=None):
    """
    takes an array of shape (t, X, Y) and plots its value in 2d for different times !
    at 8 different times
    it returns the figure
    """
    
    # to have the same scale on all plots we normalize all the response with respect to the total mean and total max
    fig = plt.figure(figsize=(7,5))
    plt.suptitle(title, fontsize=22)
    cmap = mpl.cm.jet # color map
    if zlim is not None:
        MIN, MAX = zlim[0], zlim[1]
    else:
        MIN, MAX = array.min(), array.max()

    norm = mpl.colors.Normalize(vmin=MIN, vmax=MAX) # to set the color bar

    if params is not None:
        ylim0 = params['X_extent']
    else:
        ylim0 = array.shape[1]

    ax = plt.subplot2grid((1,7), (0,0), colspan=6)
    plt.imshow(array.T, vmin=MIN, vmax=MAX, cmap=cmap,
        interpolation='none', aspect='auto', extent=[t[0],t[-1], ylim0,0])
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    if with_latency_analysis:
        TT, XX = find_latencies_over_space(t, np.linspace(0,1,array.shape[1])*ylim0,\
                                           array, phase_criteria=phase_criteria)
        plt.plot(TT, XX, 'w--', lw=3)

    if xzoom is not None:
        plt.xlim(xzoom)
    if yzoom is not None:
        plt.ylim(yzoom)
        
    ax2 = plt.subplot2grid((1, 7), (0, 6), rowspan=1)
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    cb.set_label(zlabel)
    plt.subplots_adjust(left=.12, top=0.87, bottom=0.15, right=.84)
    return ax, fig
        
