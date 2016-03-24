import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt

import sys
sys.path.append('../code/')
from my_graph import set_plot
sys.path.append('../')
from analysis.phase_latency import find_latencies_over_space

def find_latencies_over_space_simple(t, X, signal,\
                              signal_criteria=0.01,\
                              baseline=0, discard=20,\
                              amp_criteria=1./4.):
    signal2 = np.abs(signal)-np.abs(signal[1,:]).mean()
    i_discard = int(discard/(t[1]-t[0]))
    t = t[i_discard:]
    signal2 = signal2[i_discard:,:]-baseline
    XX, TT = [], []
    for i in range(signal2.shape[1]):
        imax = np.argmax(signal2[:,i])
        if signal2[imax,i]>=signal_criteria*signal2.max():
            ii = np.argmin(np.abs(signal2[:imax,i]-amp_criteria*signal2[imax,i]))
            XX.append(X[i])
            TT.append(t[ii]+t[i_discard])
    return TT, XX

def space_time_vsd_style_plot(t, array, zlabel='rate (Hz)',\
                              xlabel='time (ms)', ylabel='cortical space', title='',
                              zlim=None, with_latency_analysis=False,
                              bar_mm=5, cmap = mpl.cm.jet,
                              params={'pixels_per_mm':1.},
                              xzoom=None, yzoom=None):
    """
    takes an array of shape (t, X, Y) and plots its value in 2d for different times !
    at 8 different times
    it returns the figure
    """
    
    # to have the same scale on all plots we normalize all the response with respect to the total mean and total max
    fig = plt.figure(figsize=(7,5))
    plt.suptitle(title, fontsize=22)
    
    ax = plt.subplot2grid((1,8), (0,0), colspan=7)
    
    ax.annotate(str(bar_mm)+'mm', (0.,0.), xycoords='axes fraction') # bar
    if yzoom is None:
        yzoom = [0, array.shape[1]]
    if xzoom is None:
        xzoom = [t[0], t[-1]]
    # now in pixel coordinates
    yzoom = np.array(yzoom)*params['pixels_per_mm']
    bar_mm *= params['pixels_per_mm']
    
    
    if zlim is not None:
        MIN, MAX = zlim
    else:
        MIN = np.min(array[int(xzoom[0]/(t[1]-t[0])):int(xzoom[1]/(t[1]-t[0])), int(yzoom[0]):int(yzoom[1])])
        MAX = np.max(array[int(xzoom[0]/(t[1]-t[0])):int(xzoom[1]/(t[1]-t[0])), int(yzoom[0]):int(yzoom[1])])

    norm = mpl.colors.Normalize(vmin=MIN, vmax=MAX) # to set the color bar

    # in pixel coordinates
    plt.imshow(array.T, vmin=MIN, vmax=MAX, cmap=cmap,
               interpolation='none', aspect='auto', extent=[t[0],t[-1], array.shape[1],0])

    if with_latency_analysis:
        # in pixel coordinates
        TT, XX = find_latencies_over_space(t, np.linspace(0,1,array.shape[1])*array.shape[1],\
                                           array, phase_criteria=-np.pi/2.+np.pi/4.)
        ax.plot(TT, XX, 'w--', lw=3)

    ax.plot([xzoom[0],xzoom[0]], [yzoom[0], yzoom[0]+bar_mm], 'k-', lw=5)
    set_plot(ax, ['bottom'], xlabel=xlabel, ylabel=ylabel, xlim=xzoom, ylim=yzoom, yticks=[])
    
    ax2 = plt.subplot2grid((1, 8), (0, 7), rowspan=1)
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    cb.set_label(zlabel)
    plt.subplots_adjust(left=.12, top=0.87, bottom=0.15, right=.84)
    return ax, fig
        
