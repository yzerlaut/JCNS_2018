from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
from scipy.signal import convolve2d
import matplotlib.cm as cm


fig, ax = plt.subplots(1, figsize=(4.5,3))
plt.subplots_adjust(bottom=.23, top=.97, right=.85)

f = loadmat('matx2_br131023_NL_dw_1020_23.mat')
data = 1e3*f['matNL'][0]['stim1'][0]
time = f['matNL'][0]['time'][0].flatten()+100
space = f['matNL'][0]['space'][0].flatten()
Nsmooth = 2
smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
smooth_data = convolve2d(data, smoothing, mode='same')

cond = (time>-10) & (time<250)
c = ax.contourf(time[cond], space, smooth_data[:,cond], cmap=cm.viridis)
plt.colorbar(c, label='VSD signal ($\perthousand$)', ticks=-.3+.3*np.arange(5))

ax.plot([0,0], [0,2], '-', color='gray', lw=4)
ax.annotate('2mm', (0,2), rotation=90, fontsize=14)

def find_latencies_over_space_simple(t, X, signal,\
                                     signal_criteria=0.4,\
                                     amp_criteria=.6):
    signal2 = signal.T
    # i_discard = int(discard/(t[1]-t[0]))
    # t = t[i_discard:]
    # signal2 = signal2[i_discard:,:]-baseline
    XX, TT = [], []
    for i in range(signal2.shape[1]):
        imax = np.argmax(signal2[:,i])
        if signal2[imax,i]>=signal_criteria:
            ii = np.argmin(np.abs(signal2[:imax,i]-amp_criteria*signal2[imax,i]))
            XX.append(X[i])
            TT.append(t[ii])
    return np.array(TT), np.array(XX)

tt, xx = find_latencies_over_space_simple(time, space, smooth_data[:,cond])

plt.plot(tt+100, xx, 'o', lw=0, ms=1, color='k')

for intervals in [[0,2.3], [2.5,5.7], [5.9,8.5]]:
    cond = (xx>intervals[0]) & (xx<intervals[1]) & (tt<20)
    pol = np.polyfit(xx[cond], tt[cond]+100, 1)
    xxx = np.linspace(xx[cond][0], xx[cond][-1])
    plt.plot(np.polyval(pol, xxx), xxx, 'w--', lw=2)


set_plot(ax, ['bottom'], yticks=[], xlabel='time (ms)')
fig.savefig('/Users/yzerlaut/Desktop/3.svg')
plt.show()



