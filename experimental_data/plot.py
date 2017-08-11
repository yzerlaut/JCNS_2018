from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
from scipy.signal import convolve2d
import matplotlib.cm as cm


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

def plot_response(args):
    
    fig, ax = plt.subplots(1, figsize=(4.5,3))
    plt.subplots_adjust(bottom=.23, top=.97, right=.85)

    f = loadmat(args.filename)
    data = 1e3*f['matNL'][0]['stim1'][0]
    time = f['matNL'][0]['time'][0].flatten()+args.tshift
    space = f['matNL'][0]['space'][0].flatten()
    smoothing = np.ones((args.Nsmooth, args.Nsmooth))/args.Nsmooth**2
    smooth_data = convolve2d(data, smoothing, mode='same')

    cond = (time>args.t0) & (time<args.t1)
    c = ax.contourf(time[cond], space, smooth_data[:,cond], cmap=cm.viridis)
    plt.colorbar(c, label='VSD signal ($\perthousand$)',
                 ticks=-.3+.3*np.arange(5))

    x1, x2 = ax.get_xlim()
    ax.plot([x1,x1], [0,2], '-', color='gray', lw=4)
    ax.annotate('2mm', (x1,2), rotation=90, fontsize=14)

    if args.with_onset_propag:
        tt, xx = find_latencies_over_space_simple(time, space,
                                                  smooth_data[:,cond],
                                                  signal_criteria=args.signal_criteria,\
                                                  amp_criteria=args.amp_criteria)
        plt.plot(tt+args.tshift, xx, 'o', lw=0, ms=1, color='k')

        # for intervals in [[0,2.3], [2.5,5.7], [5.9,8.5]]:
        #     cond = (xx>intervals[0]) & (xx<intervals[1]) & (tt<20)
        #     pol = np.polyfit(xx[cond], tt[cond]+100, 1)
        #     xxx = np.linspace(xx[cond][0], xx[cond][-1])
        #     plt.plot(np.polyval(pol, xxx), xxx, 'w--', lw=2)

    set_plot(ax, ['bottom'], yticks=[], xlabel='time (ms)')
    if args.SAVE:
        fig.savefig('/Users/yzerlaut/Desktop/'+args.filename.replace('/matx','.matx')+'.png')
    else:
        plt.show()

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--STIM",help="Choose a network model", default='CENTER')
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    parser.add_argument("--no_sim", help="plot only", action="store_true")
    parser.add_argument("--with_onset_propag", action="store_true")
    parser.add_argument("-f", "--filename",help="filename for saving",
                        default='matx2_br131023_NL_dw_1020_23.mat')
    parser.add_argument("--Nsmooth", type=int, default=2)
    parser.add_argument("--tshift", type=float, default=0)
    parser.add_argument("--signal_criteria", type=float, default=0.4)
    parser.add_argument("--amp_criteria", type=float, default=0.6)
    parser.add_argument("--t0", type=float, default=-1000.)
    parser.add_argument("--t1", type=float, default=1000.)

    args = parser.parse_args()
    plot_response(args)
