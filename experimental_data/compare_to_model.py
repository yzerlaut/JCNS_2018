from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
# from data_analysis.processing.signanalysis import * #gaussian_smoothing
from scipy.signal import convolve2d
import matplotlib.cm as cm
from dataset import get_dataset
from scipy.optimize import minimize

def get_onset_time(t, data, fraction_of_max_criteria=0.3, debug=False):
    spatial_average = np.mean(data, axis=0)
    cond = (spatial_average>fraction_of_max_criteria*spatial_average.max())
    if debug:
        plt.plot(t, spatial_average)
        plt.plot([t[cond][0]], spatial_average[cond][0], 'D')
        plt.show()
    return t[cond][0]

def gaussian(x, amp, mean, std):
    return amp*np.exp(-(x-mean)**2/2./std**2)

from scipy.ndimage.filters import gaussian_filter1d

def gaussian_smoothing(signal, idt_sbsmpl=10):
    """Gaussian smoothing of the data"""
    return gaussian_filter1d(signal, idt_sbsmpl)

def get_stim_center(space, data, Nsmooth=4, debug=False):
    """ we smoothe the average over time and take the x position of max signal"""
    temporal_average = np.mean(data, axis=1)
    smoothed = gaussian_smoothing(temporal_average, Nsmooth)[:-int(Nsmooth)]
    i0 = np.argmax(smoothed)
    x0 = space[:-int(Nsmooth)][i0]
    if debug:
        plt.plot(space, temporal_average)
        plt.plot(space[:-int(Nsmooth)], smoothed)
        plt.plot([x0], [smoothed[i0]], 'D')
        plt.show()
        # def to_minimize(X):
        #     return np.sum((temporal_average-gaussian(space, *X))**2)
        # res = minimize(to_minimize, [5., 1.], options={'maxiter':1e3})
        # plt.plot(space, temporal_average)
        # plt.plot(space, gaussian(space, *res.x))
        # plt.plot(space, gaussian(space, temporal_average.max(), *res.x))
        # plt.show()
        # return res.x[1]
        # x0 = res.x[0]
    return x0

def get_data(dataset_index, Nsmooth=2, t0=-150, t1=100):

    # loading data
    print(get_dataset()[dataset_index])
    f = loadmat(get_dataset()[dataset_index]['filename'])
    data = 1e3*f['matNL'][0]['stim1'][0]
    data[np.isnan(data)] = 0 # blanking infinite data
    time = f['matNL'][0]['time'][0].flatten()
    space = f['matNL'][0]['space'][0].flatten()
    smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
    smooth_data = convolve2d(data, smoothing, mode='same')
    # apply time conditions
    cond = (time>t0) & (time<t1)
    new_time, new_data = np.array(time[cond]), np.array(smooth_data[:,cond])
    # get onset time
    t_onset = get_onset_time(new_time, new_data, debug=args.debug)
    x_center = get_stim_center(space, new_data, debug=args.debug)
    return new_time-t_onset, space-x_center, new_data


def get_residual(args,
                 new_time, space, new_data,
                 Nsmooth=2,
                 fn='../ring_model/data/example_data.npy',
                 with_plot=False):

    # loading model and centering just like in the model
    args2, t, X, Fe_aff, Fe, Fi, muVn = np.load(fn) # we just load a file
    t*=1e3 # bringing to ms
    X -= args2.X_extent/2.+args2.X_extent/args2.X_discretization/2.
    Xcond = (X>=space.min()) & (X<=space.max())
    new_X, new_muVn = X[Xcond], muVn.T[Xcond,:]
    t -= get_onset_time(t, new_muVn)  # centering over time in the same than for data
    print(new_X)
    print(new_muVn)
    # let's construct the spatial subsampling of the data that
    # matches the spatial discretization of the model
    new_data_common_sampling = np.zeros((len(new_X), len(new_time)))
    for i, nx in enumerate(new_X):
        i0 = np.argmin(np.abs(nx-space)**2)
        new_data_common_sampling[i, :] = new_data[i0,:]

    # let's construct the temporal subsampling of the model that
    # matches the temporal discretization of the data
    new_muVn_common_sampling = np.zeros((len(new_X), len(new_time)))
    for i, nt in enumerate(new_time):
        i0 = np.argmin((np.abs(t-nt))**2)
        new_muVn_common_sampling[:, i] = new_muVn[:, i0]

        
    if with_plot:

        fig, AX = plt.subplots(2, figsize=(4.5,5))
        plt.subplots_adjust(bottom=.23, top=.97, right=.85, left=.3)
        plt.axes(AX[0])
        c = AX[0].contourf(new_time, new_X, new_data_common_sampling, cmap=cm.viridis)
        plt.colorbar(c, label='VSD signal ($\perthousand$)', ticks=.5*np.arange(3))
        set_plot(AX[0], [], xticks=[], ylabel='space (mm)')
        plt.axes(AX[1])
        c2 = AX[1].contourf(new_time, new_X, new_muVn_common_sampling, cmap=cm.viridis)
        plt.colorbar(c2, label='VSD signal ($\perthousand$)')
        set_plot(AX[1], ['bottom'], xlabel='time (ms)', ylabel='space (mm)')

        plt.show()

    return np.sum((new_data_common_sampling/new_data_common_sampling.max()-\
                   new_muVn_common_sampling/new_muVn_common_sampling.max())**2)


if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description=
            """
            """,
            formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--Nsmooth", help="for data plots", type=int, default=1)
    parser.add_argument("-a", "--analyze", help="analyze", action="store_true")
    parser.add_argument("-p", "--plot", help="plot analysis", action="store_true")
    parser.add_argument("-d", "--debug", help="with debugging", action="store_true")
    parser.add_argument("--model_filename", '-f', type=str,
                        default='../ring_model/data/example_data.npy')
    parser.add_argument("--data_index", '-df', type=int,
                        default=1)
    parser.add_argument("--t0", type=float, default=-100.)
    parser.add_argument("--t1", type=float, default=150.)
    args = parser.parse_args()

    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)
    print(get_residual(args,
                       new_time, space, new_data,
                       Nsmooth=args.Nsmooth,
                       fn=args.model_filename,
                       with_plot=True))
