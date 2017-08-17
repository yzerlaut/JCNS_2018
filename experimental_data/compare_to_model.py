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

from scipy.ndimage.filters import gaussian_filter1d
def gaussian_smoothing(signal, idt_sbsmpl=10):
    """Gaussian smoothing of the data"""
    return gaussian_filter1d(signal, idt_sbsmpl)

def get_time_max(t, data, debug=False, Nsmooth=1):
    spatial_average = np.mean(data, axis=0)
    smoothed = gaussian_smoothing(spatial_average, Nsmooth)[:-int(Nsmooth)]
    i0 = np.argmax(smoothed)
    t0 = t[:-int(Nsmooth)][i0]
    if debug:
        plt.plot(t, spatial_average)
        plt.plot(t[:-int(Nsmooth)], smoothed)
        plt.plot([t0], [smoothed[i0]], 'D')
        plt.show()
    return t0

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
    return x0

def get_data(dataset_index,
             t0=-150, t1=100, debug=False,\
             Nsmooth=2,
             smoothing=None):

    # loading data
    print(get_dataset()[dataset_index])
    delay = get_dataset()[dataset_index]['delay']
    f = loadmat(get_dataset()[dataset_index]['filename'])
    data = 1e3*f['matNL'][0]['stim1'][0]
    data[np.isnan(data)] = 0 # blanking infinite data
    time = f['matNL'][0]['time'][0].flatten()
    space = f['matNL'][0]['space'][0].flatten()
    if smoothing is None:
        smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
    smooth_data = convolve2d(data, smoothing, mode='same')
    # apply time conditions
    cond = (time>t0-delay) & (time<t1-delay)
    new_time, new_data = np.array(time[cond]), np.array(smooth_data[:,cond])
    # get onset time
    t_onset = get_time_max(new_time, new_data, debug=debug)
    x_center = get_stim_center(space, new_data, debug=debug)
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
    t -= get_time_max(t, new_muVn)  # centering over time in the same than for data
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
        c = AX[0].contourf(new_time, new_X, new_data_common_sampling,
                           cmap=cm.viridis,
                           levels=\
                           np.linspace(new_data_common_sampling.min(),
                                       new_data_common_sampling.max(), 30))
        plt.colorbar(c, label='VSD signal ($\perthousand$)',
                     ticks=.5*np.arange(3))
        # set_plot(AX[0], xticks_labels=[], ylabel='space (mm)')
        set_plot(AX[0], ylabel='space (mm)')
        plt.axes(AX[1])
        c2 = AX[1].contourf(new_time, new_X, new_muVn_common_sampling,
                            cmap=cm.viridis,
                            levels=\
                            np.linspace(0, new_muVn_common_sampling.max(), 30))
        # plt.colorbar(c2, label='$\\delta V_N$',
        #              ticks=np.linspace(0, new_muVn_common_sampling.max(), 5))
        set_plot(AX[1], xlabel='time (ms)', ylabel='space (mm)')

        if args.save:
            fig.savefig('/Users/yzerlaut/Desktop/temp.svg')
        else:
            plt.show()

    return np.sum((new_data_common_sampling-new_muVn_common_sampling)**2)

def get_space_residual(args,
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
    t -= get_time_max(t, new_muVn)  # centering over time in the same than for data
    # let's construct the spatial subsampling of the data that
    # matches the spatial discretization of the model
    new_data_common_sampling = np.zeros((len(new_X), len(new_time)))
    for i, nx in enumerate(new_X):
        i0 = np.argmin(np.abs(nx-space)**2)
        # new_data_common_sampling[i, :] = new_data[i0,:]
        # normalizing by local maximum over time
        new_data_common_sampling[i, :] = new_data[i0,:]/new_data[i0,:].max()

    # let's construct the temporal subsampling of the model that
    # matches the temporal discretization of the data
    new_muVn_common_sampling = np.zeros((len(new_X), len(new_time)))
    for i, nt in enumerate(new_time):
        i0 = np.argmin((np.abs(t-nt))**2)
        new_muVn_common_sampling[:, i] = new_muVn[:, i0]
    # normalizing by local maximum over time
    for i, nx in enumerate(new_X):
        new_muVn_common_sampling[i, :] /= new_muVn_common_sampling[i,:].max()
    
    if with_plot:

        fig, AX = plt.subplots(2, figsize=(4.5,5))
        plt.subplots_adjust(bottom=.23, top=.97, right=.85, left=.3)
        plt.axes(AX[0])
        c = AX[0].contourf(new_time, new_X, new_data_common_sampling, cmap=cm.viridis)
        plt.colorbar(c, label='norm. VSD signal',
                     ticks=.5*np.arange(3))
        # set_plot(AX[0], xticks_labels=[], ylabel='space (mm)')
        set_plot(AX[0], ylabel='space (mm)')
        plt.axes(AX[1])
        c2 = AX[1].contourf(new_time, new_X, new_muVn_common_sampling, cmap=cm.viridis)
        plt.colorbar(c2, label='norm. $\\delta V_N$',
                     ticks=.5*np.arange(3))
        set_plot(AX[1], xlabel='time (ms)', ylabel='space (mm)')

        if args.save:
            fig.savefig('/Users/yzerlaut/Desktop/temp.svg')
        else:
            plt.show()


    return np.sum((new_data_common_sampling/new_data_common_sampling.max()-\
                   new_muVn_common_sampling/new_muVn_common_sampling.max())**2)

def heaviside(x):
    return 0.5*(1+np.sign(x))
def double_gaussian(t, t0, T1, T2, amplitude):
    return amplitude*(np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                      np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))

def get_time_residual(args,
                      new_time, space, new_data,
                      Nsmooth=2,
                      fraction_of_max_criteria = 0.8,
                      fn='../ring_model/data/example_data.npy',
                      with_plot=False, return_fit_directly=False):

    Xcond_data = np.argwhere(\
            new_data.max(axis=1)>fraction_of_max_criteria*new_data.max()).flatten()
    if args.debug:
        print(Xcond_data)
        
    new_data_common_sampling = new_data[Xcond_data,:].mean(axis=0)
    # normalizing by local maximum
    new_data_common_sampling /= new_data_common_sampling.max()

    if return_fit_directly:
        
        def to_minimize(X):
            """ X are the parameters """
            return np.sum((double_gaussian(new_time, *X)-new_data_common_sampling)**2)
        res = minimize(to_minimize,
                 x0= [0., 10., 100., 1.])
        return res.x[1], res.x[2]

    else:
        # loading model and centering just like in the model
        args2, t, X, Fe_aff, Fe, Fi, muVn = np.load(fn) # we just load a file
        t*=1e3 # bringing to ms
        X -= args2.X_extent/2.+args2.X_extent/args2.X_discretization/2.
        Xcond = (X>=space.min()) & (X<=space.max())
        new_X, new_muVn = X[Xcond], muVn.T[Xcond,:]
        t -= get_time_max(t, new_muVn)  # centering over time in the same than for data
        imodel = np.argmax(new_muVn.max(axis=1))
        if args.debug:
            print(imodel)
        new_muVn_common_sampling0 = new_muVn[imodel, :]
        # # let's construct the temporal subsampling of the model that
        # # matches the temporal discretization of the data
        new_muVn_common_sampling = np.zeros(len(new_time))
        for i, nt in enumerate(new_time):
            i0 = np.argmin((np.abs(t-nt))**2)
            new_muVn_common_sampling[i] = new_muVn_common_sampling0[i0]
        new_muVn_common_sampling /= new_muVn_common_sampling.max()


        if with_plot:

            fig, ax = plt.subplots(figsize=(4,2.7))
            plt.subplots_adjust(left=.3, bottom=.3)
            tt = np.linspace(new_time.min(), new_time.max(), 1e3)
            if return_fit_directly:
                ax.plot(tt, double_gaussian(tt, *res.x), lw=4, alpha=.8, label='Model')
                ax.annotate('$\\tau_1$='+str(round(res.x[1]))+'ms\n$\\tau_2$='+str(round(res.x[2]))+'ms', (0., 0.8))
            ax.plot(new_time, new_data_common_sampling, 'k-', lw=2, label='data')
            ax.plot(new_time, new_muVn_common_sampling, 'k:', lw=2, label='Model')
            # ax.plot(new_time, double_gaussian(new_time, *res.x), lw=2, label='fit')
            ax.legend()
            set_plot(ax, xlabel='time (ms)', ylabel='norm. signal', yticks=[0, 0.5, 1.])
            if args.save:
                fig.savefig('/Users/yzerlaut/Desktop/temp.svg')
            else:
                plt.show()

        else:
            return np.sum((new_data_common_sampling-new_muVn_common_sampling)**2)


if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description=
            """
            """,
            formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--Nsmooth", help="for data plots", type=int, default=1)
    parser.add_argument("-s", "--save", help="save fig", action="store_true")
    parser.add_argument("-a", "--analyze", help="analyze", action="store_true")
    parser.add_argument("-p", "--plot", help="plot analysis", action="store_true")
    parser.add_argument("-d", "--debug", help="with debugging", action="store_true")
    parser.add_argument("--space", help="space residual", action="store_true")
    parser.add_argument("--time", help="temporal residual", action="store_true")
    parser.add_argument("--model_filename", '-f', type=str,
                        default='../ring_model/data/example_data.npy')
    parser.add_argument("--data_index", '-df', type=int,
                        default=1)
    parser.add_argument("--t0", type=float, default=-100.)
    parser.add_argument("--t1", type=float, default=200.)
    args = parser.parse_args()

    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1,
                                         debug=args.debug)
    if args.space:
        print(get_space_residual(args,
                             new_time, space, new_data,
                             Nsmooth=args.Nsmooth,
                             fn=args.model_filename,
                             with_plot=True))
    elif args.time:
        new_time, space, new_data = get_data(args.data_index,
                                             smoothing=np.ones((1, 4))/4**2,
                                             t0=args.t0, t1=args.t1,
                                             debug=args.debug)
        print(get_time_residual(args,
                                new_time, space, new_data,
                                Nsmooth=2,
                                fn=args.model_filename,
                                with_plot=True))
    else:
        print(get_residual(args,
                             new_time, space, new_data,
                             Nsmooth=args.Nsmooth,
                             fn=args.model_filename,
                             with_plot=True))
