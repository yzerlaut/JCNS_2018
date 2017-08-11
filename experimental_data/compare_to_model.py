from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
from scipy.signal import convolve2d
import matplotlib.cm as cm



f = loadmat('matx2_br131023_NL_dw_1020_23.mat')
data = 1e3*f['matNL'][0]['stim1'][0]
time = f['matNL'][0]['time'][0].flatten()+100
space = f['matNL'][0]['space'][0].flatten()
space -= space.max()/2.
Nsmooth = 2
smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
smooth_data = convolve2d(data, smoothing, mode='same')
cond = (time>0) & (time<250)
new_time, new_data = np.array(time[cond]), np.array(smooth_data[:,cond])
# fig, ax = plt.subplots(1, figsize=(4.5,3))
# plt.subplots_adjust(bottom=.23, top=.97, right=.85)
# c = ax.contourf(new_time, space, new_data, cmap=cm.viridis)
# plt.colorbar(c, label='VSD signal ($\perthousand$)', ticks=-.3+.3*np.arange(5))

def get_residual(args, fn='../ring_model/data/example_data.npy', with_plot=False):


    args2, t, X, Fe_aff, Fe, Fi, muVn = np.load(fn) # we just load a file
    print(args2)
    t*=1e3 # bringing to ms
    X -= args2.X_extent/2.+args2.X_extent/args2.X_discretization
    # restricting to the same spatio-temporal space than data
    Xcond = (X>=space.min()) & (X<=space.max())
    # Tcond = (t>=space.min()) and (X<=space.max())
    
    new_X, new_muVn = X[Xcond], muVn.T[Xcond,:]

    # let's construct the spatial subsampling of the data that
    # matches the spatial discretization of the model
    new_data_common_sampling = np.zeros((len(new_X), len(new_time)))
    for i, nx in enumerate(new_X):
        i0 = np.argmin(np.abs(space-nx)**2)
        new_data_common_sampling[i, :] = new_data[i0, :]

    # shifting the time axis to match the onsets ! (by minimizing the residual)
    Residuals = []
    for j in range(len(new_time)-1):
        
        # let's construct the temporal subsampling of the model that
        # matches the temporal discretization of the data
        new_muVn_common_sampling = np.zeros((len(new_X), len(new_time)))
        for i, nt in enumerate(new_time):
            i0 = np.argmin(np.abs(t-nt-new_time[j])**2)
            new_muVn_common_sampling[:, i] = new_muVn[:, i0]
        Residuals.append(np.sum((new_data_common_sampling/new_data_common_sampling.max()-\
                                 new_muVn_common_sampling/new_muVn_common_sampling.max())**2))

    j0 = np.argmin(np.array(Residuals))
    
    if with_plot:
        new_muVn_common_sampling = np.zeros((len(new_X), len(new_time)))
        for i, nt in enumerate(new_time):
            i0 = np.argmin(np.abs(t-nt-new_time[j0])**2)
            new_muVn_common_sampling[:, i] = new_muVn[:, i0]

        fig, AX = plt.subplots(3, figsize=(4.5,8))
        plt.subplots_adjust(bottom=.23, top=.97, right=.85, left=.3)
        plt.axes(AX[0])
        c = AX[0].contourf(new_time, new_X, new_data_common_sampling, cmap=cm.viridis)
        plt.colorbar(c, label='VSD signal ($\perthousand$)', ticks=-.3+.3*np.arange(5))
        set_plot(AX[0], [], xticks=[], ylabel='space (mm)')
        plt.axes(AX[1])
        c2 = AX[1].contourf(new_time, new_X, new_muVn_common_sampling, cmap=cm.viridis)
        plt.colorbar(c2, label='VSD signal ($\perthousand$)')
        set_plot(AX[1], ['bottom'], xlabel='time (ms)', ylabel='space (mm)')
        AX[2].plot(new_time[1:], Residuals, 'o')
        set_plot(AX[2], xlabel='time shift (ms)', ylabel='Residual')

        plt.show()

    return Residuals[j0]



print(get_residual({}, fn='../ring_model/data/example_data.npy', with_plot=True))
