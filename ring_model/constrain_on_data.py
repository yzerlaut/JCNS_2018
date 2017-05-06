import numpy as np
import matplotlib.pylab as plt
import sys
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
from graphs.plot_export import put_list_of_figs_to_svg_fig

from model import Euler_method_for_ring_model
from ring_model.ring_models import pixels_per_mm
import matplotlib.cm as cm
from scipy.optimize import minimize

from plotting_tools import space_time_vsd_style_plot

baseline_params = {\
                  'sX':1.5, # extension of the stimulus (gaussian in space)
                  'dX1':-1., # extension of the stimulus (gaussian in space)
                  'dX2':1., # extension of the stimulus (gaussian in space)
                  'dt':5e-4,
                  'BIN':5e-3, # for markovian formalism
                  'tstop':400e-3,
                  'tstart':150e-3,
                  'amp':15.,
                  'Tau1':50e-3,
                  'Tau2':150e-3}


""" LOADING EXPERIMENTAL DATA """
from scipy.io import loadmat
from scipy.signal import convolve2d
f = loadmat('../experimental_data/matx2_br131023_NL_dw_1020_23.mat')
data = 1e3*f['matNL'][0]['stim1'][0]
time = f['matNL'][0]['time'][0].flatten()+100
space = f['matNL'][0]['space'][0].flatten()
space  -= space.mean()
cond = (time>-10) & (time<300)
Nsmooth = 4
smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
smooth_data = convolve2d(data, smoothing, mode='same')
time, smooth_data = time[cond], smooth_data[:,cond]
smooth_data_norm = smooth_data/smooth_data.max() # normalize to maximum

def compare_data_and_model(t, X, Fe_aff, Fe, Fi, muVn, with_plot=False):

    muVn_norm = muVn/muVn.max()
    X2 = X-X.mean()

    space_cond = (X2>=space.min()) & (X2<=space.max())
    temp_cond = (1e3*t>=time.min()) & (1e3*t<=time.max())
    SC, TC = np.meshgrid(np.arange(len(X2))[space_cond], np.arange(len(t))[temp_cond])

    time_model, space_model, data_model = 1e3*t[temp_cond], X2[space_cond], muVn_norm[TC, SC].T

    # limiting factor for temporal resolution -> data: find coincident indexes
    icomp_time = []
    for tt in time:
        icomp_time.append(np.argmin((tt-time_model)**2))
    # limiting factor for spatial resolution -> model: find coincident indexes
    icomp_space = []
    for xx in space_model:
        icomp_space.append(np.argmin((xx-space)**2))
        
    icomp_space, icomp_time = np.array(icomp_space), np.array(icomp_time)
    time_model_comp = time_model[icomp_time]
    space_comp = space[icomp_space]
    smooth_data_comp, data_model_comp = smooth_data[icomp_space,:], data_model[:,icomp_time]

    def to_minimize(x):
        return np.sum((x[0]+x[1]*data_model_comp-smooth_data_comp)**2)
    res = minimize(to_minimize, [-0.1, 0.1])
    difference, factors = to_minimize(res.x), res.x
    
    if with_plot:
        fig1, AX = plt.subplots(1, 2)

        c = AX[0].contourf(time, space_comp, smooth_data_comp)
        plt.colorbar(c, label='VSD signal ($\perthousand$)')

        c = AX[1].contourf(time_model_comp, space_model, data_model_comp*factors[1]+factors[0])
        plt.colorbar(c, label='VSD-like signal ($\perthousand$)')

        fig2, AX = plt.subplots(1, 2, figsize=(7,3))

        DT = 50
        for tt in np.arange(-50,200,DT):
            cond = (time>tt) & (time<tt+DT)
            AX[0].plot(space_comp, smooth_data_comp[:,cond].mean(axis=1)+factors[0])
            AX[0].plot(space_model, data_model_comp[:,cond].mean(axis=1)*factors[1])
        
        plt.show()
        return fig1, fig2, difference, factors
    
    else:
        return difference, factors
    
def minimization_procedure(args, N=2):

    t, X, Fe_aff, Fe, Fi, muVn = Euler_method_for_ring_model(\
                                                             args.NRN1, args.NRN2,\
                                                             args.NTWK, args.RING, baseline_params)

    
    baseline_params = {\
                  'sX':1.5, # extension of the stimulus (gaussian in space)
                  'dX1':-1., # extension of the stimulus (gaussian in space)
                  'dX2':1., # extension of the stimulus (gaussian in space)
                  'dt':5e-4,
                  'BIN':5e-3, # for markovian formalism
                  'tstop':400e-3,
                  'tstart':150e-3,
                  'amp':15.,
                  'Tau1':50e-3,
                  'Tau2':150e-3}


    def to_minimize(X):
        X = Tstart, dX1, dX2, sX, Tau1, Tau2
        params = baseline_params.copy()
        params['tstart'] = Tstart
        params['dX1'], params['dX2'] = dX1, dX2
        params['sX'] = sX
        params['Tau1'], params['Tau2'] = Tau1, Tau2
        t, X, Fe_aff, Fe, Fi, muVn = Euler_method_for_ring_model(\
                                                                 args.NRN1, args.NRN2,\
                                                                 args.NTWK, args.RING, params)
        difference, factors = compare_data_and_model(t, X, Fe_aff, Fe, Fi, muVn, with_plot=False)    
        return difference
    
    res = minimize(to_minimize, [3, 6, 2, .9])
    
    return t, X, Fe_aff, Fe, Fi, muVn


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--NRN1",help="Choose a cell model", default='RS-cell')
    parser.add_argument("--NRN2",help="Choose a cell model", default='FS-cell')
    parser.add_argument("--NTWK",help="Choose a network model", default='CONFIG1')
    parser.add_argument("--RING",help="Choose a ring model", default='RING1')
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    parser.add_argument("--no_sim", help="plot only", action="store_true")
    parser.add_argument("-f", "--file",help="filename for saving", default='data/example_data.npy')
    
    args = parser.parse_args()
    
    if not args.no_sim:
        print('simulation [...]')
        t, X, Fe_aff, Fe, Fi, muVn = minimization_procedure(args)
        np.save(args.file, [args, t, X, Fe_aff, Fe, Fi, muVn])
        args2 = args
    else:
        args2, t, X, Fe_aff, Fe, Fi, muVn = np.load(args.file) # we just load a file
        compare_data_and_model(t, X, Fe_aff, Fe, Fi, muVn)
        
    # params = {'pixels_per_mm':pixels_per_mm(args2.RING)}

    # ax, fig1 = space_time_vsd_style_plot(t*1e3, Fe_aff,\
    #                                      title='$\\nu_e^{aff}(x, t)$',\
    #                                      params=params,
    #                                      xlabel='time (ms)', with_latency_analysis=True)
    # ax, fig2 = space_time_vsd_style_plot(t*1e3, .8*Fe+.2*Fi,\
    #                                      title='$\\nu(x, t)$',\
    #                                      params=params,
    #                                      xlabel='time (ms)', with_latency_analysis=True)
    # ax, fig3 = space_time_vsd_style_plot(t*1e3, 1e2*muVn,\
    #                                      xlabel='time (ms)', title='$\delta V / V_0 (x, t)$',\
    #                                      params=params,
    #                                      zlabel='%', with_latency_analysis=True)
    # if args.SAVE:
    #     put_list_of_figs_to_svg_fig([fig1, fig2, fig3], visualize=False)
    #     for i in range(1,4):
    #         # exec("fig"+str(i)+".savefig('fig"+str(i)+".png', dpi=300)")
    #         exec("fig"+str(i)+".savefig('fig"+str(i)+".svg')")
    # else:
    #     plt.show()
