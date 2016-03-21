import sys
sys.path.append('../')
from mean_field.euler_method import run_mean_field
from network_simulations.waveform_input import double_gaussian, smooth_heaviside

sys.path.append('../code')
from my_graph import set_plot

import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

def heaviside(x):
    return .5*(1+np.sign(x))
def step_input(t, T0, amp, T1=0.02):
    return amp*heaviside(t-T0)*(1-np.exp(-(t-T0)/T1))
    # return 0*amp*np.exp(-(t-T0)**2/2./T1**2)


t0, T1, T2, tstop = 250e-3, 50e-3, 70e-3, 500e-3

if sys.argv[-1]=='full':

    fig1, [ax1, ax2] = plt.subplots(2, figsize=(5,6))
    plt.subplots_adjust(left=.25, bottom=.25 )

    N=100
    
    max_f_amp, max_vm_amp = np.zeros(N), np.zeros(N)
    amplitudes = np.linspace(0, 15, N)
    
    for i in range(N):
        amp = amplitudes[i]
        def func(t):
            # return step_input(t, 0.02, 1.)*amp
            return double_gaussian(t, t0, T1, T2, amp)
        t, fe, fi, muV, sV, muG, Tv = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3,\
                                                     ext_drive_change=0., PURE_EXC_AFF=False,
                                                     extended_output=True, tstop=tstop)
        max_f_amp[i] = np.max(.8*fe+.2*fi)
        max_vm_amp[i] = np.max(1e2*np.abs((muV-muV[0])/muV[0]))

    ax1.plot(amplitudes, max_f_amp, 'k-', lw=3)
    pol = np.polyfit(amplitudes[:3], max_f_amp[:3], 1)
    ax1.plot(amplitudes[:int(2*N/3.)], np.polyval(pol, amplitudes[:int(2*N/3.)]), 'k--')
    set_plot(ax1, ['left'], ylabel='max. $\\nu$ (Hz)', xticks=[])
    ax2.plot(amplitudes, max_vm_amp, 'k-', lw=3)
    pol = np.polyfit(amplitudes[:3], max_vm_amp[:3], 1)
    ax2.plot(amplitudes[:int(2*N/3.)], np.polyval(pol, amplitudes[:int(2*N/3.)]), 'k--')
    set_plot(ax2, ylabel=r'$\| \delta V/V_0 \| $ %', xlabel='max. $\\nu_e^{aff}$ (Hz)')
    plt.show()

else:

    fig1, [ax1, ax2, ax3] = plt.subplots(3, figsize=(5,6))
    plt.subplots_adjust(left=.25, bottom=.25 )

    for amp in np.linspace(0, 15, 10):
        def func(t):
            # return step_input(t, 0.02, 1.)*amp
            return double_gaussian(t, t0, T1, T2, amp)

        t, fe, fi, muV, sV, muG, Tv = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3,\
                                                     ext_drive_change=0., PURE_EXC_AFF=False,
                                                     extended_output=True, tstop=tstop)
        ax1.plot(1e3*t, func(t), 'k')
        ax2.plot(1e3*t, .8*fe+.2*fi, 'k')
        ax3.plot(1e3*t, 1e2*np.abs((muV-muV[0])/muV[0]), 'k')

    ax1.annotate('external input', (0,0))    
    set_plot(ax1, ['left'], ylabel='$\\nu_e^{aff}$', xticks=[])
    ax2.annotate('network response', (0,0))    
    set_plot(ax2, ['left'], ylabel='$\\nu$', xticks=[])
    ax3.annotate('mean membrane potential', (0,0))    
    set_plot(ax3, ylabel=r'$\| \delta V/V_0 \| $ %', xlabel='time (ms)')
    plt.show()



"""

import sys
sys.path.append('../')

# theoretical tools
import transfer_functions.theoretical_tools as th

# in one pop ring model
from one_pop_ring_model import model as oneD_ring_model 

# list of models
from transfer_functions import neuronal_models # in one pop ring model
from transfer_functions import network_models # in one pop ring model
from one_pop_ring_model import ring_models
from plotting_tools import plot_input_output_response

# just graphics
sys.path.append('/home/yann/work/python_library/')
from my_graph import set_plot

from scipy.optimize import curve_fit

dt, tstop = 5e-4, 150e-3
t = np.arange(int(tstop/dt))*dt
t0_aff = 10e-3

def Euler_method_for_rate_model(t, params,
            exc_aff_level=0., inh_aff_level=0., t0_aff=10e-3,\
            nu_0=10., W_0=0e-12, BIN=5e-3):

    # initial values
    nu_aff_exc = np.array([0 if tt<t0_aff else exc_aff_level for tt in t])
    nu_aff_inh = np.array([0 if tt<t0_aff else inh_aff_level for tt in t])

    muV0, sV0, muGn0, TvN0 = th.get_fluct_regime_vars(nu_0, nu_0, 0, *th.pseq_params(params))

    vm = 0*nu_aff_exc+muV0 # array initialization, membrane potential !
    nu = 0*nu_aff_exc+nu_0
    W = 0*nu+W_0
    
    nu, vm, W = th.make_loop(t, nu, vm, W, nu_aff_exc, nu_aff_inh, BIN, *th.pseq_params(params))
    
    return nu, vm, W, nu_aff_exc, nu_aff_inh


### for network properties prop
def compute_network_response(params, savefig=False,
            levels_exc = np.linspace(0.5, 20., 10),\
            levels_inh = np.linspace(0.5, 10., 10),\
            t0_aff=t0_aff):
    fig1 = plt.figure(figsize=(8,7))
    plt.subplots_adjust(bottom=.3, left=.3, hspace=.3, wspace=.3)
    ax11 = plt.subplot2grid((3, 2), (0,0), rowspan=2)
    ax12 = plt.subplot2grid((3, 2), (2, 0), rowspan=1)
    ax21 = plt.subplot2grid((3, 2), (0,1), rowspan=2)
    ax22 = plt.subplot2grid((3, 2), (2, 1), rowspan=1)

    # first trial for fixed point
    nu, vm, W, nu_aff_exc, nu_aff_inh = Euler_method_for_rate_model(t, params,
                exc_aff_level=0, t0_aff=t0_aff,\
                nu_0=10., BIN=5e-3)
    nu_0 = nu[:-50].mean()
    W_0 = W[:-50].mean()

    def func(t_fit, a, b, c):
        return a*(1-np.exp(-b*t_fit))#+c

    mymap = mpl.colors.LinearSegmentedColormap.from_list(\
                                'mycolors',['black', 'blue'])
    norm = mpl.colors.BoundaryNorm(levels_exc, mymap.N)
    for i in range(len(levels_exc)):
        r = (levels_exc[i]-levels_exc.min())/(levels_exc.max()-levels_exc.min())

        nu, vm, W, nu_aff_exc, nu_aff_inh = Euler_method_for_rate_model(t, params,
                exc_aff_level=levels_exc[i], t0_aff=t0_aff,\
                nu_0=nu_0, W_0=W_0, BIN=5e-3)
        ax12.plot(1e3*t, nu_aff_exc, color=mymap(r,1))
        ax11.plot(1e3*t, nu, color=mymap(r,1))
        imax = np.argmax(nu)
        nu_fit = nu[int(t0_aff/dt):imax]-nu[int(t0_aff/dt)]
        t_fit = t[int(t0_aff/dt):imax]-t0_aff
        popt, pcov = curve_fit(func, t_fit, nu_fit, [10.,1./1e-3,10.])
        ax11.plot(1e3*(t_fit+t0_aff), func(t_fit, *popt)+nu[int(t0_aff/dt)],\
                  ':', color=mymap(r,1))
    set_plot(ax11, ylabel='recurrent activity (Hz)')
    set_plot(ax12, xlabel='time (ms)', ylabel='excitatory \n input (Hz)')

    mymap = mpl.colors.LinearSegmentedColormap.from_list(\
                                'mycolors',['black', 'red'])
    norm = mpl.colors.BoundaryNorm(levels_inh, mymap.N)
    for i in range(len(levels_inh)):
        r = (levels_inh[i]-levels_inh.min())/(levels_inh.max()-levels_inh.min())

        nu, vm, W, nu_aff_exc, nu_aff_inh = Euler_method_for_rate_model(t, params,
                inh_aff_level=levels_inh[i], t0_aff=t0_aff,\
                nu_0=nu_0, W_0=W_0, BIN=5e-3)
        ax22.plot(1e3*t, nu_aff_inh, color=mymap(r,1))
        ax21.plot(1e3*t, nu, color=mymap(r,1))
        imax = np.argmin(nu)
        nu_fit = nu[int(t0_aff/dt):imax]-nu[int(t0_aff/dt)]
        t_fit = t[int(t0_aff/dt):imax]-t0_aff
        popt, pcov = curve_fit(func, t_fit, nu_fit, [10.,1./1e-3,10.])
        ax21.plot(1e3*(t_fit+t0_aff), func(t_fit, *popt)+nu[int(t0_aff/dt)],\
                  ':', color=mymap(r,1))
    set_plot(ax21, xlabel='time (ms)', ylabel='recurrent activity (Hz)')
    set_plot(ax22, ylabel='inhibitory \n input (Hz)')

    if savefig:
        fig1.savefig('../figures/step_reponse_in_network.svg')
    return fig1
    
def compute_input_output_response(params,\
            levels_exc = np.linspace(0.1, 20., 100),\
            levels_inh = np.linspace(0.1, 10., 100),\
            t0_aff=t0_aff, fit=True):
            
    t_fit = t[int(t0_aff/dt):]-t0_aff
    def func(t_fit, a, b, c):
        return a*(1-np.exp(-b*t_fit))+c

    nu, vm, W, nu_aff_exc, nu_aff_inh = Euler_method_for_rate_model(t, params,
                t0_aff=t0_aff, nu_0=10., BIN=5e-3)
    nu_0 = nu[:-50].mean()
    W_0 = W[:-50].mean()
    
    # now more detailed for the function
    Tntwk_exc, Fstat_exc = 0*levels_exc, 0*levels_exc
    
    for i in range(len(levels_exc)):
        nu, vm, W, nu_aff_exc, nu_aff_inh = Euler_method_for_rate_model(t, params,
                exc_aff_level=levels_exc[i], t0_aff=t0_aff,\
                nu_0=nu_0, W_0=W_0, BIN=5e-3)
        imax = np.argmax(nu)
        nu_fit = nu[int(t0_aff/dt):imax]
        t_fit = t[int(t0_aff/dt):imax]-t0_aff
        if fit:
            popt, pcov = curve_fit(func, t_fit, nu_fit, [10.,1./1e-3,10.])
            Tntwk_exc[i] = 1e3/popt[1]
        else:
            Tntwk_exc[i] = 1
        Fstat_exc[i] = nu[-50:].mean()

    Tntwk_inh, Fstat_inh = Tntwk_exc+0, 0*levels_inh
    
    for i in range(len(levels_inh)):
        nu, vm, W, nu_aff_inh, nu_aff_inh = Euler_method_for_rate_model(t, params,
                inh_aff_level=levels_inh[i], t0_aff=t0_aff,\
                nu_0=nu_0, W_0=W_0, BIN=5e-3)
        imax = np.argmin(nu)
        nu_fit = nu[int(t0_aff/dt):imax]
        t_fit = t[int(t0_aff/dt):imax]-t0_aff
        Fstat_inh[i] = nu[-50:].mean()
        if fit:
            popt, pcov = curve_fit(func, t_fit, nu_fit, [10.,1./1e-3,10.])
            if 1e3/popt[1]>=Tntwk_inh[max(i-1,0)]:
                Tntwk_inh[i] = 1e3/popt[1]
            else:
                Tntwk_inh[i] = -np.nan
        else:
            Tntwk_inh[i] = -np.nan
    
    return levels_exc, Fstat_exc, Tntwk_exc, levels_inh,  Tntwk_inh, Fstat_inh


"""


# import argparse

# if __name__=='__main__':
#     # First a nice documentation 
#     parser=argparse.ArgumentParser(description=
#      """ 
#      ----------------------------------------------------------------------
#      We stimulate the rate model with a step input of excitatory firing frequencies
#      In the network response, we investigate the temporal dynamics and
#      the stationary value as a function of the amplitude of the step stimulation
#      ----------------------------------------------------------------------
#      """
#     ,formatter_class=argparse.RawTextHelpFormatter)

#     parser.add_argument("Neuron_Model",help="Choose a neuronal model"+\
#                         "\n      from '../transfer_functions/neuronal_models.py'")
#     parser.add_argument("Network_Model",help="Choose a network model (synaptic and connectivity properties)"+\
#                         "\n      from '../transfer_functions/network_models'.py")

#     args = parser.parse_args()

#     MODEL, NTWK = args.Neuron_Model, args.Network_Model
#     Pnrn = neuronal_models.get_model_params(MODEL)
#     Pntwk = network_models.get_model_params(NTWK)
#     Pring = ring_models.get_model_params('RING1') # but useless
#     params = dict(Pntwk.items()+Pring.items()+Pnrn.items())
#     params['P'] = th.load_TF_coeff(MODEL, NTWK) # coefficents
    
#     fig1 = compute_network_response(params)
#     levels_exc, Fstat_exc, Tntwk_exc, levels_inh,  Tntwk_inh, Fstat_inh =\
#        compute_input_output_response(params)
#     levels = np.concatenate([-levels_inh[::-1], levels_exc])
#     Fstat = np.concatenate([Fstat_inh[::-1], Fstat_exc])
#     Tntwk = np.concatenate([Tntwk_inh[::-1], Tntwk_exc])
#     fig2 = plot_input_output_response(levels, [Fstat], [Tntwk], [MODEL])
#     plt.show()
#     np.save('data/'+params['name_for_tf']+'_'+\
#             params['ntwk_name_for_tf']+'.npy',\
#             np.array([params, levels, Fstat, Tntwk]))

#     # elif sys.argv[-1]=='Adapt':
#     #     VAR = ['', '+b20']
#     #     SET_OF_F, SET_OF_T, SET_OF_LABEL = [], [], []
#     #     for v in VAR:
#     #         MODEL, NTWK = 'IaF'+v, 'CONFIG1'+v
#     #         Pnrn = neuronal_models.get_model_params(MODEL)
#     #         Pntwk = network_models.get_model_params(NTWK)
#     #         Pring = ring_models.get_model_params('RING1') # but useless
#     #         params = dict(Pnrn.items()+Pntwk.items()+Pring.items())
#     #         P_tf = np.load('../transfer_functions/data/fit_'+\
#     #                        params['name_for_tf']+'_'+\
#     #                        params['ntwk_name_for_tf']+'.npy')
#     #         SET_OF_LABEL.append(MODEL)
#     #         levels, Fstat, Tntwk = \
#     #                 compute_input_output_response(params)
#     #         SET_OF_F.append(Fstat) # 
#     #         SET_OF_T.append(Tntwk)
#     #     fig2 = plot_input_output_response(levels, SET_OF_F, SET_OF_T, SET_OF_LABEL) # 
#     #     plt.show()
        

