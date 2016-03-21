import numpy as np
import matplotlib.pylab as plt

import sys
sys.path.append('../')
from single_cell_models.cell_library import get_neuron_params
from single_cell_models.cell_construct import get_membrane_equation
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
from synapses_and_connectivity.syn_and_connec_construct import build_up_recurrent_connections_for_2_pop,\
    build_up_recurrent_connections, build_up_poisson_group_to_pop
sys.path.append('../code')
from signanalysis import gaussian_func
from scipy.special import erf
from ntwk_sim_demo import *
from my_graph import set_plot

def heaviside(x):
    return 0.5*(1+np.sign(x))

def smooth_heaviside(x):
    return 0.5*(1+erf(x))

def smooth_double_gaussian(t, t0, T1, T2, amplitude, smoothing=1e-2):
    return amplitude*(\
                      np.exp(-(t-t0)**2/2./T1**2)*smooth_heaviside(-(t-t0)/smoothing)+\
                      np.exp(-(t-t0)**2/2./T2**2)*smooth_heaviside((t-t0)/smoothing))

def double_gaussian(t, t0, T1, T2, amplitude):
    return amplitude*(\
                      np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                      np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))

def run_simulation_with_input(args, filename='data/1.npy'):
    
    t = np.arange(int(args.tstop/args.DT))*args.DT
    
    input_rate = double_gaussian(t, args.t0, args.T1, args.T2, args.amp)

    M = get_connectivity_and_synapses_matrix(args.CONFIG.split('--')[2])
    ext_drive = M[0,0]['ext_drive']

    run_simulation(\
                   NRN_exc=args.CONFIG.split('--')[0],\
                   NRN_inh=args.CONFIG.split('--')[1],\
                   NTWK=args.CONFIG.split('--')[2],
                   kick_value=0, kick_duration=args.kick_duration,
                   DT=args.DT, tstop=args.tstop, SEED=args.SEED,\
                   ext_drive=ext_drive,input_rate=input_rate, \
                   full_recording=True, n_rec=args.n_rec,\
                   filename=filename)

if __name__=='__main__':


    import argparse
    # First a nice documentation 
    parser=argparse.ArgumentParser(description=
     """ 
     ----------------------------------------------------------------------
     Run the a network simulation using brian2

     Choose CELLULAR and NTWK PARAMETERS from the available libraries
     see  ../synapses_and_connectivity.syn_and_connec_library.py for the CELLS
     see ../synapses_and_connectivity.syn_and_connec_library.py for the NTWK

     Then construct the input as "NRN_exc--NRN_inh--NTWK"
     example: "LIF--LIF--Vogels-Abbott"
     ----------------------------------------------------------------------
     """
    ,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--CONFIG",help="Cell and Network configuration !", default='RS-cell--FS-cell--CONFIG1')
    parser.add_argument("--amp",help="stimulation amplitude in Hz", type=float, default=1.5)
    parser.add_argument("--t0",help="stimulation middle point in ms", type=float, default=600.)
    parser.add_argument("--T1",help="stimulation rise time in ms", type=float, default=40.)
    parser.add_argument("--T2",help="stimulation rise time in ms", type=float, default=100.)
    parser.add_argument("--DT",help="time steps in ms", type=float, default=0.1)
    parser.add_argument("--tstop",help="time of simulation in ms", type=float, default=1200.)
    parser.add_argument("--kick_duration",help=" stimulation duration (ms) for the initial kick", type=float, default=100.)
    parser.add_argument("--SEED",help="SEED for the simulation", type=int, default=5)
    parser.add_argument("-f", "--file",help="filename for saving", default='data/waveform_input_example.npy')
    parser.add_argument("--n_rec",help="number of recorded neurons", type=int, default=3)
    parser.add_argument('-S', "--sim",action='store_true') # FOR SIMULATION

    args = parser.parse_args()

    if args.sim:
        run_simulation_with_input(args, filename=args.file)
    else: # plot
        from plot_single_sim import *
        t0 = args.t0-4*args.T1
        AX, FIG = plot_ntwk_sim_output(*np.load(args.file),\
                                       zoom_conditions=[t0, args.tstop],\
                                       raster_number=400)
        # adding the theoretical eval
        from mean_field.euler_method import run_mean_field, run_mean_field_extended
        def rate_func(t):
            return double_gaussian(t, 1e-3*args.t0, 1e-3*args.T1, 1e-3*args.T2, args.amp)
        
        t, fe, fi = run_mean_field(args.CONFIG.split('--')[0],\
                                   args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
                                   rate_func,
                                   tstop=args.tstop*1e-3)
        AX.plot(1e3*t[1e3*t>t0], fe[1e3*t>t0], 'g-', lw=5, alpha=.4, label='mean field \n pred.')
        AX.plot(1e3*t[1e3*t>t0], fi[1e3*t>t0], 'r-', lw=5, alpha=.4, label='num. sim.')
        AX.plot(1e3*t[1e3*t>t0], .8*fe[1e3*t>t0]+.2*fi[1e3*t>t0], 'k-', lw=5, alpha=.4, label='..')

        if False: # second order mean field
            t, fe, fi, sfe, sfei, sfi = run_mean_field_extended(args.CONFIG.split('--')[0],\
                                       args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
                                       rate_func,
                                       tstop=args.tstop*1e-3)
            AX.fill_between(1e3*t[1e3*t>t0], fe[1e3*t>t0]-sfe[1e3*t>t0], fe[1e3*t>t0]+sfe[1e3*t>t0],\
                            color='g', alpha=.3)
            AX.fill_between(1e3*t[1e3*t>t0], fi[1e3*t>t0]-sfi[1e3*t>t0], fi[1e3*t>t0]+sfi[1e3*t>t0],\
                            color='r', alpha=.3)
        AX.legend(prop={'size':'xx-small'})
        plt.show()
        put_list_of_figs_to_svg_fig(FIG, visualize=False)
