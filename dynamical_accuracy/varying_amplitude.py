import numpy as np
import matplotlib.pylab as plt

import sys
sys.path.append('../')
from single_cell_models.cell_library import get_neuron_params
from single_cell_models.cell_construct import get_membrane_equation
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
from synapses_and_connectivity.syn_and_connec_construct import build_up_recurrent_connections_for_2_pop,\
    build_up_recurrent_connections, build_up_poisson_group_to_pop
from network_simulations.ntwk_sim_demo import *
from network_simulations.plot_single_sim import bin_array
sys.path.append('../code')
from signanalysis import gaussian_func
from mean_field.euler_method import run_mean_field

from my_graph import set_plot

def heaviside(x):
    return 0.5*(1+np.sign(x))

def double_gaussian(t, t0, T1, T2, amplitude):
    return amplitude*(\
                      np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                      np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))

def run_simulation_with_input(args, BIN=5):

    t = np.arange(int(args.tstop/args.DT))*args.DT
    
    input_rate = double_gaussian(t, args.t0, args.T1, args.T2, args.amp)

    M = get_connectivity_and_synapses_matrix(args.CONFIG.split('--')[2])
    ext_drive = M[0,0]['ext_drive']

    time_array, rate_array, fe, fi = run_simulation(\
                    NRN_exc=args.CONFIG.split('--')[0],\
                    NRN_inh=args.CONFIG.split('--')[1],\
                    NTWK=args.CONFIG.split('--')[2],
                    kick_value=0, kick_duration=args.kick_duration,
                    DT=args.DT, tstop=args.tstop, SEED=args.SEED,\
                    ext_drive=ext_drive,input_rate=input_rate, \
                    full_recording=False)

    # binning it !
    fe = bin_array(fe, BIN, time_array)
    fi = bin_array(fi, BIN, time_array)
    rate_array = bin_array(rate_array, BIN, time_array)
    time_array = bin_array(time_array, BIN, time_array)
    
    # now simulation
    def rate_func(t):
        return double_gaussian(t, 1e-3*args.t0, 1e-3*args.T1, 1e-3*args.T2, args.amp)

    t_th, fe_th, fi_th = run_mean_field(args.CONFIG.split('--')[0],\
                               args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
                               rate_func,
                               tstop=args.tstop*1e-3)
    
    np.save(args.file, [time_array, rate_array, fe, fi, t_th, fe_th, fi_th])

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
    parser.add_argument("--amp",help="stimulation amplitude in Hz", type=float, default=1.)
    parser.add_argument("--t0",help="stimulation middle point in ms", type=float, default=500.)
    parser.add_argument("--T1",help="stimulation rise time in ms", type=float, default=40.)
    parser.add_argument("--T2",help="stimulation rise time in ms", type=float, default=100.)
    parser.add_argument("--DT",help="time steps in ms", type=float, default=0.1)
    parser.add_argument("--tstop",help="time of simulation in ms", type=float, default=1000.)
    parser.add_argument("--kick_duration",help=" stimulation duration (ms) for the initial kick", type=float, default=100.)
    parser.add_argument("--SEED",help="SEED for the simulation", type=int, default=5)
    parser.add_argument("-f", "--file",help="filename for saving", default='data/waveform_input_example.npy')
    parser.add_argument('-S', "--sim",action='store_true') # FOR SIMULATION

    args = parser.parse_args()

    if args.sim:
        args.file = 'data/varying_amp_'+str(round(args.amp))+'.npy'
        run_simulation_with_input(args)
    else: # plot
        fig, ax1 = plt.subplots(figsize=(5,3))
        fig, ax3 = plt.subplots(figsize=(5,3))
        fig, ax2 = plt.subplots(figsize=(5,2))
        import os
        for file in os.listdir("./data/"):
            split = file.split("varying_amp_")
            if (len(split)>1):
                if (float(split[1].split('.npy')[0]) in 2+np.arange(8)*4):
                    t, rate_input, fe, fi, t_th, fe_th, fi_th = np.load('data/'+file)
                    ax1.plot(t[t>300], .8*fe[t>300]+.2*fi[t>300], 'k-')
                    ax3.plot(t[t>300], fi[t>300], 'r-')
                    ax3.plot(t[t>300], fe[t>300], 'g-')
                    ax2.plot(t[t>300], rate_input[t>300], 'k-')
                    t_th *=1e3
                    ax1.plot(t_th[t_th>300], .8*fe_th[t_th>300]+.2*fi_th[t_th>300], 'k-', lw=3, alpha=.5)
        plt.show()
