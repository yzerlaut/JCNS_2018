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

from ntwk_sim_demo import *
from my_graph import set_plot

def double_gaussian(t, t0, T1, T2, amplitude):
    return amplitude*np.array([np.exp(-(tt-t0)**2/2./T1**2) if tt<t0 else np.exp(-(tt-t0)**2/2./T2**2) for tt in t])

def run_simulation_with_input(args):
    
    t = np.arange(int(args.tstop/args.DT))*args.DT
    
    input_rate = double_gaussian(t, args.t0, args.T1, args.T2, args.amp)

    M = get_connectivity_and_synapses_matrix(args.CONFIG.split('--')[2])
    ext_drive = M[0,0]['ext_drive']

    trace_Vm_exc, trace_Vm_inh, trace_Ge_exc, trace_Gi_exc,\
        trace_Ge_inh, trace_Gi_inh, raster_exc,\
        raster_inh, time_array, rate_array, rate_exc, rate_inh, M = run_simulation(\
                    NRN_exc=args.CONFIG.split('--')[0],\
                    NRN_inh=args.CONFIG.split('--')[1],\
                    NTWK=args.CONFIG.split('--')[2],
                    kick_value=0, kick_duration=args.kick_duration,
                    DT=args.DT, tstop=args.tstop, SEED=args.SEED,\
                    ext_drive=ext_drive,input_rate=input_rate, \
                    full_recording=True, n_rec=args.n_rec)

    Raster_exc, Raster_inh, Vm_exc, Vm_inh, Ge_exc, Ge_inh, Gi_exc, Gi_inh =\
       transform_to_simple_arrays(trace_Vm_exc, trace_Vm_inh, trace_Ge_exc, trace_Gi_exc,\
                                  trace_Ge_inh, trace_Gi_inh, raster_exc, raster_inh,\
                                  M, n_rec=args.n_rec)
    
    np.save(args.file,
            [time_array, rate_array, rate_exc, rate_inh, Raster_exc, Raster_inh, Vm_exc, Vm_inh, Ge_exc, Ge_inh, Gi_exc, Gi_inh])
    

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
    parser.add_argument("--amp",help="stimulation amplitude in Hz", type=float, default=10.)
    parser.add_argument("--t0",help="stimulation middle point in ms", type=float, default=500.)
    parser.add_argument("--T1",help="stimulation rise time in ms", type=float, default=40.)
    parser.add_argument("--T2",help="stimulation rise time in ms", type=float, default=100.)
    parser.add_argument("--DT",help="time steps in ms", type=float, default=0.1)
    parser.add_argument("--tstop",help="time of simulation in ms", type=float, default=1000.)
    parser.add_argument("--kick_duration",help=" stimulation duration (ms) for the initial kick", type=float, default=100.)
    parser.add_argument("--SEED",help="SEED for the simulation", type=int, default=5)
    parser.add_argument("-f", "--file",help="filename for saving", default='data/example.npy')
    parser.add_argument("--n_rec",help="number of recorded neurons", type=int, default=3)

    args = parser.parse_args()

    run_simulation_with_input(args)
