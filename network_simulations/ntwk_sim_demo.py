from __future__ import print_function

from brian2 import *
from time_varying_input import *
import numpy as np

def run_simulation(NRN_exc='LIF', NRN_inh='LIF', NTWK='Vogels-Abbott', DT=0.1, tstop=300,\
                   kick_value=50., kick_duration=30., SEED=1, ext_drive=0., input_rate=None,\
                   n_rec=2):

    seed(SEED)
    
    M = get_connectivity_and_synapses_matrix(NTWK, number=2)
    # number of neurons
    Ne, Ni= int(M[0,0]['Ntot']*(1-M[0,0]['gei'])), int(M[0,0]['Ntot']*M[0,0]['gei'])
    
    exc_neurons, eqs = get_membrane_equation(get_neuron_params(NRN_exc, number=Ne), M[:,0], return_equations=True)
    inh_neurons, eqs = get_membrane_equation(get_neuron_params(NRN_inh, number=Ni), M[:,1], return_equations=True)

    ## INITIAL CONDITIONS
    exc_neurons.Gee, exc_neurons.Gie, exc_neurons.V = '0.*nS', '0.*nS', '-65*mV'
    inh_neurons.Gei, inh_neurons.Gii, inh_neurons.V = '0.*nS', '0.*nS', '-65*mV'
    
    ## FEEDFORWARD EXCITSTORY CONNECTIONS
    time_array = np.arange(int(tstop/DT))*DT
    rate_array = np.array([kick_value*tt/kick_duration+(tt/kick_duration-1)*ext_drive\
                           if tt<kick_duration else 0. for tt in time_array])
    if input_array is None:
        input_array = 0.*rate_array
        
    rate_array += ext_drive+input_array
    
    input_exc, fdfrwd_to_exc, input_inh, fdfrwd_to_inh = \
        build_up_excitatory_feedforward_connections_for_2_pop(\
                            [exc_neurons, inh_neurons], M,
                            time_array, rate_array, SEED=(SEED+1)**2)

    ## RECURRENT CONNECTIONS
    exc_exc, exc_inh, inh_exc, inh_inh = \
      build_up_recurrent_connections_for_2_pop([exc_neurons, inh_neurons], M,\
                                               SEED=(SEED+2)**2) # only for 2 pop !

    # recording
    trace_Vm_exc = StateMonitor(exc_neurons, 'V', record=range(n_rec))
    trace_Vm_inh = StateMonitor(inh_neurons, 'V', record=range(n_rec))
    trace_Ge_exc = StateMonitor(exc_neurons, 'Gee', record=range(n_rec))
    trace_Gi_exc = StateMonitor(exc_neurons, 'Gie', record=range(n_rec))
    trace_Ge_inh = StateMonitor(exc_neurons, 'Gei', record=range(n_rec))
    trace_Gi_inh = StateMonitor(exc_neurons, 'Gii', record=range(n_rec))
    raster_exc = SpikeMonitor(exc_neurons)
    raster_inh = SpikeMonitor(inh_neurons)

    defaultclock.dt = DT*ms
    run(tstop*ms)

    return trace_Vm_exc, trace_Vm_inh, trace_Ge_exc, trace_Gi_exc,\
        trace_Ge_inh, trace_Gi_inh, raster_exc, raster_inh

def plot_ntwk_sim_output(trace_Vm_exc, trace_Vm_inh, trace_Ge_exc, trace_Gi_exc,\
                         trace_Ge_inh, trace_Gi_inh, raster_exc, raster_inh):
    
    # plotting 
    fig1 = figure(figsize=(10,6))
    plot(raster_exc.t/ms, raster_exc.i, '.g', raster_inh.t/ms, raster_inh.i+len(exc_neurons), '.r')
    xlabel('Time (ms)');ylabel('Neuron index')
    fig2 = figure(figsize=(10,5))
    for i in range(n_rec):
        plot(trace_exc.t[1:] / ms, trace_exc[i].V[1:] / mV, 'g')
        plot(trace_inh.t[1:] / ms, trace_inh[i].V[1:] / mV, 'r')

    show()

if __name__=='__main__':

    import sys
    sys.path.append('../')
    from single_cell_models.cell_library import get_neuron_params
    from single_cell_models.cell_construct import get_membrane_equation
    from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
    from synapses_and_connectivity.syn_and_connec_construct import build_up_recurrent_connections_for_2_pop,\
        build_up_recurrent_connections, build_up_poisson_group_to_pop

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

    parser.add_argument("--CONFIG",help="Cell and Network configuration !", default='LIF--LIF--Vogels-Abbott')
    parser.add_argument("--Ne",help="number of excitatory neurons", type=int, default=4000)
    parser.add_argument("--Ni",help="number of inhibitory neurons", type=int, default=1000)
    parser.add_argument("--DT",help="time steps in ms", type=float, default=0.1)
    parser.add_argument("--tstop",help="time of simulation in ms", type=float, default=300)
    parser.add_argument("--kick_value",help=" stimulation (Hz) for the initial kick", type=float, default=0.)
    parser.add_argument("--kick_duration",help=" stimulation duration (ms) for the initial kick", type=float, default=50.)
    parser.add_argument("--ext_drive",help=" stimulation duration (ms) for the initial kick", type=float, default=0.)
    parser.add_argument("--SEED",help="SEED for the simulation", type=int, default=1)

    args = parser.parse_args()
    
    # print args.CONFIG
    run_simulation(NRN_exc=args.CONFIG.split('--')[0], NRN_inh=args.CONFIG.split('--')[1],\
                   NTWK=args.CONFIG.split('--')[2],
                   kick_value=args.kick_value, kick_duration=args.kick_duration,
                   DT=args.DT, tstop=args.tstop, SEED=args.SEED, ext_drive=args.ext_drive)










