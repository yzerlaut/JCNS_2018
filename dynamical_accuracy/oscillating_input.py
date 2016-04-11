import numpy as np
import sys
sys.path.append('../')
from single_cell_models.cell_library import get_neuron_params
from single_cell_models.cell_construct import get_membrane_equation
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
from synapses_and_connectivity.syn_and_connec_construct import build_up_recurrent_connections_for_2_pop,\
    build_up_recurrent_connections, build_up_poisson_group_to_pop
from network_simulations.ntwk_sim_demo import run_simulation
from network_simulations.plot_single_sim import bin_array
from mean_field.euler_method import run_mean_field

sys.path.append('../code')
from signanalysis import gaussian_func
from my_graph import set_plot
from scipy.optimize import minimize

def heaviside(x):
    return 0.5*(1+np.sign(x))

def sinewave(t, t0, freq, amplitude, phase=0):
    return amplitude*(1-np.cos(2.*np.pi*freq*(t-t0)+phase))*heaviside(t-t0)/2.

def run_simulation_with_input(args, BIN=5):

    t = np.arange(int(args.tstop/args.DT))*args.DT
    
    input_rate = sinewave(t, args.t0, args.freq/1e3, args.amp)

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
    # fe = bin_array(fe, BIN, time_array)
    # fi = bin_array(fi, BIN, time_array)
    # rate_array = bin_array(rate_array, BIN, time_array)
    # time_array = bin_array(time_array, BIN, time_array)
    ## NO NEED TO BIN !!
    
    # now simulation
    
    np.save(args.file, [args, time_array, rate_array, fe, fi])

def find_modulus_and_phase_shift(t_array, f_array, t0, freq, amp0=1., base_amp0=3.,\
                                 full_output=False):

    t, F = t_array[t_array>t0], f_array[t_array>t0]
    def to_minimize(X):
        [base_amp, phase, amp]  = X
        return np.mean((F-base_amp-sinewave(t, t0, freq, amp, phase=phase))**2)
    res = minimize(to_minimize, [amp0, np.pi/2., base_amp0], method='Nelder-Mead')
    
    if full_output:
        new_f = res.x[0] + np.array([0 if tt<t0 else sinewave(tt, t0, freq, res.x[2], phase=res.x[1])\
                          for tt in t_array])
        return t_array, new_f, res.x
    else:
        return res.x

    
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
    parser.add_argument("--freq",help="stimulation amplitude in Hz", type=float, default=2.)
    parser.add_argument("--t0",help="stimulation middle point in ms", type=float, default=500.)
    parser.add_argument("--DT",help="time steps in ms", type=float, default=0.1)
    parser.add_argument("--tstop",help="time of simulation in ms", type=float, default=1000.)
    parser.add_argument("--kick_duration",help=" stimulation duration (ms) for the initial kick", type=float, default=100.)
    parser.add_argument("--SEED",help="SEED for the simulation", type=int, default=5)
    parser.add_argument("-f", "--file",help="filename for saving", default='data/oscillation_example.npy')
    parser.add_argument("--n_rec",help="number of recorded neurons", type=int, default=3)
    parser.add_argument('-S', "--sim",action='store_true') # FOR SIMULATION

    args = parser.parse_args()

    if args.sim:
        run_simulation_with_input(args)
    else: # plot
        import matplotlib.pylab as plt
        fig, ax = plt.subplots(2)
        
        # loading the numerical data
        args, t, rate_input, fe, fi = np.load(args.file)
        # now building theoretical estimate
        def rate_func(t):
            return sinewave(t, 1e-3*args.t0, args.freq, args.amp)
        t_th, fe_th, fi_th = run_mean_field(args.CONFIG.split('--')[0],\
                 args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
                 rate_func, tstop=args.tstop*1e-3)

        
        _, f_fitted, coeffs = find_modulus_and_phase_shift(1e-3*t, .8*fe+.2*fi,\
                                                           1e-3*args.t0, args.freq,\
                                                           full_output=True)
        ax[0].plot(t, f_fitted, 'r-', lw=1)

        # ax[0].plot(t[t>300], .8*fe[t>300]+.2*fi[t>300], 'k-')
        # ax[1].plot(t[t>300], rate_input[t>300], 'k-')
        ax[0].plot(t, .8*fe+.2*fi, 'k-', lw=.5, alpha=.4)
        ax[1].plot(t, rate_input, 'k-')
        
        ax[0].plot(1e3*t_th, .8*fe_th+.2*fi_th, 'k-', lw=3)
        ax[1].plot(1e3*t_th, rate_func(t_th), 'k-', lw=3)
        t_th *=1e3
        # ax[0].plot(t_th[t_th>300],\
        #            .8*fe_th[t_th>300]+.2*fi_th[t_th>300],\
        #            'k-', lw=3, alpha=.5)
        
        plt.show()
