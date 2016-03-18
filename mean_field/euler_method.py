from master_equation import *
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
from single_cell_models.cell_library import get_neuron_params
from transfer_functions.theoretical_tools import get_fluct_regime_vars, pseq_params
from transfer_functions.tf_simulation import reformat_syn_parameters
import numpy as np

def func(t):
    return 2.*np.exp(-(t-1)**2/.1)

def run_mean_field(NRN1, NRN2, NTWK, array_func,\
                   T=5e-3, dt=1e-4, tstop=2, extended_output=False,\
                   ext_drive_change=0.):

    # find external drive
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    ext_drive = M[0,0]['ext_drive']
    
    X0 = find_fixed_point_first_order(NRN1, NRN2, NTWK, exc_aff=ext_drive+ext_drive_change,\
                                      verbose=False)

    TF1, TF2 = load_transfer_functions(NRN1, NRN2, NTWK)

    t = np.arange(int(tstop/dt))*dt
    
    def dX_dt_scalar(X, t=0):
        return build_up_differential_operator_first_order(TF1, TF2, T=T)(X,\
                               exc_aff=ext_drive+ext_drive_change, pure_exc_aff=array_func(t))
    fe, fi = odeint(dX_dt_scalar, X0, t).T         # we don't need infodict here

    if extended_output:
        params = get_neuron_params(NRN2, SI_units=True)
        reformat_syn_parameters(params, M)
        muV, sV, muGn, TvN = get_fluct_regime_vars(fe+ext_drive+ext_drive_change,\
                                                   fi,\
                                                   *pseq_params(params))
        return t, fe, fi, muV, sV, muGn, TvN
    else:
        return t, fe, fi

if __name__=='__main__':
    
    import matplotlib.pylab as plt
    t, fe, fi, muV, sV, _, _ = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3,\
                               extended_output=True)
    plt.figure()
    plt.plot(t, fe, 'g')
    plt.plot(t, fi, 'r')
    plt.figure()
    plt.plot(t, 1e3*muV)
    plt.show()

