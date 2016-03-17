from master_equation import *
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
import numpy as np
import matplotlib.pylab as plt

def func(t):
    return 20.*np.exp(-(t-1)**2/.1)

def run_mean_field(NRN1, NRN2, NTWK, array_func,\
                   T=5e-3, dt=1e-4, tstop=2):

    # find external drive
    M = get_connectivity_and_synapses_matrix(NTWK)
    ext_drive = M[0,0]['ext_drive']
    
    X0 = find_fixed_point_first_order(NRN1, NRN2, NTWK, exc_aff=ext_drive,\
                                      verbose=False)

    TF1, TF2 = load_transfer_functions(NRN1, NRN2, NTWK)

    t = np.arange(int(tstop/dt))*dt
    
    fe, fi = 0*t+X0[0], 0*t+X0[1]

    def dX_dt_scalar(X, t=0):
        return build_up_differential_operator_first_order(TF1, TF2, T=T)(X,\
                            exc_aff=array_func(t)+ext_drive)
    X = odeint(dX_dt_scalar, X0, t)         # we don't need infodict here

    return t, X[:,0], X[:,1]

if __name__=='__main__':
    
    t, fe, fi = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3)
    plt.plot(t, fe)
    plt.plot(t, fi)

    plt.show()


