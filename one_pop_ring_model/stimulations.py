"""
Loads the parameters of the stimulus and experiment !!
"""

import numpy as np

default_params = {\
                  'sX':3., # extension of the stimulus (gaussian in space)
                  'dt':4e-4,
                  'BIN':5e-3, # for markovian formalism
                  'tstop':400e-3,
                  'tstart':150e-3,
                  'RATE_MAX':10.,
                  'Tau1':40e-3,
                  'Tau2':90e-3}

def heaviside(x):
    return 0.5*(1+np.sign(x))

def triple_gaussian(t, X, t0, T1, T2, X0, sX, amplitude):
    return amplitude*(\
                      np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                      np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))*\
                      np.exp(-(X-X0)**2/2./sX**2)


def get_stimulation(X, MODEL, return_print=False):

    params = default_params
    
    if MODEL=='CENTER':
        X0 = X[int(len(X)/2.)]
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff = triple_gaussian(\
                                   t1, X1, params['tstart'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['RATE_MAX'])
    if return_print:
        return params
    else:
        return t, nu_e_aff


all_models = ['CENTER']

import pprint                   
if __name__=='__main__':
    for m in all_models:
        p = get_stimulation(MODEL, return_print=True)                
        print "=============================================="
        print "===----", p['name'], "-----==========="
        print "=============================================="
        pprint.pprint(p)

