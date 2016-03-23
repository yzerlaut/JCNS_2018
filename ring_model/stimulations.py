"""
Loads the parameters of the stimulus and experiment !!
"""

import numpy as np

default_params = {\
                  'sX':1.2, # extension of the stimulus (gaussian in space)
                  'dt':5e-4,
                  'BIN':5e-3, # for markovian formalism
                  'tstop':600e-3,
                  'tstart':50e-3,
                  'amp':10.,
                  'Tau1':70e-3,
                  'Tau2':150e-3}

am_params = {
    'stimuli_shift':8.1, # 2deg by default
    'delay':50e-3
}


def heaviside(x):
    return 0.5*(1+np.sign(x))

def triple_gaussian(t, X, t0, T1, T2, X0, sX, amplitude):
    return amplitude*(\
                      np.exp(-(t-t0-2*T1)**2/2./T1**2)*heaviside(-(t-t0-2*T1))+\
                      np.exp(-(t-t0-2*T1)**2/2./T2**2)*heaviside(t-t0-2*T1))*\
                      np.exp(-(X-X0)**2/2./sX**2)


def get_stimulation(X, MODEL, return_print=False):

    params = default_params

    BASE = MODEL.split('-')[0]
    if len(MODEL.split('-'))==2:
        ARG1, ARG2 = MODEL.split('-')[1], ''
    elif len(MODEL.split('-'))==3:
        ARG1, ARG2 = MODEL.split('-')[1], MODEL.split('-')[2]
    else:
        ARG1, ARG2 = '', ''
    
    if BASE=='CENTER':
        X0 = X[int(len(X)/2.)]
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff = triple_gaussian(\
                                   t1, X1, params['tstart'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['amp'])
        
    elif BASE=='FIRST_STIM':

        if ARG1=='1deg':
            am_params['stimuli_shift'] /= 2
            
        X0 = X[int(len(X)/2.)]-am_params['stimuli_shift']/2.
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff = triple_gaussian(\
                                   t1, X1, params['tstart'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['amp'])
        
    elif BASE=='SECOND_STIM':
        
        if ARG1=='1deg':
            am_params['stimuli_shift'] /= 2
            
        X0 = X[int(len(X)/2.)]+am_params['stimuli_shift']/2.
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff = triple_gaussian(\
                                   t1, X1, params['tstart']+am_params['delay'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['amp'])
        
    elif BASE=='AM':

        if ARG1=='1deg':
            am_params['stimuli_shift'] /= 2
        
        # first stimulus
        X0 = X[int(len(X)/2.)]-am_params['stimuli_shift']/2.
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff1 = triple_gaussian(\
                                   t1, X1, params['tstart'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['amp'])
        # second stimulus
        X0 = X[int(len(X)/2.)]+am_params['stimuli_shift']/2.
        t = np.arange(int((params['tstop'])/params['dt']))*params['dt'] # time array
        X1, t1 = np.meshgrid(X, t)
        nu_e_aff2 = triple_gaussian(\
                                   t1, X1, params['tstart']+am_params['delay'],\
                                   params['Tau1'], params['Tau2'],\
                                   X0, params['sX'], params['amp'])
        nu_e_aff = nu_e_aff1+nu_e_aff2
        
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

