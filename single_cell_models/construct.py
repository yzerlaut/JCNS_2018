"""
This file construct the equations for brian2
"""
from __future__ import print_function
import numpy as np
import brian2

def get_membrane_equation(neuron_params, synaptic_array,\
                          return_equations=False):

    ## pure membrane equation
    if neuron_params['deltaV']==0:
        # if hard threshold : Integrate and Fire
        eqs = """
        dV/dt = (%(gL)f*nS*(%(El)f*mV - V) + I - w_adapt)/(%(Cm)f*pF) : volt (unless refractory) """ % neuron_params
    else:
        eqs = """
        dV/dt = (%(gL)f*nS*(%(El)f*mV - V) + %(gL)f*nS*%(deltaV)f*mV*exp(-(%(Vthre)f*mV-V)/(%(deltaV)f*mV)) + I - w_adapt)/(%(Cm)f*pF) : volt (unless refractory) """ % neuron_params

    ## Adaptation current
    if neuron_params['tauw']>0: # adaptation current or not ?
        eqs += """
        dw_adapt/dt = ( -%(a)f*nS*( %(El)f*mV - V) - w_adapt )/(%(tauw)f*ms) : amp  """ % neuron_params
    else:
        eqs += """
        w_adapt : amp  """

    ## synaptic currents, 1) adding all synaptic currents to the membrane equation via the I variable
    eqs += """
        I = I0 """
    for synapse in synaptic_array:
        # loop over each presynaptic element onto this target
        Gsyn = 'G'+synapse['name']
        eqs += '+'+Gsyn+'*(%(Erev)f*mV - V)' % synapse
    eqs += ' : amp'
    
    ## synaptic currents, 2) constructing the temporal dynamics of the synaptic conductances
    ## N.B. VALID ONLY FOR EXPONENTIAL SYNAPSES UNTIL NOW !!!!
    for synapse in synaptic_array:
        # loop over each presynaptic element onto this target
        Gsyn = 'G'+synapse['name']
        eqs += """
        """+'d'+Gsyn+'/dt = -'+Gsyn+'*(1./(%(Tsyn)f*ms)) : siemens' % synapse
    eqs += """
        I0 : amp """
    # adexp, pratical detection threshold Vthre+5*deltaV
    neurons = brian2.NeuronGroup(neuron_params['N'], model=eqs,\
                                 refractory=str(neuron_params['Trefrac'])+'*ms',
                                 threshold='V>'+str(neuron_params['Vthre']+5.*neuron_params['deltaV'])+'*mV',
                                 reset='V='+str(neuron_params['Vreset'])+'*mV; w_adapt+='+str(neuron_params['b'])+'*pA')
                                 

    if return_equations:
        return neurons, eqs
    else:
        return neurons

        
if __name__=='__main__':

    print(__doc__)
    
    # starting from an example

    from brian2 import *
    from library import get_neuron_params

    
    fig, AX = plt.subplots(1, 3, figsize=(10,3))
    for model, ax in zip(['LIF', 'EIF', 'AdExp'], AX):
        neurons, eqs =  get_membrane_equation(get_neuron_params(model), [],\
                                              return_equations=True)
        print('------------- NEURON model :', model)
        print(eqs)
        # V value initialization
        neurons.V = -70.*mV
        trace = StateMonitor(neurons, 'V', record=0)
        spikes = SpikeMonitor(neurons)
        run(100 * ms)
        neurons.I0 = 250*pA
        run(500 * ms)
        neurons.I0 = 0*pA
        run(100 * ms)
        # We draw nicer spikes
        V = trace[0].V[:]
        for t in spikes.t: V[int(t/defaultclock.dt)] = 20*mV
        ax.plot(trace.t / ms, V / mV, 'k')
        ax.set_title(model)
    xlabel('time (ms)')
    ylabel('membrane potential (mV)')
    show()

    

    
