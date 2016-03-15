from brian2 import *


def set_stim_from_time_varying_rate(time_array, rate_array, Neuron_Group, PRE='Gee_post += w',\
                                    seed=1):
    
    np.random.seed(1) # setting the seed !
    
    ## time_array in ms !!
    # so multplying rate array
    rate_array *= 1e-3
    
    indices, times = [], []
    DT = (time_array[1]-time_array[0])
    
    # trivial way to generate inhomogeneous poisson events
    for it in range(len(time_array)):
        rdm_num = np.random.random(Neuron_Group.N)
        for ii in np.arange(Neuron_Group.N)[rdm_num<DT*rate_array[it]]:
            indices.append(ii) # all the indicces
            times.append(time_array[it]) # all the same time !

    indices, times = array(indices), array(times)*ms
    
    input_spikes = SpikeGeneratorGroup(Neuron_Group.N, indices, times)
    feedforward = Synapses(input_spikes, Neuron_Group,\
                           pre=PRE,
                           model='w:siemens', connect='i==j')
    return input_spikes, feedforward


def build_up_excitatory_feedforward_connections_for_2_pop(Pops, syn_conn_matrix,\
                                                          time_array, rate_array):

    exc_neurons, inh_neurons = Pops
    P = syn_conn_matrix

    input_exc, fdfrwd_to_exc = set_stim_from_time_varying_rate(time_array, rate_array, exc_neurons,\
                                                    PRE='Gee_post += w')
    fdfrwd_to_exc.w = P[0,0]['Q']*nS

    input_inh, fdfrwd_to_inh = set_stim_from_time_varying_rate(time_array, rate_array, inh_neurons,\
                                                    PRE='Gei_post += w')
    fdfrwd_to_inh.w = P[0,1]['Q']*nS
    
    return input_exc, fdfrwd_to_exc, input_inh, fdfrwd_to_inh


def build_up_inhibitory_feedforward_connections_for_2_pop(Pops, syn_conn_matrix,\
                                                          time_array, rate_array):

    exc_neurons, inh_neurons = Pops
    P = syn_conn_matrix

    inh_input_exc, inh_fdfrwd_to_exc = set_stim_from_time_varying_rate(time_array, rate_array, exc_neurons,\
                                                    PRE='Gie_post += w')
    inh_fdfrwd_to_exc.w = P[1,0]['Q']*nS

    inh_input_inh, inh_fdfrwd_to_inh = set_stim_from_time_varying_rate(time_array, rate_array, inh_neurons,\
                                                    PRE='Gii_post += w')
    inh_fdfrwd_to_inh.w = P[1,1]['Q']*nS
    
    return inh_input_exc, inh_fdfrwd_to_exc, inh_input_inh, inh_fdfrwd_to_inh
