"""
Some configuration of neuronal properties so that we pick up
within this file
"""
from __future__ import print_function
import numpy as np


def get_connectivity_and_synapses_matrix(NAME, number=2):

    exc_pop = {'p_conn':0.02, 'Q':7., 'Tsyn':5., 'Erev':0.}
    inh_pop = {'p_conn':0.02, 'Q':67., 'Tsyn':10., 'Erev':-80.}

    # creating empty arry of objects (future dictionnaries)
    M = np.empty((number, number), dtype=object)

    if NAME=='CONFIG1':
        M[:,0] = [exc_pop.copy(), inh_pop.copy()] # post-synaptic : exc
        M[:,1] = [exc_pop.copy(), inh_pop.copy()] # post-synaptic : inh
        M[0,0]['name'], M[1,0]['name'] = 'ee', 'ie'
        M[0,1]['name'], M[1,1]['name'] = 'ei', 'ii'
        
    return M

if __name__=='__main__':

    print(__doc__)

    M = get_connectivity_and_synapses_matrix('CONFIG1')

    print('synapses of the exc. pop. (pop. 0) : M[:,0]')
    print(M[:,0])
    print('synapses of the inh. pop. (pop. 1) : M[:,1]')
    print(M[:,1])
    
