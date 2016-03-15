"""
Some configuration of neuronal properties so that we pick up
within this file
"""
from __future__ import print_function

def get_neuron_params(NAME, name='', number=1):

    if NAME=='LIF':
        params = {'name':name, 'N':number,\
                  'Gl':10., 'Cm':150.,'Trefrac':5.,\
                  'El':-60., 'Vthre':-50., 'Vreset':-60., 'delta_v':0.,\
                  'a':0., 'b': 0., 'tauw':0.}
    elif NAME=='EIF':
        params = {'name':name, 'N':number,\
                  'Gl':10., 'Cm':150.,'Trefrac':5.,\
                  'El':-60., 'Vthre':-50., 'Vreset':-60., 'delta_v':2.,\
                  'a':0., 'b':0., 'tauw':0.}
    elif NAME=='AdExp':
        params = {'name':name, 'N':number,\
                  'Gl':10., 'Cm':150.,'Trefrac':5.,\
                  'El':-60., 'Vthre':-50., 'Vreset':-60., 'delta_v':2.,\
                  'a':4., 'b':20., 'tauw':0.}
        
    return params.copy()

if __name__=='__main__':

    print(__doc__)
