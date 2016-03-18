import numpy as np
import matplotlib.pylab as plt
import sys

import neuronal_models as models

## comparison with python simulations
from single_neuron_simulations import *
t = np.arange(0,700,1e-2)*1e-3


def make_model_figure(MODEL, I0=300e-12,
                      savefig=False, for_title=None):

    Inorm = np.array([1. if (tt>100e-3 and tt<500e-3) else 0 for tt in t])
    I = I0*Inorm
    # so that the 30 mV bar corresponds to I0
    
    params = models.get_model_params(MODEL)
    params['Ee'], params['Ei'] = 0, 0 # not used anyway

    v1, spikes = adexp_sim(t, I, 0*I, 0*I, *pseq_adexp(params))
    v2, spikes2 = adexp_sim(t, -I, 0*I, 0*I, *pseq_adexp(params))

    fig = plt.figure(figsize=(6,4))
    ax = plt.subplot(111, frameon=False)
    plt.title(params['name']+'\n', fontsize=25)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.plot(1e3*t, 1e3*v1, 'k')
    plt.plot(1e3*t, -85+10*Inorm, 'k')

    for s in spikes:
        plt.plot([1e3*s,1e3*s], [5e3*params['delta_v']+\
                                     1e3*params['Vthre'],20], 'k:')
    for s in spikes2:
        plt.plot([1e3*s,1e3*s], [5e3*params['delta_v']+\
                                     1e3*params['Vthre'],20], 'k:')
    plt.plot(1e3*t, 1e3*v2, 'k--')
    plt.plot(1e3*t, -85-10*Inorm, 'k--')
    plt.tight_layout()
    plt.plot([10,10],[-25,-15], 'gray', lw=3)
    plt.plot([10,60],[-25,-25], 'gray', lw=3)
    plt.annotate('10mV', (16,-10), textcoords='data', size=13)
    plt.annotate(str(int(1e12*I0))+'pA', (16,-20), textcoords='data',size=13)
    plt.annotate('50ms', (17,-40), textcoords='data', size=13)
    if savefig==True:
        fig.savefig('../figures/'+MODEL+'_step_response.svg',\
                    format='svg', transparent=True)
    return fig
    
if len(sys.argv)>2:
    fig = make_model_figure(sys.argv[1], I0=300e-12, savefig=True)
else:
    fig = make_model_figure(sys.argv[1], I0=300e-12, savefig=False)
    plt.show()



