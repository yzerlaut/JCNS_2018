import matplotlib.pylab as plt
import numpy as np
from oscillating_input import *

from multiple_freq import * # all parameters inside !!
from oscillating_input import *
from scipy.optimize import minimize

modulus, phase_shift = 0*freqs, 0*freqs
exp_modulus, exp_phase_shift = 0*freqs, 0*freqs


for i in range(len(freqs)):
    filename = 'data/varying_freq_'+str(freqs[i])+'.npy'
    args, t, rate_input, fe, fi = np.load(filename)

    ## fitting the response
    _, f_fitted, coeffs = find_modulus_and_phase_shift(1e-3*t, .8*fe+.2*fi,\
                                                       1e-3*args.t0, args.freq,\
                                                       full_output=True)
    
    exp_phase_shift[i] = (-coeffs[1])%(2.*np.pi)
    exp_modulus[i] = np.abs(coeffs[2])

    ## theoretical estimate
    def rate_func(t):
        return sinewave(t, 1e-3*args.t0, args.freq, args.amp)
    
    t_th, fe_th, fi_th = run_mean_field(args.CONFIG.split('--')[0],\
             args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
             rate_func, tstop=args.tstop*1e-3)
    
    _, f_fitted, coeffs = find_modulus_and_phase_shift(t_th, .8*fe_th+.2*fi_th,\
                                                       1e-3*args.t0, args.freq,\
                                                       full_output=True)
    
    phase_shift[i] = (-coeffs[1])%(2.*np.pi)
    modulus[i] = np.abs(coeffs[2])
    
plt.subplots(figsize=(5,3))
plt.semilogx(freqs, modulus, 'k-', lw=3, alpha=.5)
plt.semilogx(freqs, exp_modulus, 'kD')
set_plot(plt.gca(), xticks=[1, 10, 100], xticks_labels=['1', '10', '100'])
plt.subplots(figsize=(5,3))
plt.semilogx(freqs, phase_shift, 'k-', lw=3, alpha=.5)
plt.semilogx(freqs, exp_phase_shift, 'kD')
set_plot(plt.gca(), yticks=[0, np.pi, 2.*np.pi], yticks_labels=['0', '$\pi$', '2$\pi$'],\
         xticks=[1, 10, 100], xticks_labels=['1', '10', '100'])


plt.show()