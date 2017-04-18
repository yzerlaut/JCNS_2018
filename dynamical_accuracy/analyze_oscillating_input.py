import matplotlib.pylab as plt
import numpy as np

from multiple_freq import * # all parameters inside !!
from oscillating_input import *
from scipy.optimize import minimize

modulus, phase_shift = 0*freqs, 0*freqs
exp_modulus, exp_phase_shift = 0*freqs, 0*freqs
Vexp_modulus, Vexp_phase_shift = 0*freqs, 0*freqs


fig1, ax1 = plt.subplots(figsize=(3,2))
plt.subplots_adjust(bottom=.25, left=.25)
fig2, ax2 = plt.subplots(figsize=(3,2))
plt.subplots_adjust(bottom=.25, left=.25)

for i in range(len(freqs)):
    # loading_data
    mod, phase = [], []
    for seed in seeds:
        filename = 'data/varying_freq_'+str(freqs[i])+'_seed_'+str(seed)+'.npy'
        args, t, rate_input, fe, fi = np.load(filename)

        ## fitting the response
        _, f_fitted, coeffs = find_modulus_and_phase_shift(1e-3*t, .8*fe+.2*fi,\
                                                           1e-3*args.t0, args.freq,\
                                                           full_output=True)
        mod.append(np.abs(coeffs[2]))
        phase.append((-coeffs[1]+np.pi/2.)%(2.*np.pi)-np.pi/2.)
    
    exp_modulus[i] = np.array(mod).mean()
    exp_phase_shift[i] = np.array(phase).mean()
    Vexp_modulus[i] = np.array(mod).std()
    Vexp_phase_shift[i] = np.array(phase).std()

    ## theoretical estimate
    def rate_func(t):
        return sinewave(t, 1e-3*args.t0, args.freq, args.amp)
    
    t_th, fe_th, fi_th = run_mean_field(args.CONFIG.split('--')[0],\
             args.CONFIG.split('--')[1],args.CONFIG.split('--')[2],\
             rate_func, tstop=max([2*args.t0, args.t0+1e3*5/args.freq])*1e-3, dt=5e-5)
    
    _, f_fitted, coeffs = find_modulus_and_phase_shift(t_th, .8*fe_th+.2*fi_th,\
                                                       1e-3*args.t0, args.freq,\
                                                       full_output=True, base_amp0=30./args.freq)

    phase_shift[i] = (-coeffs[1]+np.pi)%(2.*np.pi)-np.pi
    modulus[i] = np.abs(coeffs[2])

ax1.errorbar(freqs, exp_modulus, yerr=Vexp_modulus, color='k', fmt='.')
ax1.plot(freqs[:-1], modulus[:-1], 'k-', lw=3, alpha=.5)
ax1.set_xscale('log')
ax1.set_yscale('log')
set_plot(ax1, xticks=[1, 10, 100], xticks_labels=['1', '10', '100'],\
         ylabel='amplitude (Hz)', yticks=[1,3,10,30], yticks_labels=['1', '3', '10', '30'])

ax2.errorbar(freqs, exp_phase_shift, color='k', yerr=Vexp_phase_shift, fmt='.')
ax2.plot(freqs, phase_shift, 'k-', lw=3, alpha=.5)
ax2.set_xscale('log')
set_plot(ax2, yticks=[0, np.pi/2., np.pi], yticks_labels=['0', '$\pi$/2', '$\pi$'],\
         xticks=[1, 10, 100], xticks_labels=['1', '10', '100'],\
         ylabel='phase shift (Rd)', xlabel='freq. (Hz)')

plt.show()
# from my_graph import put_list_of_figs_to_svg_fig
# put_list_of_figs_to_svg_fig([fig1, fig2], visualize=False)
