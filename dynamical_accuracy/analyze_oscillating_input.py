import matplotlib.pylab as plt
import numpy as np
from oscillating_input import *

from multiple_freq import * # all parameters inside !!
from scipy.optimize import minimize

modulus, phase_shift = [], []
exp_modulus, exp_phase_shift = [], []

for freq in freqs:


    ### MEAN FIELD MODEL 
    def rate_func(t):
        return sinewave(t, 500e-3, freq, amp)

    t_th, fe_th, fi_th = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1',\
                                        rate_func, tstop=tstop*1e-3)
    
    F = .8*fe_th+.2*fi_th
    F, t_th = F[t_th>1e-3*t0], t_th[t_th>1e-3*t0]
    
    def to_minimize(X):
        [amp, phase, amp0]  = X
        return np.mean((F-amp0-sinewave(t_th, 1e-3*t0, freq, amp, phase=phase))**2)
    res = minimize(to_minimize, [2., np.pi/2., 0.], method='Nelder-Mead')
    
    # plt.plot(t_th, F, 'k-', lw=3)
    # plt.plot(t_th, res.x[2]+sinewave(t_th, 1e-3*t0, freq, res.x[0], phase=res.x[1]), 'r-')

    phase_shift.append(res.x[1]%(2.*np.pi))
    modulus.append(np.abs(res.x[0]))

    ### NUMERICAL SIM
    t, rate_input, fe, fi, t_th, fe_th, fi_th = np.load('data/varying_freq_'+str(round(freq))+'.npy')
    F = .8*fe+.2*fi
    
    def to_minimize(X):
        [amp, phase, amp0]  = X
        return np.mean((F-amp0-sinewave(1e-3*t, 1e-3*t0, freq, amp, phase=phase))**2)
    res = minimize(to_minimize, [2., np.pi/2., 0.], method='Nelder-Mead')
    

    exp_phase_shift.append(res.x[1]%(2.*np.pi))
    exp_modulus.append(np.abs(res.x[0]))


plt.subplots(figsize=(5,3))
plt.semilogx(freqs, modulus, 'k-', lw=3, alpha=.5)
plt.semilogx(freqs, exp_modulus, 'kD')
set_plot(plt.gca())
plt.subplots(figsize=(5,3))
plt.semilogx(freqs, phase_shift, 'k-', lw=3, alpha=.5)
plt.semilogx(freqs, exp_phase_shift, 'kD')
set_plot(plt.gca())


plt.show()
