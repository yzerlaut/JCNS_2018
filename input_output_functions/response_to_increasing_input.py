import sys
sys.path.append('../')
from mean_field.euler_method import run_mean_field
from network_simulations.waveform_input import double_gaussian, smooth_heaviside

sys.path.append('../code')
from my_graph import set_plot

import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

def heaviside(x):
    return .5*(1+np.sign(x))
def step_input(t, T0, amp, T1=0.02):
    return amp*heaviside(t-T0)*(1-np.exp(-(t-T0)/T1))
    # return 0*amp*np.exp(-(t-T0)**2/2./T1**2)


t0, T1, T2, tstop = 250e-3, 50e-3, 70e-3, 500e-3
amp_max = 12
N=25

if sys.argv[-1]=='full':

    fig1, [ax1, ax2] = plt.subplots(2, figsize=(5,6))
    plt.subplots_adjust(left=.25, bottom=.25 )

    
    max_f_amp, max_vm_amp = np.zeros(N), np.zeros(N)
    amplitudes = np.linspace(0, amp_max, N)
    
    for i in range(N):
        amp = amplitudes[i]
        def func(t):
            # return step_input(t, 0.02, 1.)*amp
            return double_gaussian(t, t0, T1, T2, amp)
        t, fe, fi, muV, sV, muG, Tv = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3,\
                                                     ext_drive_change=amp,
                                                     afferent_exc_fraction=None,
                                                     extended_output=True, tstop=tstop)
        max_f_amp[i] = np.max(.8*fe+.2*fi)
        max_vm_amp[i] = np.max(1e2*np.abs((muV-muV[0])/muV[0]))

    ax1.plot(amplitudes, max_f_amp, 'k-', lw=3)
    pol = np.polyfit(amplitudes[:3], max_f_amp[:3], 1)
    ax1.plot(amplitudes[:int(2*N/3.)], np.polyval(pol, amplitudes[:int(2*N/3.)]), 'k--')
    set_plot(ax1, ['left'], ylabel='max. $\\nu$ (Hz)', xticks=[])
    ax2.plot(amplitudes, max_vm_amp, 'k-', lw=3)
    pol = np.polyfit(amplitudes[:3], max_vm_amp[:3], 1)
    ax2.plot(amplitudes[:int(2*N/3.)], np.polyval(pol, amplitudes[:int(2*N/3.)]), 'k--')
    set_plot(ax2, ylabel=r'$\| \delta V/V_0 \| $ %', xlabel='max. $\\nu_e^{aff}$ (Hz)')
    plt.show()

else:

    fig1, [ax1, ax2, ax3] = plt.subplots(3, figsize=(5,6))
    plt.subplots_adjust(left=.25, bottom=.25 )

    for amp in np.linspace(0, amp_max, 10):
        def func(t):
            # return step_input(t, 0.02, 1.)*amp
            return double_gaussian(t, t0, T1, T2, amp)

        t, fe, fi, muV, sV, muG, Tv = run_mean_field('RS-cell', 'FS-cell', 'CONFIG1', func, T=5e-3,\
                                                     ext_drive_change=amp,
                                                     afferent_exc_fraction=None,
                                                     extended_output=True, tstop=tstop)
        ax1.plot(1e3*t, func(t), 'k')
        ax2.plot(1e3*t, .8*fe+.2*fi, 'k')
        ax3.plot(1e3*t, 1e2*np.abs((muV-muV[0])/muV[0]), 'k')

    ax1.annotate('external input', (0,0))    
    set_plot(ax1, ['left'], ylabel='$\\nu_e^{aff}$', xticks=[])
    ax2.annotate('network response', (0,0))    
    set_plot(ax2, ['left'], ylabel='$\\nu$', xticks=[])
    ax3.annotate('mean membrane potential', (0,0))    
    set_plot(ax3, ylabel=r'$\| \delta V/V_0 \| $ %', xlabel='time (ms)')
    plt.show()

