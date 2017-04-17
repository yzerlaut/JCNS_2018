import os
import numpy as np

N = 2
freqs = np.logspace(-1, np.log(350)/np.log(10), N)
seeds = np.arange(4)
tstop = 2000.
amp = 5.
t0 = 500.

if __name__=='__main__':
    baseline_cmd = 'python oscillating_input.py -S --amp '+str(amp)+' --tstop '+str(tstop)
    for freq in freqs:
        for seed in seeds:
            filename = 'data/varying_freq_'+str(freq)+'_seed_'+str(seed)+'.npy'
            os.system(baseline_cmd+' --freq '+str(freq)+' --SEED '+str(seed)+' --file '+filename)
