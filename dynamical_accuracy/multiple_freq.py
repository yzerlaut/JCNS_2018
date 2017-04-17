import os
import numpy as np

N = 10
freqs = np.logspace(-1, np.log(350)/np.log(10), N)
seeds = np.arange(3)
amp = 5.
t0 = 500.

if __name__=='__main__':
    for freq in freqs:
        ts = max([2*t0, t0+1e3*5/freqs])
        baseline_cmd = 'python oscillating_input.py -S --amp '+str(amp)+' --tstop '+str(ts)
        for seed in seeds:
            filename = 'data/varying_freq_'+str(freq)+'_seed_'+str(seed)+'.npy'
            os.system(baseline_cmd+' --freq '+str(freq)+' --SEED '+str(seed)+' --file '+filename)
