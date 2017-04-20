import os
import numpy as np

N = 20
freqs = np.logspace(0, np.log(300)/np.log(10), N)
seeds = np.arange(3)
amp = 4.5
t0 = 500.

if __name__=='__main__':
    for freq in freqs:
        ts = min([max([2*t0, t0+1e3*5/freq]), 10000])
        baseline_cmd = 'python oscillating_input.py -S --amp '+str(amp)+' --tstop '+str(ts)+' --DT 0.05'
        for seed in seeds:
            filename = 'data/varying_freq_'+str(freq)+'_seed_'+str(seed)+'.npy'
            os.system(baseline_cmd+' --freq '+str(freq)+' --SEED '+str(seed)+' --file '+filename)
