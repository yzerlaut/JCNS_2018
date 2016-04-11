from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../code/')
from my_graph import set_plot


fig, AX = plt.subplots(2, 2)

f = loadmat('Amp_4yann.mat')
AX[0,0].plot(1e3*f['amp1_time'])
AX[0,0].set_title('Stim 1')
set_plot(AX[0,0], ['left'], ylabel=r'amp. ($\perthousand$)', xticks=[])
AX[0,1].plot(1e3*f['amp2_time'])
AX[0,1].set_title('Stim 2')
set_plot(AX[0,1], ['left'], ylabel=r'amp. ($\perthousand$)', xticks=[], ylim=[0,1.8])

f = loadmat('sigma_4yann.mat')
AX[1,0].plot(f['sigma1_time'])
set_plot(AX[1,0], ylabel='$\sigma$ (mm)', xlabel='time (frames)')
AX[1,1].plot(f['sigma2_time'])
set_plot(AX[1,1], ylabel='$\sigma$ (mm)', xlabel='time (frames)')

plt.show()



