from scipy.io import loadmat
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../../')
from graphs.my_graph import set_plot
from scipy.signal import convolve2d
import matplotlib.cm as cm

def gaussian(x, m, s, w):
    y = np.exp(-(x-m)**2/2/s**2)
    return w*y

def heaviside(x):
    return .5*(1+np.sign(x))

def double_gaussian(x, m1, m2, s, w):
    y = heaviside(m1-x)*np.exp(-(x-m1)**2/2/s**2)+\
        heaviside(m2-x)*heaviside(x-m1)+\
        heaviside(x-m2)*np.exp(-(x-m2)**2/2/s**2)
    return w*y

f = loadmat('matx2_br131023_NL_dw_1020_23.mat')
data = 1e3*f['matNL'][0]['stim1'][0]
time = f['matNL'][0]['time'][0].flatten()+100
space = f['matNL'][0]['space'][0].flatten()
Nsmooth = 3
smoothing = np.ones((Nsmooth, Nsmooth))/Nsmooth**2
smooth_data = convolve2d(data, smoothing, mode='same')

cond = (time>-10) & (time<300)

fig, ax = plt.subplots(1, figsize=(4.5,3))
plt.subplots_adjust(bottom=.23, top=.97, right=.85)
c = ax.contourf(time[cond], space, smooth_data[:,cond], cmap=cm.viridis)
plt.colorbar(c, label='VSD signal ($\perthousand$)', ticks=-.3+.3*np.arange(5))
ax.plot([0,0], [0,2], '-', color='gray', lw=4)
ax.annotate('2mm', (0,2), rotation=90, fontsize=14)

from scipy.optimize import minimize

cond = (time>70) & (time<120)
mean_data = smooth_data[:,cond].mean(axis=1)

fig, AX = plt.subplots(1, 2, figsize=(8,3))
plt.subplots_adjust(bottom=.23, top=.97, right=.85)
def to_minimize(x):
    return np.sum((mean_data-double_gaussian(space, *x))**2)
res = minimize(to_minimize, [3, 6, 2, .9])
print(res.x)
AX[0].plot(space, double_gaussian(space, *res.x), label='double gaussian')

def to_minimize(x):
    return np.sum((mean_data-gaussian(space, *x))**2)
res = minimize(to_minimize, [5, 2, .9])
AX[1].plot(space, gaussian(space, *res.x), label='simple gaussian')
for ax in AX:
    ax.plot(space, mean_data, label='data [70,120] ms')
    ax.legend()
set_plot(AX[0], xlabel='space (mm)', ylabel='VSDi signal ($\perthousand$)')
set_plot(AX[1], xlabel='space (mm)')
fig.savefig('/Users/yzerlaut/Desktop/1.png')
plt.show()



