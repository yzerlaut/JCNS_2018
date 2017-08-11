import sys
sys.path.append('../code')
from signanalysis import gaussian_func
from my_graph import set_plot
import matplotlib.pylab as plt
from scipy.optimize import minimize
import numpy as np

def heaviside(x):
    return 0.5*(1+np.sign(x))

def double_gaussian(t, t0, T1, T2, amplitude, amp0):
    return amp0+amplitude*(\
                           np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                           np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))

def simple_gaussian(X, X0, SX, amplitude, amp0):
    return amp0+amplitude*np.exp(-(X-X0)**2/2./SX**2)

def get_mean_extent_and_temporal_quant(t, X, Y, amp_criteria=0.02):

    ## temporal quantities
    XX, T1, T2 = [], [], []
    for i in range(Y.shape[1]):
        imax = np.argmax(Y[:,i])
        if Y[imax,i]>=amp_criteria*Y.max():
            def to_minimize(X):
                [t0, t1, t2, amp, amp0]  = X
                return np.mean((Y[:,i]-double_gaussian(t, t0, t1, t2, amp, amp0))**2)
            X0 = [20e-2, 50e-3, 100e-3, Y[:,i].max()-Y[0,i], Y[0,i]]
            res = minimize(to_minimize, X0, method='Nelder-Mead')
            XX.append(X[i])
            T1.append(res.x[1])
            T2.append(res.x[2])

    ## spatial quantities
    TT, AMP, SIGMA = [], [], []
    for i in range(Y.shape[0]):
        imax = np.argmax(Y[i,:])
        if Y[i, imax]>=amp_criteria*Y.max():
            def to_minimize(X0):
                [x0, sx, amp, amp0]  = X0
                return np.mean((Y[i,:]-simple_gaussian(X, x0, sx, amp, amp0))**2)
            X0 = [20, 5, Y[i,:].max()-Y[0,0], Y[0,0]]
            res = minimize(to_minimize, X0, method='Nelder-Mead')
            TT.append(t[i])
            AMP.append(res.x[2])
            SIGMA.append(res.x[1])

    
    return np.array(XX), np.array(T1), np.array(T2), np.array(TT), np.array(AMP), np.array(SIGMA)

    
if __name__=='__main__':

    args, t, X, Fe_aff, Fe, Fi, muVn = np.load('data/example_data.npy')

    XX, T1, T2, TT, AMP, SIGMA = \
            get_mean_extent_and_temporal_quant(t, X, muVn)

    fig, AX = plt.subplots(2, 2)
    plt.subplots_adjust(bottom=.2,wspace=.4)
    
    AX[0,0].set_title('Spatial profile in time \n')
    AX[0,0].plot(1e3*TT, 1e2*AMP, 'k-', lw=2)
    set_plot(AX[0,0], ['left'], ylabel=r'amp. %', xticks=[])
    AX[1,0].plot(1e3*TT, SIGMA, 'k-', lw=2)
    set_plot(AX[1,0], ylabel='$\sigma$ (mm)', xlabel='time (ms)')

    AX[0,1].set_title('Temporal profile in space \n')
    AX[0,1].plot(XX, 1e3*T1, 'k-', lw=2)
    set_plot(AX[0,1], ['left'], ylabel='$\\tau_1$ (ms)', xticks=[])
    AX[1,1].plot(XX, 1e3*T2, 'k-', lw=2)
    set_plot(AX[1,1], ylabel='$\\tau_2$ (ms)', xlabel='space (mm)')

    fig.savefig('fig.svg')
    plt.show()
