import sys
sys.path.append('../code/')
from my_graph import set_plot, put_list_of_figs_to_svg_fig
from signanalysis import gaussian_func

import numpy as np
import matplotlib.pylab as plt
sys.path.append('../')
from mean_field.master_equation import find_fixed_point
from transfer_functions.theoretical_tools import mean_and_var_conductance, get_fluct_regime_vars, pseq_params
from scipy.signal import gaussian
from single_cell_models.cell_library import get_neuron_params
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
from transfer_functions.tf_simulation import reformat_syn_parameters

def skew(x):
    return np.abs(np.mean((x.mean()-x)**3))**(1./3.)#/x.std()**3)

def plot_ntwk_sim_output(time_array, rate_array, rate_exc, rate_inh,\
                         Raster_exc, Raster_inh,\
                         Vm_exc, Vm_inh, Ge_exc, Ge_inh, Gi_exc, Gi_inh,\
                         BIN=5, min_time=200):
    
    
    cond_t = (time_array>min_time) # transient behavior after 400 ms

    params = get_neuron_params('RS-cell', SI_units=True)
    M = get_connectivity_and_synapses_matrix('CONFIG1', SI_units=True)
    EXC_AFF = M[0,0]['ext_drive']
    
    print 'starting fixed point'
    fe0, fi0, sfe, sfie, sfi = find_fixed_point('RS-cell', 'FS-cell', 'CONFIG1',\
                                                exc_aff=EXC_AFF, Ne=8000, Ni=2000, verbose=False)
    print 'end fixed point'
    
    reformat_syn_parameters(params, M) # merging those parameters
    
    xfe = fe0+np.linspace(-4,4)*sfe
    fe_pred = gaussian_func(xfe, fe0, sfe)
    xfi = fi0+np.linspace(-4,4)*sfi
    fi_pred = gaussian_func(xfi, fi0, sfi)

    mGe, mGi, sGe, sGi = mean_and_var_conductance(fe0+EXC_AFF, fi0, *pseq_params(params))
    muV, sV, muGn, TvN = get_fluct_regime_vars(fe0+EXC_AFF, fi0, *pseq_params(params))

    FE, FI = np.meshgrid(xfe, xfi)
    pFE, pFI = np.meshgrid(fe_pred, fi_pred)
    MUV, SV, _, _ = get_fluct_regime_vars(FE+EXC_AFF, FI, *pseq_params(params))*pFE*pFI/np.sum(pFE*pFI)
    ### MEMBRANE POTENTIAL
    MEAN_VM, STD_VM, KYRT_VM = [], [], []
    for i in range(len(Vm_exc)):
        MEAN_VM.append(Vm_exc[i][(time_array>min_time) & (Vm_exc[i]!=-65) & (Vm_exc[i]<-50)].mean())
        MEAN_VM.append(Vm_inh[i][(time_array>min_time) & (Vm_inh[i]!=-65) & (Vm_inh[i]<-50)].mean())
        for vv in [Vm_exc[i][(time_array>min_time)], Vm_inh[i][(time_array>min_time)]]:
            i0 = np.where((vv[:-1]>-52) & (vv[1:]<-60))[0]
            print(i0)
            sv = []
            if len(i0)==0:
                STD_VM.append(vv.std())
            elif len(i0)==1:
                STD_VM.append(vv[i0:].std())
            else:
                for i1, i2 in zip(i0[:-1], i0[1:]):
                    if i2-i1>60:
                        sv.append(vv[i1+30:i2-30].std())
                STD_VM.append(np.array(sv).mean())
        # STD_VM.append(Vm_inh[i][(time_array>min_time) & (Vm_inh[i]<-50)].std())

    fig1, AX1 = plt.subplots(1, 3, figsize=(3,2)) # for means
    fig2, AX2 = plt.subplots(1, 3, figsize=(3,2)) # for std
    
    AX1[0].bar([0], np.array(MEAN_VM).mean()+65, yerr=np.array(MEAN_VM).std(), color='w', edgecolor='k', lw=3, error_kw=dict(elinewidth=3,ecolor='k'))
    AX2[0].bar([0], np.array(STD_VM).mean(), yerr=np.array(STD_VM).std(), color='w', edgecolor='k', lw=3, error_kw=dict(elinewidth=3,ecolor='k'), label='$V_m$')
    AX1[0].bar([1], [1e3*muV+65], color='gray', alpha=.5, label='$V_m$')
    AX2[0].bar([1], [1e3*sV], color='gray', alpha=.5)

    set_plot(AX1[0], ['left'], xticks=[], ylim=[0,11], yticks=[0, 5, 10], yticks_labels=['-65', '-60', '-55'], ylabel='mean (mV)')
    set_plot(AX2[0], ['left'], xticks=[], ylim=[0,5], yticks=[0, 2, 4], ylabel='std. dev. (mV)')

    
    ### EXCITATORY CONDUCTANCE
    MEAN_GE, STD_GE, KYRT_GE = [], [], []
    for i in range(len(Ge_exc)):
        MEAN_GE.append(Ge_exc[i][(time_array>min_time)].mean())
        MEAN_GE.append(Ge_inh[i][(time_array>min_time)].mean())
        STD_GE.append(Ge_exc[i][(time_array>min_time)].std())
        STD_GE.append(Ge_inh[i][(time_array>min_time)].std())

    AX1[1].bar([0], np.array(MEAN_GE).mean(), yerr=np.array(MEAN_GE).std(), color='w', edgecolor='g', lw=3, error_kw=dict(elinewidth=3,ecolor='g'), label='num. sim.')
    AX2[1].bar([0], np.array(STD_GE).mean(), yerr=np.array(STD_GE).std(), color='w', edgecolor='g', lw=3, error_kw=dict(elinewidth=3,ecolor='g'), label='exc.')
    AX1[1].bar([1], [1e9*mGe], color='g', label='mean field \n pred.')
    AX2[1].bar([1], [1e9*sGe], color='g', label='exc.')
    
    set_plot(AX1[1], ['left'], xticks=[], yticks=[0,15,30], ylabel='mean (nS)')
    set_plot(AX2[1], ['left'], xticks=[], yticks=[0,5,10], ylabel='std. dev. (nS)')
    
    ### INHIBITORY CONDUCTANCE
    MEAN_GI, STD_GI, KYRT_GI = [], [], []
    for i in range(len(Gi_exc)):
        MEAN_GI.append(Gi_exc[i][(time_array>min_time)].mean())
        MEAN_GI.append(Gi_inh[i][(time_array>min_time)].mean())
        STD_GI.append(Gi_exc[i][(time_array>min_time)].std())
        STD_GI.append(Gi_inh[i][(time_array>min_time)].std())

    AX1[2].bar([0], np.array(MEAN_GI).mean(), yerr=np.array(MEAN_GI).std(), color='w', edgecolor='r', lw=3, error_kw=dict(elinewidth=3,ecolor='r'), label='num. sim.')
    AX2[2].bar([0], np.array(STD_GI).mean(), yerr=np.array(STD_GI).std(), color='w', edgecolor='r', lw=3, error_kw=dict(elinewidth=3,ecolor='r'), label='inh.')
    AX1[2].bar([1], [1e9*mGi], color='r', label='mean field \n pred.')
    AX2[2].bar([1], [1e9*sGi], color='r', label='inh.')
    
    set_plot(AX1[2], ['left'], xticks=[], yticks=[0,15,30], ylabel='mean (nS)')
    set_plot(AX2[2], ['left'], xticks=[], yticks=[0,5,10], ylabel='std. dev. (nS)')

    ### POPULATION RATE ###
    
    fig, ax = plt.subplots(figsize=(4,3))
    # we bin the population rate
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[cond_t][-1]-time_array[cond_t][0])/BIN)
    time_array = time_array[cond_t][:N0*N1].reshape((N1,N0)).mean(axis=1)
    rate_exc = rate_exc[cond_t][:N0*N1].reshape((N1,N0)).mean(axis=1)
    rate_inh = rate_inh[cond_t][:N0*N1].reshape((N1,N0)).mean(axis=1)
    rate_array = rate_array[cond_t][:N0*N1].reshape((N1,N0)).mean(axis=1)

    hh, bb = np.histogram(rate_exc, bins=8, normed=True)
    ax.bar(.5*(bb[:-1]+bb[1:]), hh, color='w', width=bb[1]-bb[0], edgecolor='g', lw=3, label='exc.', alpha=.7)
    hh, bb = np.histogram(rate_inh, bins=8, normed=True)
    ax.bar(.5*(bb[:-1]+bb[1:]), hh, color='w', width=bb[1]-bb[0], edgecolor='r', lw=3, label='inh', alpha=.7)

    ax.fill_between(xfe, 0*fe_pred, fe_pred, color='g')
    ax.fill_between(xfi, 0*fi_pred, fi_pred, color='r')

    set_plot(plt.gca(), ['bottom', 'left'], xlabel='pop. activity (Hz)', yticks=[], ylabel='density')

    for ax in AX1: ax.legend()
    for ax in AX2: ax.legend()
    
    return [fig, fig1, fig2]

if __name__=='__main__':


    import argparse
    # First a nice documentation 
    parser=argparse.ArgumentParser(description=
     """ 
     ----------------------------------------------------------------------
     Run the a network simulation using brian2

     Choose CELLULAR and NTWK PARAMETERS from the available libraries
     see  ../synapses_and_connectivity.syn_and_connec_library.py for the CELLS
     see ../synapses_and_connectivity.syn_and_connec_library.py for the NTWK

     Then construct the input as "NRN_exc--NRN_inh--NTWK"
     example: "LIF--LIF--Vogels-Abbott"
     ----------------------------------------------------------------------
     """
    ,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-f", "--file",help="filename for saving", default='data/example.npy')
    parser.add_argument("-z", "--zoom",help="zoom for activity", type=float, nargs=2)
    parser.add_argument("-b", "--bar_ms",help="bar for legend", type=int, default=100)
    parser.add_argument("-r", "--raster_number",help="max neuron number", type=int, default=10000)
    parser.add_argument("-t_after_transient", type=float, default=500)
    parser.add_argument("-s", "--save",action='store_true')

    args = parser.parse_args()
    
    time_array, rate_array, rate_exc, rate_inh,\
        Raster_exc, Raster_inh, Vm_exc, Vm_inh,\
        Ge_exc, Ge_inh, Gi_exc, Gi_inh = np.load(args.file)

    FIG = plot_ntwk_sim_output(time_array, rate_array, rate_exc, rate_inh,\
                               Raster_exc, Raster_inh,\
                               Vm_exc, Vm_inh, Ge_exc, Ge_inh, Gi_exc, Gi_inh,\
                               min_time=args.t_after_transient)

    if args.save:
        put_list_of_figs_to_svg_fig(FIG, visualize=False)
    else:
        plt.show()







        
