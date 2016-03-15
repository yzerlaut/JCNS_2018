import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../code')
from my_graph import set_plot
import matplotlib

from theoretical_tools import *

def make_exc_inh_fig(DATA, P=None):
    
    MEANfreq, SDfreq, Fe_eff, fiSim, params = np.load(DATA)
    Fe_eff, Fout = np.array(Fe_eff), np.array(MEANfreq)
    fiSim = np.meshgrid(np.zeros(Fe_eff.shape[1]), fiSim)[1]
    levels = np.unique(fiSim) # to store for colors
    
    if P is not None:
        params['P']=P
        
    # # #### FIGURE AND COLOR GRADIENT STUFF
    
    fig1 = plt.figure(figsize=(6,4))
    plt.subplots_adjust(bottom=.2, left=.15, right=.85)
    # -- Setting up a colormap that's a simple transtion
    mymap = matplotlib.colors.LinearSegmentedColormap.from_list(\
                    'mycolors',['blue','red'])
    # -- Using contourf to provide my colorbar info, then clear the figure
    Z = [[0,0],[0,0]]
    CS3 = plt.contourf(Z, levels, cmap=mymap)
    plt.clf()

    ax = plt.subplot(111)
    cb = plt.colorbar(CS3,use_gridspec=True) ## TO BE ADDED
    cb.set_label('$\\nu_i$ inh. freq. (Hz)')
    
    for i in range(levels.size):

        SIMvector = MEANfreq[i][:]
        SDvector = SDfreq[i][:]
        feSim = Fe_eff[i][:]
        feth = np.linspace(feSim.min(), feSim.max(), 1e2)
        fi = fiSim[i][0]

        r = (float(levels[i])-levels.min())/(levels.max()-levels.min())
        ax.errorbar(feSim, SIMvector, yerr=SDvector,\
                    color=mymap(r,1),marker='D',ms=5, capsize=3, elinewidth=1, lw=0)
        if params.has_key('P'):
            Fout_th = TF_my_template(feth, fi, *pseq_params(params))
            ax.plot(feth, Fout_th, color=mymap(r,1), lw=5, alpha=.5)

    set_plot(plt.gca(), ['bottom', 'left'], xlabel='$\\nu_e$ exc. freq. (Hz)',\
             ylabel='$\\nu_{out}$   output. freq. (Hz)')

if __name__=='__main__':

    # First a nice documentation 
    parser=argparse.ArgumentParser(description=
     """ Runs two types of protocols on a given neuronal and network model
        1)  ==> Preliminary transfer function protocol ===
           to find the fixed point (with possibility to add external drive)
        2)  =====> Full transfer function protocol ==== 
           i.e. scanning the (fe,fi) space and getting the output frequency""",
              formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument("Protocol",help="Two types of protocols : PRE or FULL")
    # parser.add_argument("Neuron_Model",help="Choose a neuronal model from 'neuronal_models.py'")
    # parser.add_argument("Network_Model",help="Choose a network model (synaptic and connectivity properties)"+\
    #                     "\n      from 'network_models'.py")
    
    args = parser.parse_args()

    # if args.Protocol=='PRE':
    #     make_fixed_point_plot('data/preliminary_tf.npy')
    # else:
    make_exc_inh_fig('data/example_data.npy')
    plt.show()
