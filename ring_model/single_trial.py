import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../code/')
from my_graph import put_list_of_figs_to_svg_fig

from model import Euler_method_for_ring_model
from ring_model.ring_models import pixels_per_mm


from plotting_tools import space_time_vsd_style_plot

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--NRN1",help="Choose a cell model", default='RS-cell')
    parser.add_argument("--NRN2",help="Choose a cell model", default='FS-cell')
    parser.add_argument("--NTWK",help="Choose a network model", default='CONFIG1')
    parser.add_argument("--RING",help="Choose a ring model", default='RING1')
    parser.add_argument("--STIM",help="Choose a network model", default='CENTER')
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    parser.add_argument("--no_sim", help="plot only", action="store_true")
    parser.add_argument("-f", "--file",help="filename for saving", default='data/example_data.npy')
    
    args = parser.parse_args()
    
    if not args.no_sim:
        print('simulation [...]')
        t, X, Fe_aff, Fe, Fi, muVn = Euler_method_for_ring_model(\
                                                                 args.NRN1, args.NRN2,\
                                                                 args.NTWK, args.RING, args.STIM)
        np.save(args.file, [args, t, X, Fe_aff, Fe, Fi, muVn])
        args2 = args
    else:
        args2, t, X, Fe_aff, Fe, Fi, muVn = np.load(args.file) # we just load a file
        
    params = {'pixels_per_mm':pixels_per_mm(args2.RING)}

    ax, fig1 = space_time_vsd_style_plot(t*1e3, Fe_aff,\
                                         title='$\\nu_e^{aff}(x, t)$',\
                                         params=params,
                                         xlabel='time (ms)', with_latency_analysis=True)
    ax, fig2 = space_time_vsd_style_plot(t*1e3, .8*Fe+.2*Fi,\
                                         title='$\\nu(x, t)$',\
                                         params=params,
                                         xlabel='time (ms)', with_latency_analysis=True)
    ax, fig3 = space_time_vsd_style_plot(t*1e3, 1e2*muVn,\
                                         xlabel='time (ms)', title='$\delta V / V_0 (x, t)$',\
                                         params=params,
                                         zlabel='%', with_latency_analysis=True)
    if args.SAVE:
        put_list_of_figs_to_svg_fig([fig1, fig2, fig3], visualize=False)
        for i in range(1,4):
            exec("fig"+str(i)+".savefig('fig"+str(i)+".png', dpi=300)")
            # exec("fig"+str(i)+".savefig('fig"+str(i)+".svg')")
    else:
        plt.show()
