import numpy as np
import matplotlib.pylab as plt
import sys

from model import Euler_method_for_ring_model
############ WE RUN IT OF MAIN FILE

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
    args = parser.parse_args()
    
    
    ### default configuration !
    print 'simulation [...]'
    t, X, Fe_aff, Fe, Fi, muVn = Euler_method_for_ring_model(\
                                                             args.NRN1, args.NRN2,\
                                                             args.NTWK, args.RING, args.STIM)
    print 'done !, now plotting...'
    ax, fig1 = space_time_vsd_style_plot(t*1e3, Fe_aff,\
                                         title='$\\nu_e^{aff}(x, t)$',\
                                         xlabel='time (ms)')
    ax, fig2 = space_time_vsd_style_plot(t*1e3, .8*Fe+.2*Fi,\
                                         title='$\\nu(x, t)$',\
                                         xlabel='time (ms)')
    ax, fig3 = space_time_vsd_style_plot(t*1e3, muVn,\
                                         xlabel='time (ms)', title='$V_m$ (x, t)',\
                                         zlabel='(mV)')
    
    # ax, fig1 = space_time_vsd_style_plot(t*1e3, muVn,\
    #                                      title='$\\nu_e^{aff}(x, t)$',\
    #                                      xlabel='time (ms)', with_latency_analysis=True)
    # ax, fig2 = space_time_vsd_style_plot(t*1e3, nu1,\
    #                title='$\\nu(x, t)$',\
    #                xlabel='time (ms)',params=params)
    # ax, fig3 = space_time_vsd_style_plot(t*1e3, 1e3*vm1,\
    #                xlabel='time (ms)', title='$V_m$ (x, t)',\
    #                zlabel='(mV)',params=params)
    # ax, fig4 = space_time_vsd_style_plot(t*1e3, vsd_signal1,\
    #                xlabel='time (ms)', title='vsd-like signal',\
    #                zlabel='$\|(V - \mu_V^\mathrm{rest}) / \mu_V^\mathrm{rest} \| (\%) \sim \Delta F / F$',\
    #                with_latency_analysis=True,params=params)

    # if args.SAVE:
    #     import sys
    #     sys.path.append('/home/yann/work/yns_python_libray/')
    #     from my_graph import put_list_of_figs_to_svg_fig
    #     put_list_of_figs_to_svg_fig([fig1, fig2, fig3, fig4], visualize=False)
    #     print 'saved as fig.svg'
    # else:
    #     plt.show()
    plt.show()
