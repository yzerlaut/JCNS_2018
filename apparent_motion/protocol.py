import sys
sys.path.append('../')
import numpy as np
import matplotlib.pylab as plt

from ring_model.plotting_tools import space_time_vsd_style_plot as vsd_plot
from ring_model.model import Euler_method_for_ring_model


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs the apparent motion protocol
            """,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--NRN1",help="Choose a cell model", default='RS-cell')
    parser.add_argument("--NRN2",help="Choose a cell model", default='FS-cell')
    parser.add_argument("--NTWK",help="Choose a network model", default='CONFIG1')
    parser.add_argument("--RING",help="Choose a ring model", default='RING1')
    parser.add_argument("--xzoom",help="zoom on the x axis", nargs=2, default=[0.,400.], type=float)
    parser.add_argument("--xzoom_suppression",help="zoom on the x axis", nargs=2, default=[50.,300.], type=float)
    parser.add_argument("--yzoom",help="zoom on the y axis", nargs=2, default=[18.,32.], type=float)
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    args = parser.parse_args()
    # we perform one experiment with the default and we store the figs

    print 'FIRST STIM simulation [...]'
    t, X, Fe_aff1, Fe1, Fi1, muVn1 = Euler_method_for_ring_model(\
                                                             args.NRN1, args.NRN2,\
                                                             args.NTWK, args.RING, 'FIRST_STIM')
    print 'SECOND STIM simulation [...]'
    t, X, Fe_aff2, Fe2, Fi2, muVn2 = Euler_method_for_ring_model(\
                                                             args.NRN1, args.NRN2,\
                                                             args.NTWK, args.RING, 'SECOND_STIM')
    print 'APPARENT MOTION simulation [...]'
    t, X, Fe_aff3, Fe3, Fi3, muVn3 = Euler_method_for_ring_model(\
                                                             args.NRN1, args.NRN2,\
                                                             args.NTWK, args.RING, 'AM')

    vsd_label = r'$\|(\mu_V(x,t) - \mu_V^\mathrm{rest}) / \mu_V^\mathrm{rest} \| \sim \Delta F / F$ ($\perthousand$)'

    F1, F2, F3 = .8*Fe1+.2*Fi1, .8*Fe2+.2*Fi2, .8*Fe3+.2*Fi3
    
    suppression = muVn3-(muVn1+muVn2)
    suppression[0,0] = -suppression.min() #

    FIGS = []
    ax, fig = vsd_plot(t*1e3,\
                       suppression,\
                       title=r'suppression signal',\
                       phase_criteria=-np.pi/2.+np.pi/5.,
                       zlabel=vsd_label, with_latency_analysis=True,\
                       xzoom=args.xzoom_suppression, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, Fe_aff1,\
                   title='$\\nu_e^{aff}(x, t)$',\
                  xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, F1,\
                   title='$\\nu(x, t)$', xlabel='time (ms)',\
                  xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, muVn1,\
                   title='vsd-like signal',zlabel=vsd_label,\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, Fe_aff2,\
                   title='$\\nu_e^{aff}(x, t)$',\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, F2,\
                   title='$\\nu(x, t)$', xlabel='time (ms)',\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, muVn2,\
                   title='vsd-like signal',zlabel=vsd_label,\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, Fe_aff3,\
                   title='$\\nu_e^{aff}(x, t)$',\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    
    ax, fig = vsd_plot(t*1e3, F3,\
                   title='$\\nu(x, t)$',\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, muVn3,\
                   title='vsd-like signal', zlabel=vsd_label,\
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    if args.SAVE:
        sys.path.append('/home/yann/work/yns_python_libray/')
        from my_graph import put_list_of_figs_to_svg_fig
        put_list_of_figs_to_svg_fig(FIGS, visualize=False)
        print 'saved as fig.svg'
    else:
        plt.show()













