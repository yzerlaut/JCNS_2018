import sys
sys.path.append('../')
import numpy as np
import matplotlib.pylab as plt

from ring_model.plotting_tools import space_time_vsd_style_plot as vsd_plot
from ring_model.model import Euler_method_for_ring_model
from ring_model.ring_models import pixels_per_mm

sys.path.append('../code/')
from my_graph import get_linear_colormap

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
    parser.add_argument("--yzoom",help="zoom on the y axis", nargs=2, default=[8.,26.], type=float)
    parser.add_argument("--yzoom_suppression",help="zoom on the y axis", nargs=2, default=[8.,26.], type=float)
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    parser.add_argument("--no_sim", help="plot only", action="store_true")
    parser.add_argument("-f", "--file",help="filename for saving", default='data/example_data.npy')
    
    args = parser.parse_args()
    # we perform one experiment with the default and we store the figs

    if not args.no_sim:
        print 'simulation [...]'
        print '====================== FIRST STIM simulation [...]'
        t, X, Fe_aff1, Fe1, Fi1, muVn1 = Euler_method_for_ring_model(\
                                                                 args.NRN1, args.NRN2,\
                                                                 args.NTWK, args.RING, 'FIRST_STIM')
        print '====================== SECOND STIM simulation [...]'
        t, X, Fe_aff2, Fe2, Fi2, muVn2 = Euler_method_for_ring_model(\
                                                                 args.NRN1, args.NRN2,\
                                                                 args.NTWK, args.RING, 'SECOND_STIM')
        print '====================== APPARENT MOTION simulation [...]'
        t, X, Fe_aff3, Fe3, Fi3, muVn3 = Euler_method_for_ring_model(\
                                                                 args.NRN1, args.NRN2,\
                                                                 args.NTWK, args.RING, 'AM')
        np.save(args.file, [args, t, X, Fe_aff1, Fe1, Fi1, muVn1,\
                            Fe_aff2, Fe2, Fi2, muVn2, Fe_aff3, Fe3, Fi3, muVn3])
        args2 = args
    else:
        args2, t, X, Fe_aff1, Fe1, Fi1, muVn1,\
          Fe_aff2, Fe2, Fi2, muVn2, Fe_aff3, Fe3, Fi3, muVn3 = np.load(args.file)

    vsd_label = '%'

    F1, F2, F3 = .8*Fe1+.2*Fi1, .8*Fe2+.2*Fi2, .8*Fe3+.2*Fi3
    
    suppression = muVn3-(muVn1+muVn2)
    
    zlim_vsd = [0, 1e2*np.max(np.abs(muVn1)+np.abs(muVn2))]
    
    params = {'pixels_per_mm':pixels_per_mm(args2.RING)}
    
    FIGS = []
    ax, fig = vsd_plot(t*1e3,\
                       1e2*suppression,\
                       zlim=[-1e2*np.abs(suppression).max(),1e2*np.abs(suppression).max()],
                       title=r'suppression signal',\
                       zlabel=vsd_label, with_latency_analysis=True,\
                       params=params,
                       xzoom=args.xzoom_suppression, yzoom=args.yzoom_suppression)
    FIGS.append(fig)
    
    ax, fig = vsd_plot(t*1e3,\
                       1e2*(np.abs(muVn1)+np.abs(muVn2)),\
                       zlim=zlim_vsd,
                       title=r'linear prediction',\
                       zlabel=vsd_label, with_latency_analysis=True,\
                       params=params,
                       xzoom=args.xzoom_suppression, yzoom=args.yzoom_suppression)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, Fe_aff1,\
                   title='$\\nu_e^{aff}(x, t)$',\
                       params=params,
                       xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, F1,\
                   title='$\\nu(x, t)$', xlabel='time (ms)',\
                       params=params,
                  xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, 1e2*muVn1,\
                       zlim=zlim_vsd,
                   title='vsd-like signal',zlabel=vsd_label,\
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, Fe_aff2,\
                   title='$\\nu_e^{aff}(x, t)$',\
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    ax, fig = vsd_plot(t*1e3, F2,\
                   title='$\\nu(x, t)$', xlabel='time (ms)',\
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, 1e2*muVn2,\
                   title='vsd-like signal',zlabel=vsd_label,\
                       zlim=zlim_vsd,
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, Fe_aff3,\
                   title='$\\nu_e^{aff}(x, t)$',\
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)
    
    ax, fig = vsd_plot(t*1e3, F3,\
                   title='$\\nu(x, t)$',\
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    ax, fig = vsd_plot(t*1e3, 1e2*muVn3,\
                   title='vsd-like signal', zlabel=vsd_label,\
                       zlim=zlim_vsd,
                       params=params,
                   xzoom=args.xzoom, yzoom=args.yzoom)
    FIGS.append(fig)

    
    if args.SAVE:
        for i in range(1,len(FIGS)+1):
            exec("FIGS["+str(i-1)+"].savefig('fig"+str(i)+".png', dpi=300)")
    else:
        plt.show()













