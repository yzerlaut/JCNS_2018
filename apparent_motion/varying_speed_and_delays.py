import sys
import numpy as np
from apparent_motion_protocol import compute_response_to_two_stimuli, get_vsd_signals_and_suppression


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs the apparent motion protocol
            """,
            formatter_class=argparse.RawTextHelpFormatter)
    # model params
    parser.add_argument("--CELL",help="Choose a cell model", default='cell1')
    parser.add_argument("--NTWK",help="Choose a network model", default='CONFIG2')
    parser.add_argument("--RING",help="Choose a ring model", default='RING1')
    parser.add_argument("--STIM",help="Choose a network model", default='STIM1')
    # graphics
    parser.add_argument("--xzoom",help="zoom on the x axis", nargs=2, default=[0.,400.], type=float)
    parser.add_argument("--yzoom",help="zoom on the y axis", nargs=2, default=[20.,30.], type=float)
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    args = parser.parse_args()

    params, t, X, vm0, nu0, w0, nu_e_aff1, vm1, nu1, w1, nu_e_aff2,\
        vm2, nu2, w2, nu_e_aff12, vm12, nu12, w12 = \
         compute_response_to_two_stimuli(args.CELL, args.NTWK, args.RING, args.STIM)
         
    vsd_signal1, vsd_signal2, vsd_signal12, suppression = \
      get_vsd_signals_and_suppression(vm0, vm1, vm2, vm12)
      
    suppression[0,0] = -suppression.min() #

    FIGS = []
    ax, fig = vsd_plot(t*1e3,\
                       suppression,\
                       title=r'suppression signal',\
                       phase_criteria=-np.pi/2.+np.pi/5.,
                       zlabel=vsd_label, with_latency_analysis=True,\
                       xzoom=args.xzoom_suppression, yzoom=args.yzoom, params=params)
