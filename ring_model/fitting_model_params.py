import numpy as np
import matplotlib.pylab as plt
import itertools
# everything stored within a zip file
import zipfile, sys
sys.path.append("../experimental_data")
from dataset import get_dataset
from compare_to_model import get_data, get_residual
from model import Euler_method_for_ring_model
from scipy.optimize import minimize

def run_sim(X, args):
    t, X, Fe_aff, Fe, Fi, muVn =\
                                 Euler_method_for_ring_model(\
                                                             'RS-cell', 'FS-cell',\
                                                             'CONFIG1', 'RING1', 'CENTER',\
                                        custom_ring_params={\
                                                            'X_discretization':args.X_discretization,
                                                            'X_extent':args.X_extent,
                                                            'conduction_velocity_mm_s':X[0],
                                                            'exc_connect_extent':X[1],
                                                            'inh_connect_extent':X[1]/5.},
                                        custom_stim_params={\
                                                            'sX':X[4], 'amp':15.,
                                                            'Tau1':X[2], 'Tau2':X[3]})
    np.save('../ring_model/data/temp.npy', [args, t, X, Fe_aff, Fe, Fi, muVn])
                                                                    

                                                     
def run_fitting(args):
    
    ## loading data
    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)

    def to_minimize(X):
        """ X are the parameters """
        run_sim(X, args)
        return get_residual(args,
                           new_time, space, new_data,
                           Nsmooth=args.Nsmooth,
                           fn='../ring_model/data/temp.npy')
    BOUNDS = [args.vc, args.Econn_radius, args.Tau1, args.Tau2, args.stim_extent]
    X0 = np.array([args.vc, args.Econn_radius, args.Tau1, args.Tau2, args.stim_extent]).mean(axis=1)
    
    res = minimize(to_minimize, method='L-BFGS-B',
             x0=X0,
             bounds=BOUNDS,
             options={'maxiter':args.N})
        

    print(res)

    np.save('../ring_model/data/analyzed_scan_data_'+str(args.data_index)+'.npy',
            res.x)

def plot_analysis(args):

    X = np.load('../ring_model/data/analyzed_scan_data_'+str(args.data_index)+'.npy')

    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)
    # run_sim(X, args)
    res = get_residual(args,
                       new_time, space, new_data,
                       Nsmooth=args.Nsmooth,
                       fn='../ring_model/data/temp.npy', with_plot=True)
    plt.show()


def full_analysis(args):

    DATA = get_dataset()
    for i in range(len(DATA)):
        args.data_index = i
        run_fitting(args)


# def full_plot(args):

#     DATA = get_dataset()
#     VC, ECR, TAU2, TAU1, DUR = [], [], [], [], []
#     for i in range(len(DATA)):
#         args.data_index = i
#         params = get_minimum_params(args)
#         for vec, VEC in zip(params, [VC, ECR, TAU2, TAU1]):
#             VEC.append(vec)
#         DUR.append(DATA[i]['duration'])

#     fig, AX = plt.subplots(1, 3, figsize=(4.8,2.3))
#     plt.subplots_adjust(bottom=.3, left=.25, wspace=3.)
#     for ax, vec, label, ylim in zip(AX, [VC, ECR, np.array(ECR)/5.],
#                                     ['$v_c$ (mm/s)', '$s_{exc}$ (mm)', '$s_{inh}$ (mm)'],
#                                     [[0,500], [0,6], [0,6]]):
#         ax.plot([0, 0], ylim, 'w.', ms=0.1)
#         ax.bar([0], [np.array(vec).mean()], yerr=[np.array(vec).std()],
#                color='lightgray', edgecolor='k', lw=3)
#         set_plot(ax, ['left'], xticks=[], ylabel=label)

#     fig2, AX = plt.subplots(1, 2, figsize=(3.2,2.3))
#     plt.subplots_adjust(bottom=.3, left=.3, wspace=1.8)
#     AX[0].plot(DUR, 1e3*np.array(TAU1), 'o')
#     AX[0].plot(DUR, 1e3*np.array(TAU2), 'o')
#     plt.show()
    
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--vc", nargs=2, type=float, default=[50., 600.])
    parser.add_argument("--stim_extent", nargs=2, type=float, default=[0.1, 4.])
    parser.add_argument("--Econn_radius", nargs=2, type=float, default=[1., 10.])
    parser.add_argument("--ratio_Iconn_Econn", nargs=2, type=float, default=[0.05, 1.5])
    parser.add_argument("--Tau1", nargs=2, type=float, default=[5e-3, 100e-3])
    parser.add_argument("--Tau2", nargs=2, type=float, default=[50e-3, 400e-3])
    parser.add_argument("--N", type=int, default=20)
    parser.add_argument("--X_discretization", type=int, default=30) # PUT 100 for HD
    parser.add_argument("--X_extent", type=float, default=36.)
    parser.add_argument("--zip_filename", '-f', type=str, default='../ring_model/data/data.zip')
    # data
    parser.add_argument("--data_index", '-df', type=int,
                        default=7)
    parser.add_argument("--t0", type=float, default=-50.)
    parser.add_argument("--t1", type=float, default=200.)
    parser.add_argument("--Nsmooth", help="for data plots", type=int, default=2)
    # script function
    parser.add_argument("--fitting", help="fitting", action="store_true")
    parser.add_argument("-p", "--plot", help="plot analysis", action="store_true")
    parser.add_argument("-z", "--zip", help="zip datafiles", action="store_true")
    parser.add_argument("-uz", "--unzip", help="unzip datafiles", action="store_true")
    parser.add_argument("--full", help="full analysis", action="store_true")
    parser.add_argument("--full_plot", help="plot of full analysis", action="store_true")
    
    args = parser.parse_args()
    if args.fitting:
        run_fitting(args)
    elif args.plot:
        plot_analysis(args)
    # elif args.zip:
    #     zip_data(args)
    # elif args.unzip:
    #     unzip_data(args)
    # elif args.full:
    #     full_analysis(args)
    # elif args.full_plot:
    #     full_plot(args)
    # else:
    #     create_grid_scan_bash_script(args)
