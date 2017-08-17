import numpy as np
import matplotlib.pylab as plt
import itertools
# everything stored within a zip file
import zipfile, sys
sys.path.append("../experimental_data")
from dataset import get_dataset
from compare_to_model import get_data, get_residual, get_time_residual, get_space_residual
from model import Euler_method_for_ring_model
from scipy.optimize import minimize
sys.path.append("../../")
from graphs.my_graph import set_plot

def run_sim(X, args, fn=None):
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
    if fn is None:
        # random filename for parrallel minimization        
        fn = str(int(\
                     np.random.randint(100000)*abs(X[0])*abs(X[1])))+'.npy' 
    np.save(fn, [args, t, X, Fe_aff, Fe, Fi, muVn])
    return fn
                                                                    
                                                     
def run_fitting(args):

    # baseline params and boundaries
    X0 = np.array([args.vc, args.Econn_radius, args.Tau1, args.Tau2, args.stim_extent]).mean(axis=1)
    BOUNDS = [args.vc, args.Econn_radius, args.Tau1, args.Tau2, args.stim_extent]

    #####################################################################
    ## first estimate of temporal features
    #####################################################################
    new_time, space, new_data = get_data(args.data_index,
                                         smoothing=np.ones((1, 4))/4**2,
                                         t0=args.t0, t1=args.t1)
    
    T1, T2 = get_time_residual(args,
                               new_time, space, new_data,
                               return_fit_directly=True)
    X0[2], X0[3] = 1e-3*T1, 1e-3*T2 # forcing temporal features to previous fitting, from ms to s
    
    #####################################################################
    ## fitting of spatial features
    #####################################################################
    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)

    def to_minimize(Xfit):
        """ X are the parameters """
        X = [Xfit[0], Xfit[1], X0[2], X0[3], Xfit[2]] # splitting between, fixed and to fit
        fn = run_sim(X, args)
        return get_space_residual(args,
                                  new_time, space, new_data,
                                  Nsmooth=args.Nsmooth,
                                  fn=fn)
    
    res = minimize(to_minimize, method='L-BFGS-B',
             x0=[X0[0], X0[1], X0[4]], # SET TEMPORAL FEATURES HERE
             bounds=[BOUNDS[0], BOUNDS[1], BOUNDS[4]],
             options={'maxiter':args.N})
    
    X0[0], X0[1], X0[4] = res.x # forcing spatial features to previous fitting

    #####################################################################
    ## fitting of temporal features
    #####################################################################
    new_time, space, new_data = get_data(args.data_index,
                                         smoothing=np.ones((1, 4))/4**2,
                                         t0=args.t0, t1=args.t1)
    
    def to_minimize(Xfit):
        """ X are the parameters """
        X = [X0[0], X0[1], Xfit[0], Xfit[1], X0[4]] # splitting between, fixed and to fit
        fn = run_sim(X, args)
        return get_time_residual(args,
                                 new_time, space, new_data,
                                 fn=fn)
    
    res = minimize(to_minimize, method='L-BFGS-B',
             x0=[X0[2], X0[3]], # SET TEMPORAL FEATURES HERE
             bounds=[BOUNDS[2], BOUNDS[3]],
             options={'maxiter':args.N})

    X0[2], X0[3] = res.x

    np.save('../ring_model/data/fitted_data_'+str(args.data_index)+'.npy', X0)

def plot_analysis(args):

    X0 = np.load('../ring_model/data/fitted_data_'+str(args.data_index)+'.npy')
    print(X0)

    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)
    fn = run_sim(X0, args, fn='../ring_model/data/model_data_'+str(args.data_index)+'.npy')
    res = get_residual(args,
                       new_time, space, new_data,
                       Nsmooth=args.Nsmooth,
                       fn=fn, with_plot=True)
    plt.show()


def full_analysis(args):

    DATA = get_dataset()
    VC, ECR, TAU2, TAU1, SEXT, DUR = [], [], [], [], [], []
    for i in range(len(DATA)):
        X = np.load('../ring_model/data/analyzed_scan_data_'+\
                    str(i)+'.npy')
        print(X.shape)
        if X.shape==(5,):
            for vec, VEC in zip(X, [VC, ECR, TAU2, TAU1, SEXT]):
                VEC.append(vec)
            DUR.append(DATA[i]['duration'])
        
    fig, AX = plt.subplots(1, 4, figsize=(6.,2.3))
    plt.subplots_adjust(bottom=.3, left=.25, wspace=3.)
    for ax, vec, label, ylim in zip(AX, [VC, ECR, np.array(ECR)/5., SEXT],
        ['$v_c$ (mm/s)', '$l_{exc}$ (mm)', '$l_{inh}$ (mm)', '$s_{ext}$ (mm)'],
                                    [[0,500], [0,6], [0,6], [0,4]]):
        ax.plot([0, 0], ylim, 'w.', ms=0.1)
        ax.bar([0], [np.array(vec).mean()], yerr=[np.array(vec).std()],
               color='lightgray', edgecolor='k', lw=3)
        set_plot(ax, ['left'], xticks=[], ylabel=label)

    fig2, AX = plt.subplots(1, 2, figsize=(3.2,2.3))
    plt.subplots_adjust(bottom=.3, left=.3, wspace=1.8)
    AX[0].plot(DUR, 1e3*np.array(TAU1), 'o')
    AX[0].plot(DUR, 1e3*np.array(TAU2), 'o')
    plt.show()


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
    parser.add_argument("-d", "--debug", help="with debugging", action="store_true")
    parser.add_argument("--full_plot", help="plot of full analysis", action="store_true")
    
    args = parser.parse_args()
    if args.fitting:
        run_fitting(args)
    elif args.plot:
        plot_analysis(args)
    elif args.full_plot:
        full_analysis(args)
    else:
        print('need arg')
