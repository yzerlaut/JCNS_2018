import numpy as np
import matplotlib.pylab as plt
import itertools
# everything stored within a zip file
import zipfile, sys, os
sys.path.append("../experimental_data")
# from compare_to_model import *
sys.path.append("../../")
from graphs.my_graph import set_plot
from dataset import get_dataset
from compare_to_model import get_data, get_residual
# from compare_to_model import get_data, get_time_residual, get_space_residual, get_residual

FACTOR_FOR_MUVN_NORM = abs((-54.+58.)/58.) # 6% corresponding to a ~5mV wrt to rest level

def to_filename(vc, se, ecr, icr, t2, t1):
    return '../ring_model/data/scan_'+str(vc)+'_'+str(se)+'_'+str(ecr)+'_'+str(icr)+'_'+str(t2)+'_'+str(t1)+'.npy'

def create_grid_scan_bash_script(args):
    
    def cmd(vc, se, ecr, icr, t2, t1):
        fn = to_filename(vc, se, ecr, icr, t2, t1)
        return fn, 'python single_trial.py '+\
            ' --conduction_velocity_mm_s '+str(vc)+\
            ' --sX '+str(se)+\
            ' --exc_connect_extent '+str(ecr)+\
            ' --inh_connect_extent '+str(icr)+\
            ' --Tau2 '+str(t2)+' --Tau1 '+str(t1)+' -f '+fn+\
            ' --no_plot --X_extent 36 --X_discretization 30'

    VC = np.linspace(args.vc[0], args.vc[1], args.N)
    SE = np.linspace(args.stim_extent[0], args.stim_extent[1], args.N)
    ECR = np.linspace(args.Econn_radius[0], args.Econn_radius[1], args.N)
    ICR = np.linspace(args.Iconn_radius[0], args.Iconn_radius[1], args.N)
    TAU2 = np.linspace(args.Tau2[0], args.Tau2[1], args.N)
    TAU1 = np.linspace(args.Tau1[0], args.Tau1[1], args.N)

    f = open('bash_parameter_scan.sh', 'w')
    FILENAMES = []
    for vc, se, ecr, icr, t2, t1 in itertools.product(VC, SE, ECR, ICR, TAU2, TAU1):
        fn, c = cmd(vc, se, ecr, icr, t2, t1)
        if (t1<TAU1[-1]) or (t2<TAU2[-1]):
            c += ' & \n'
        else:
            c += ' \n'
        if (args.force==True) or (os.path.isfile(fn)==False):
            f.write(c)
        else:
            print('existing datafile')
        FILENAMES.append(fn)
    f.close()
    np.save('../ring_model/data/scan_data.npy', [VC, SE, ECR, ICR, TAU2, TAU1, np.array(FILENAMES)])


def zip_data(args):
    zf = zipfile.ZipFile(args.zip_filename, mode='w')
    # writing the parameters
    zf.write('../ring_model/data/scan_data.npy')
    VC, SE, ECR, ICR, TAU2, TAU1, FILENAMES = np.load('../ring_model/data/scan_data.npy')
    for fn in FILENAMES:
        zf.write(fn)
    zf.close()

def unzip_data(args):
    zf = zipfile.ZipFile(args.zip_filename, mode='r')
    # writing the parameters
    data = zf.read('../ring_model/data/scan_data.npy')
    with open('../ring_model/data/scan_data.npy', 'wb') as f: f.write(data)
    VC, SE, ECR, ICR, TAU2, TAU1, FILENAMES = np.load('../ring_model/data/scan_data.npy')
    for fn in FILENAMES:
        data = zf.read(fn)
        with open(fn, 'wb') as f: f.write(data)
    zf.close()
    
def analyze_scan(args):
    
    VC, SE, ECR, ICR, TAU2, TAU1, FILENAMES = np.load('../ring_model/data/scan_data.npy')
    print(TAU1[0], FILENAMES[0])
    
    Residuals, time_Residuals, spatial_Residuals = [], [], []
    vcFull, seFull, ecrFull, icrFull, t2Full, t1Full = [], [], [], [], [], []
    
    ## loading data for time residual
    new_time, space, new_data = get_data(args.data_index,
                                         smoothing=np.ones((1, 4))/4**2,
                                         t0=args.t0, t1=args.t1)
    print(TAU1, TAU2)
    for vc, se, ecr, icr, t2, t1 in itertools.product(VC, SE, ECR, ICR, TAU2, TAU1):
        fn = to_filename(vc, se, ecr, icr, t2, t1)
        try:
            res = get_residual(args,
                               new_time, space, new_data,
                               model_normalization_factor=FACTOR_FOR_MUVN_NORM,
                               fn=fn)
            Residuals.append(res)
            t2Full.append(t2)
            t1Full.append(t1)
            vcFull.append(vc)
            seFull.append(se)
            ecrFull.append(ecr)
            icrFull.append(icr)
        except (IOError, OSError):
            print('missing', fn)
        
    np.save('../ring_model/data/residuals_data_'+str(args.data_index)+'.npy',
            [np.array(Residuals), np.array(time_Residuals), np.array(spatial_Residuals),
             np.array(vcFull), np.array(seFull),
             np.array(ecrFull), np.array(icrFull),
             np.array(t2Full), np.array(t1Full)])

def plot_analysis(args):
    
    Residuals, time_Residuals, spatial_Residuals,\
        vcFull, seFull, ecrFull, icrFull,\
        t2Full, t1Full = np.load(\
            '../ring_model/data/residuals_data_'+str(args.data_index)+'.npy')

    i0 = np.argmin(Residuals)
    Residuals/=Residuals[i0] # normalizing
    
    fig, AX = plt.subplots(1, 6, figsize=(9,2.))
    plt.subplots_adjust(bottom=.3, left=.15)
    for ax, vec, label in zip(AX,
                              [vcFull, seFull, ecrFull, icrFull, t2Full, t1Full],\
                              ['$v_c (mm/s)$','$l_{stim}$ (mm)',
                               '$l_{exc}$ (mm)', '$l_{inh}$ (mm)',
                               '$\\tau_2$ (ms)', '$\\tau_1$ (ms)']):
        ax.plot(vec, Residuals, 'o')
        ax.plot([vec[i0]], [Residuals[i0]], 'ro')
        ax.set_yscale('log')
        if ax==AX[0]:
            set_plot(ax, xlabel=label, ylabel='Residual (norm.)',
                     yticks=[1, 2, 5, 10, 20], yticks_labels=['1', '2', '5', '10', '20'])
        else:
            set_plot(ax, xlabel=label, yticks=[1, 5, 10, 20], yticks_labels=[])

    _, _, _, _, _, _, FILENAMES = np.load('../ring_model/data/scan_data.npy')
    new_time, space, new_data = get_data(args.data_index,
                                         Nsmooth=args.Nsmooth,
                                         t0=args.t0, t1=args.t1)
    res = get_residual(args,
                       new_time, space, new_data,
                       Nsmooth=args.Nsmooth,
                       fn='../ring_model/'+FILENAMES[i0], with_plot=True)
    plt.show()

def get_minimum_params(args):
    Residuals, time_Residuals, spatial_Residuals,\
        vcFull, seFull, ecrFull, icrFull,\
        t2Full, t1Full = np.load(\
            '../ring_model/data/residuals_data_'+str(args.data_index)+'.npy')

    i0 = np.argmin(Residuals)
    Residuals/=Residuals[i0] # normalizing
    return vcFull[i0], seFull[i0], ecrFull[i0], icrFull[i0], t2Full[i0], t1Full[i0]

    
def full_analysis(args):

    DATA = get_dataset()
    for i in range(len(DATA)):
        print('analyzing cell ', i, ' [...]')
        args.data_index = i
        analyze_scan(args)

def full_plot(args):

    DATA = get_dataset()
    VC, SE, ECR, ICR, TAU2, TAU1, DUR, MONKEY = [], [], [], [], [], [], [], []
    for i in range(len(DATA)):
        args.data_index = i
        params = get_minimum_params(args)
        for vec, VEC in zip(params, [VC, SE, ECR, ICR, TAU2, TAU1]):
            VEC.append(vec)
        DUR.append(DATA[i]['duration'])
        MONKEY.append(DATA[i]['Monkey'])

    fig, AX = plt.subplots(1, 3, figsize=(6.,2.3))
    plt.subplots_adjust(bottom=.3, left=.25, wspace=4.)
    # vc
    print(MONKEY)
    i0, i1 = np.array(MONKEY)=='1', np.array(MONKEY)=='2'
    AX[0].bar([0], [np.array(VC)[i0].mean()], yerr=[np.array(VC)[i0].std()])
    AX[0].bar([1], [np.array(VC)[i1].mean()], yerr=[np.array(VC)[i1].std()])
    
    # for ax, vec, label, ylim in zip(AX, [VC, ECR, ICR, SE],
    #                     ['$v_c$ (mm/s)', '$l_{exc}$ (mm)', '$l_{inh}$ (mm)', '$l_{stim}$ (mm)'],
    #                                 [[0,600], [0,10], [0,10], [0,3]]):
    #     ax.plot([0, 0], ylim, 'w.', ms=0.1)
    #     ax.bar([0], [np.array(vec).mean()], yerr=[np.array(vec).std()],
    #            color='lightgray', edgecolor='k', lw=3)
    #     set_plot(ax, ['left'], xticks=[], ylabel=label)

    # fig2, AX = plt.subplots(1, 2, figsize=(3.2,2.3))
    # plt.subplots_adjust(bottom=.3, left=.3, wspace=1.8)
    # AX[0].plot(DUR, 1e3*np.array(TAU1), 'o')
    # AX[0].plot(DUR, 1e3*np.array(TAU2), 'o')
    # plt.show()
    
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--vc", nargs=2, type=float, default=[50., 600.])
    parser.add_argument("--stim_extent", nargs=2, type=float, default=[0.2, 2.])
    parser.add_argument("--Econn_radius", nargs=2, type=float, default=[1., 7.])
    parser.add_argument("--Iconn_radius", nargs=2, type=float, default=[1., 7.])
    parser.add_argument("--Tau1", nargs=2, type=float, default=[5e-3, 50e-3])
    parser.add_argument("--Tau2", nargs=2, type=float, default=[50e-3, 200e-3])
    parser.add_argument("--N", type=int, default=2)
    parser.add_argument("--zip_filename", '-f', type=str, default='../ring_model/data/data.zip')
    # data
    parser.add_argument("--data_index", '-df', type=int,
                        default=7)
    parser.add_argument("--t0", type=float, default=-100.)
    parser.add_argument("--t1", type=float, default=300.)
    parser.add_argument("--Nsmooth", help="for data plots", type=int, default=1)
    # script function
    parser.add_argument("-s", "--save", help="save fig", action="store_true")
    parser.add_argument("-a", "--analyze", help="analyze", action="store_true")
    parser.add_argument("-p", "--plot", help="plot analysis", action="store_true")
    parser.add_argument("-z", "--zip", help="zip datafiles", action="store_true")
    parser.add_argument("-uz", "--unzip", help="unzip datafiles", action="store_true")
    parser.add_argument("-d", "--debug", help="with debugging", action="store_true")
    parser.add_argument("--force", help="force simulation", action="store_true")
    parser.add_argument("--full", help="full analysis", action="store_true")
    parser.add_argument("--full_plot", help="plot of full analysis", action="store_true")
    
    args = parser.parse_args()
    if args.analyze:
        analyze_scan(args)
    elif args.plot:
        plot_analysis(args)
    elif args.zip:
        zip_data(args)
    elif args.unzip:
        unzip_data(args)
    elif args.full:
        full_analysis(args)
    elif args.full_plot:
        full_plot(args)
    else:
        create_grid_scan_bash_script(args)
