import numpy as np
import matplotlib.pylab as plt
import itertools

def create_grid_scan_bash_script(args):

    def cmd(vc, ecr, rIE, t2):
        fn = 'data/scan_'+str(vc)+'_'+str(ecr)+'_'+str(rIE)+'_'+str(t2)+'.npy'
        return fn, 'python single_trial.py '+\
            ' --conduction_velocity_mm_s '+str(vc)+\
            ' --exc_connect_extent '+str(ecr)+\
            ' --inh_connect_extent '+str(ecr*rIE)+\
            ' --Tau2 '+str(t2)+' -f '+fn+' --no_plot --X_extent 30 X_discretization 30 \n'

    VC = np.linspace(args.vc[0], args.vc[1], args.N)
    ECR = np.linspace(args.Econn_radius[0], args.Econn_radius[1], args.N)
    RIE = np.logspace(np.log(args.ratio_Iconn_Econn[0])/np.log(10),\
                      np.log(args.ratio_Iconn_Econn[1])/np.log(10), args.N)
    TAU2 = np.linspace(args.Tau2[0], args.Tau2[1], args.N)

    f = open('bash_parameter_scan.sh', 'w')
    FILENAMES = []
    for vc, ecr, rIE, t2 in itertools.product(VC, ECR, RIE, TAU2):
        fn, c = cmd(vc, ecr, rIE, t2)
        f.write(c)
        FILENAMES.append(fn)
    f.close()
    np.save('data/scan_data.npy', [VC, ECR, RIE, TAU2, np.array(FILENAMES)])

def analyze_scan(args):
    
    VC, ECR, RIE, TAU2, FILENAMES = np.load('data/scan_data.npy')
    print(FILENAMES)

    
    
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs a single trial with all options possible
            """,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--vc", nargs=2, type=float, default=[50., 600.])
    parser.add_argument("--stim_extent", nargs=2, type=float, default=[1., 10.])
    parser.add_argument("--Econn_radius", nargs=2, type=float, default=[1., 10.])
    parser.add_argument("--ratio_Iconn_Econn", nargs=2, type=float, default=[0.05, 1.5])
    parser.add_argument("--Tau1", nargs=2, type=float, default=[5., 100.])
    parser.add_argument("--Tau2", nargs=2, type=float, default=[50., 300.])
    parser.add_argument("--N", type=int, default=2)
    parser.add_argument("-a", "--analyze", help="analyze", action="store_true")
    
    args = parser.parse_args()
    if args.analyze:
        analyze_scan(args)
    else:
        create_grid_scan_bash_script(args)
