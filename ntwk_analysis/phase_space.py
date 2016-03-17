from scipy import integrate
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
from transfer_functions.load_config import load_transfer_functions
from mean_field.master_equation import build_up_differential_operator_first_order

def draw_phase_space(args, T=5e-3):

    
    TF1, TF2 = load_transfer_functions(args.NRN1, args.NRN2, args.NTWK)

    ## now we compute the dynamical system
    ## X = [nu_e,nu_i]
    # T*d(nu_e)/dt = TF1(nu_e,nu_i) - n_e
    # T*d(nu_i)/dt = TF2(nu_e,nu_i) - n_i
    
    def dX_dt_scalar(X, t=0):
        return build_up_differential_operator_first_order(TF1, TF2, T=5e-3)(X, exc_aff=args.ext_drive)

    def dX_dt_mesh(IN, t=0):
        [x,y] = IN
        DFx = 0*x
        DFy = 0*y
        for inh in range(x.shape[0]):
            for exc in range(x.shape[1]):
                fe = x[inh][exc]
                fi = y[inh][exc]
                DFx[inh,exc] = (TF1(fe+args.ext_drive, fi)-fe)/T
                DFy[inh,exc] = (TF2(fe+args.ext_drive, fi)-fi)/T
        return DFx,DFy


    t = np.linspace(0,.04,1e5)              # time


    # values  = np.array([[3,9],[4,9.5],[7.5,7],[5.5,4],[4,4]])
    values  = np.array([[3,9],[4,12],[12,14],[13,8],[5.5,4]])
    vcolors = plt.cm.gist_heat(np.linspace(0, .8, len(values)))  # colors for each trajectory

    ff1 = plt.figure(figsize=(5,4))
    f1 = plt.subplot(111)
    plt.tight_layout()

    ff2 = plt.figure(figsize=(6,5))
    f2 = plt.subplot(111)
    plt.tight_layout()

    # X_f0 = np.array([4.,4.])
    # X_f1 = np.array([10.,10.]) # fixed points values

    #-------------------------------------------------------
    # plot trajectories
    for v, col in zip(range(len(values)), vcolors):
        X0 = values[v]      # starting point
        X = integrate.odeint(dX_dt_scalar, X0, t)         # we don't need infodict here
        l1 = f1.plot(1e3*t,X[:,0],'--',lw=2, color=col,label='excitation')
        l2 = f1.plot(1e3*t,X[:,1],'-',lw=2, color=col)
        f2.plot(X[:,0], X[:,1], lw=2, color=col, label='X0=(%.f, %.f)' % ( X0[0], X0[1]) )
    f1.legend(('exc. pop.','inh. pop.'), prop={'size':'x-small'})
    f1.set_xlabel("time (ms)")
    f1.set_ylabel("frequency (Hz)")

    ###### ============ VECTOR FIELD  ============= #######

    #-------------------------------------------------------
    # define a grid and compute direction at each point
    ymax = 15 ; ymin = 0
    xmax = 15 ; xmin = 0
    nb_points = 20

    x = np.linspace(xmin, xmax, nb_points)
    y = np.linspace(ymin, ymax, nb_points)

    X1,Y1 = np.meshgrid(x, y)        # create a grid
    DX1,DY1 = dX_dt_mesh([X1,Y1])   # compute growth rate on the gridt
    M = (np.hypot(DX1, DY1))         # Norm of the growth rate 
    M[ M==0 ] = 1.                  # Avoid zero division errors
    M /=1e3 # TO BE CHECKED -> Hz/ms
    DX1 /= M                         # Normalize each arrows
    DY1 /= M

    #-------------------------------------------------------
    # Drow direction fields, using matplotlib 's quiver function
    # I choose to plot normalized arrows and to use colors to give information on
    # the growth speed
    Q = plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=plt.cm.jet)
    cb = plt.colorbar(Q, cmap=plt.cm.jet,shrink=.5)
    cb.set_label('absolute speed (Hz/ms)')
    #cb.set_ticks([0,int(M.max()/2000)*1000,int(M.max()/1000)*1000])
    #cb.set_ticklabels([0,int(M.max()/2000),int(M.max()/1000)])
    plt.xlabel('exc. pop. freq. (Hz)')
    plt.ylabel('inh. pop. freq. (Hz)')
    # plt.legend(loc='best')
    # plt.grid()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    # now fixed point
    plt.plot([X[-1,0]], [X[-1,1]], 'ko', ms=10)
    plt.plot([X[-1,0], X[-1,0]], [0, X[-1,1]], 'k--', lw=4, alpha=.4)
    plt.plot([0, X[-1,0]], [X[-1,1], X[-1,1]], 'k--', lw=4, alpha=.4)
    print 'Fe=', X[-1,0]
    print 'Fi=', X[-1,1]

    plt.tight_layout()

    # ff2.savefig('figures/phase_space_with_fs.pdf', format='pdf')
    # ff1.savefig('figures/traject_with_fs.pdf', format='pdf')
    # ff2.savefig('figures/phase_space.pdf', format='pdf')
    # ff1.savefig('figures/traject.pdf', format='pdf')

    plt.show()


if __name__=='__main__':

    import argparse
    
    # First a nice documentation 
    parser=argparse.ArgumentParser(description=
     """ Runs two types of protocols on a given neuronal and network model
        1)  ==> Preliminary transfer function protocol ===
           to find the fixed point (with possibility to add external drive)
        2)  =====> Full transfer function protocol ==== 
           i.e. scanning the (fe,fi) space and getting the output frequency""",
              formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("NRN1",help="Choose a neuronal model")
    parser.add_argument("NRN2",help="Choose a neuronal model")
    parser.add_argument("NTWK",help="Choose a network model (synaptic and connectivity properties)")

    parser.add_argument("--ext_drive",type=float, default=0.,\
                        help="External drive (Hz)")
    
    parser.add_argument("--max_Fe",type=float, default=20.,\
                        help="Maximum excitatory frequency (default=20.)")
    parser.add_argument("--max_Fi",type=float, default=20.,\
                        help="Maximum inhibitory frequency (default=20.)")
    
    parser.add_argument("-s", "--save", help="save with the right name",
                         action="store_true")
    
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                         action="store_true")

    args = parser.parse_args()

    draw_phase_space(args)
    

