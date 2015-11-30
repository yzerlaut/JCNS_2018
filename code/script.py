""" documentation for this file is written at the end
 (use of argparse) """

import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/home/yann/work/python_library/') #

# download (i mean fork it and contribute !) my customizations for matplotlib graphics at :
# https://bitbucket.org/yzerlaut/python_library/src/7a803695fd16e6abbaa6d37c986363a92c486674/my_graph.py

try: # if you did download it
    from my_graph import set_plot
except ImportError: # if you didn't, need to define the function anyway, it will be use throughout the code
    def set_plot(ax, xlabel='', ylabel='', xticks=None, yticks=None, xticks_labels=None, yticks_labels=None):
        ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)

def generate_white_noise(mean, std, nsample=100):
    return mean+std*np.random.randn(nsample)

# in case not used as a modulus
if __name__=='__main__':
    import argparse
    # First a nice documentation 
    parser=argparse.ArgumentParser(description=
     """ 
     Generating random sample of a given distributions and
     comparing it with its theoretical value
     """
    ,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--DISTRIBUTION",\
        help="""
        Choose the random distribution type :
        WN : White Noise 
        OU : Ornstein-Uhlenbeck
        """, default='WN')
    parser.add_argument("--mean",help="mean of the random values",\
                        type=float, default=5.)
    parser.add_argument("--std",help="std of the random values",\
                        type=float, default=10.)
    parser.add_argument("--n",help="number of random events",\
                        type=int, default=2000)
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("-s", "--save", help="save the figures",
                        action="store_true")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(5,3))
    plt.subplots_adjust(left=.2,bottom=.2)
    z = np.log(np.abs(generate_white_noise(args.mean, args.std, nsample=args.n)))
    plt.hist(z, color='lightgray', edgecolor='k', lw=3)
    set_plot(ax, xlabel='xlabel (unit)', ylabel='ylabel (#)')
    fig.savefig('../figures/log_WN_hist.svg', format='svg')

    fig, ax = plt.subplots(figsize=(17,3))
    plt.subplots_adjust(bottom=.2)
    for i in range(10):
        z = generate_white_noise(\
            args.mean+10*i, args.std, nsample=args.n)
        plt.hist(z, lw=3)
    set_plot(ax, xlabel='xlabel (unit)', ylabel='ylabel (#)')
    fig.savefig('../figures/many_WN_hist.svg', format='svg')
        
    plt.show()
