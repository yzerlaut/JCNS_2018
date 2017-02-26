import numpy as np
import sys, os
sys.path.append(os.path.expanduser('~')+os.path.sep+'work')
from graphs.my_graph import set_plot

def zca_whiten(X, EPS = 1e-5):
    """
    https://gist.github.com/iborko/5d9c2c16004ce8b926ea

    X: numpy 2d array
        input data, rows are data points, columns are features

    Returns: ZCA whitened 2d array
    """
    assert(X.ndim == 2)
    
    mu = np.mean(X)
    
    #   covariance matrix
    cov = np.dot(X.T, X)
    #   d = (lambda1, lambda2, ..., lambdaN)
    d, E = np.linalg.eigh(cov)
    #   D = diag(d) ^ (-1/2)
    D = np.diag(1. / np.sqrt(d + EPS))
    #   W_zca = E * D * E.T
    whMat = np.dot(np.dot(E, D), E.T)

    X_white = np.dot(X, whMat)
    return X_white, mu, np.linalg.pinv(whMat), whMat


def dowhitefunc(X,whMat):
    
    # % Apres avoir utilise la fonction whiten.m -> recuperer la matrice whMat
    # % Pour l'appliquer sur les donnees
    
    Xwh2 = (X-np.mean(X)*np.diag(X.shape))*whMat;
    return Xwh2 


def decode_stim(t, S1, S2, AM, stim_length=100e-3, frame_norm=True, non_zero=1e-10):

    # extracting temporal quantities
    dt = t[1]-t[0] # time step
    istim = int(stim_length/dt) # duration in bins of response to determine cort. repres.
    # flatten over space S1 and S2 to extract onset times
    fS1, fS2 = S1.mean(axis=1), S2.mean(axis=1)
    iT1, iT2 = np.argmax(fS1), np.argmax(fS1)
    
    # mean cortical representations
    meanS1 = S1[iT1:iT1+istim,:].mean(axis=0)
    meanS2 = S2[iT2:iT2+istim,:].mean(axis=0)
    meanS12 = meanS1+meanS2
    # then normalization
    meanS1 /= non_zero+np.linalg.norm(meanS1)
    meanS2 /= non_zero+np.linalg.norm(meanS2)
    meanS12 /= non_zero+np.linalg.norm(meanS12)
    Blank = np.zeros(len(meanS1)) # blank is zero 

    LP = S1+S2 # linear prediction

    dAM, dLP = {}, {} # dictionaries that we will fill with the analysis values
    for key in ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL',
                'P_S1', 'P_S2', 'P_S12', 'P_BL']:
        dAM[key], dLP[key] = np.zeros(len(t)), np.zeros(len(t)) # initialized to zeros vectos
        
    # now loop over frames to get all the deviations to get an average deviation for the proba calc
    # for i in [500, 1100]:
    for i in range(len(t)):
        # --> Apparent Motion
        nAM = AM[i,:]/(non_zero+np.linalg.norm(AM[i,:]))
        # --> Linear Prediction
        nLP = LP[i,:]/(non_zero+np.linalg.norm(LP[i,:]))
        # computing deviations
        for key, stim in zip(\
                             ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL'],
                             [meanS1, meanS2, meanS12, Blank]):
            dAM[key][i] = np.std(nAM-stim)
            dLP[key][i] = np.std(nLP-stim)

    # now loop over frames to get all the deviations to get an average deviation for the proba calc
    for i in range(len(t)):

        # --> Apparent Motion
        nAM = AM[i,:]/(non_zero+np.linalg.norm(AM[i,:]))
        # --> Linear Prediction
        nLP = LP[i,:]/(non_zero+np.linalg.norm(LP[i,:]))
        # computing deviations
        Ptot_AM, Ptot_LP = 0, 0 # total proba
        for key1, key2, stim in zip(\
                             ['P_S1', 'P_S2', 'P_S12', 'P_BL'],
                             ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL'],
                             [meanS1, meanS2, meanS12, Blank]):
            dAM[key1][i] = np.exp(-0.5*np.sum((nAM-stim)**2/(1e-9+dAM[key2].std())))+1e-9
            Ptot_AM += dAM[key1][i]
            dLP[key1][i] = np.exp(-0.5*np.sum((nLP-stim)**2/(1e-9+dLP[key2].std())))+1e-9
            Ptot_LP += dLP[key1][i]

        for key in ['P_S1', 'P_S2', 'P_S12', 'P_BL']:
            dAM[key][i] /= Ptot_AM
            dLP[key][i] /= Ptot_LP

    fig, AX = plt.subplots(1, 2, figsize=(7,2.5))
    AX[0].plot(nAM)
    AX[0].plot(meanS1)
    AX[1].plot(nLP)
    AX[1].plot(meanS12)
    AX[1].plot(meanS2)
    plt.show()
    
    # figure
    fig, AX = plt.subplots(1, 2, figsize=(7,2.5))
    plt.subplots_adjust(bottom=.2, wspace=.3, hspace=.3)
    for key, label in zip(['P_S1', 'P_S2', 'P_S12', 'P_BL'],
                          ['P(1)', 'P(2)', 'P(1+2)', 'P(Blank)']):
        for ax, D in zip(AX[:2], [dAM, dLP]):
            ax.plot(1e3*t, D[key], label=label)
            
    for ax, title in zip(AX, ['App. Motion', 'Linear Pred.']):
        ax.legend(prop={'size':'x-small'})
        ax.set_title(title)
        set_plot(ax, yticks=[0,.5,1], ylabel='Proba', xlabel='time (ms)')

    desambuig = dAM['P_S2']-dLP['P_S2']
    return fig, desambuig

### revwriting without normalisation

def decode_stim(t, S1, S2, AM, stim_length=100e-3, frame_norm=True, sharpness_factor=1e-1):

    # extracting temporal quantities
    dt = t[1]-t[0] # time step
    istim = int(stim_length/dt) # duration in bins of response to determine cort. repres.
    # flatten over space S1 and S2 to extract onset times
    fS1, fS2 = S1.mean(axis=1), S2.mean(axis=1)
    iT1, iT2 = np.argmax(fS1), np.argmax(fS1)
    
    # mean cortical representations
    meanS1 = S1[iT1:iT1+istim,:].mean(axis=0)
    meanS2 = S2[iT2:iT2+istim,:].mean(axis=0)
    meanS12 = meanS1+meanS2
    Blank = np.zeros(len(meanS1)) # blank is zero 

    LP = S1+S2 # linear prediction

    dAM, dLP = {}, {} # dictionaries that we will fill with the analysis values
    for key in ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL',
                'P_S1', 'P_S2', 'P_S12', 'P_BL']:
        dAM[key], dLP[key] = np.zeros(len(t)), np.zeros(len(t)) # initialized to zeros vectos
        
    # now loop over frames to get all the deviations to get an average deviation for the proba calc
    # for i in [500, 1100]:
    for i in range(len(t)):
        # --> Apparent Motion
        nAM = AM[i,:]
        # --> Linear Prediction
        nLP = LP[i,:]
        # computing deviations
        for key, stim in zip(\
                             ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL'],
                             [meanS1, meanS2, meanS12, Blank]):
            dAM[key][i] = np.sum((nAM-stim)**2)
            dLP[key][i] = np.sum((nLP-stim)**2)

    denom = 0
    for key in ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL']:
        denom += sharpness_factor*dAM[key].std()/4.

        
    # now loop over frames to get all the deviations to get an average deviation for the proba calc
    for i in range(len(t)):

        # --> Apparent Motion
        nAM = AM[i,:]
        # --> Linear Prediction
        nLP = LP[i,:]
        # computing deviations
        Ptot_AM, Ptot_LP = 0, 0 # total proba
        for key1, key2, stim in zip(\
                             ['P_S1', 'P_S2', 'P_S12', 'P_BL'],
                             ['Dev_vs_S1', 'Dev_vs_S2', 'Dev_vs_S12', 'Dev_vs_BL'],
                             [meanS1, meanS2, meanS12, Blank]):
            dAM[key1][i] = np.exp(-np.sum((nAM-stim)**2)/denom)
            Ptot_AM += dAM[key1][i]
            dLP[key1][i] = np.exp(-np.sum((nLP-stim)**2)/denom)
            Ptot_LP += dLP[key1][i]
            # dAM[key1][i] = np.exp(-0.5*np.sum((nAM-stim)**2/(1e-9+dAM[key2].mean())))+1e-9
            # Ptot_AM += dAM[key1][i]
            # dLP[key1][i] = np.exp(-0.5*np.sum((nLP-stim)**2/(1e-9+dLP[key2].mean())))+1e-9
            # Ptot_LP += dLP[key1][i]

        for key in ['P_S1', 'P_S2', 'P_S12', 'P_BL']:
            dAM[key][i] /= Ptot_AM
            dLP[key][i] /= Ptot_LP

    # figure
    fig, AX = plt.subplots(1, 2, figsize=(7,2.5))
    plt.subplots_adjust(bottom=.2, wspace=.3, hspace=.3)
    for key, label in zip(['P_S1', 'P_S2', 'P_S12', 'P_BL'],
                          ['P(1)', 'P(2)', 'P(1+2)', 'P(Blank)']):
        for ax, D in zip(AX[:2], [dAM, dLP]):
            ax.plot(1e3*t, D[key], label=label)
            
    for ax, title in zip(AX, ['App. Motion', 'Linear Pred.']):
        ax.legend(prop={'size':'x-small'})
        ax.set_title(title)
        set_plot(ax, yticks=[0,.5,1], ylabel='Proba', xlabel='time (ms)')

    desambuig = np.abs(dAM['P_S2']-dLP['P_S2'])
    return fig, desambuig


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=
            """
            runs the apparent motion protocol
            """,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-s", "--SAVE",help="save the figures as SVG", action="store_true")
    parser.add_argument("-f", "--file",help="filename for saving",\
                        default='../apparent_motion/data/hd.npy')
    
    args = parser.parse_args()
    # we perform one experiment with the default and we store the figs

    import matplotlib.pylab as plt
    
    args2, t, X, Fe_aff1, Fe1, Fi1, muVn1,\
      Fe_aff2, Fe2, Fi2, muVn2, Fe_aff3, Fe3, Fi3, muVn3 = np.load(args.file)

    fig1, desambuigVSD = decode_stim(t, muVn1, muVn2, muVn3, stim_length=100e-3)
    fig2, desambuigFR = decode_stim(t, Fe1-Fe1[0,:].mean(), Fe2-Fe2[0,:].mean(),
                                    Fe3-Fe3[0,:].mean(), stim_length=100e-3)
    
    # fig1.savefig('data/VSD.svg')
    # fig2.savefig('data/FR.svg')

    fig, ax = plt.subplots(1, figsize=(3.5,2.5))
    plt.subplots_adjust(bottom=.3, left=.3)
    ax.plot(1e3*t, desambuigFR, label='Firing Rate')
    ax.plot(1e3*t, desambuigVSD, label='Vm')
    ax.legend(prop={'size':'x-small'})
    set_plot(ax, ylabel='Desambuiguation \n $P^{AM}(2) - P^{LP}(2)$', xlabel='time (ms)')
    
    plt.show()





