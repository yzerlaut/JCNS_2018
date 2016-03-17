import numpy as np


def build_up_differential_operator_first_order(TF1, TF2, params, T=5e-3):
    """
    simple first order system
    """
    def A0(V, exc_aff=0, inh_aff=0):
        return 1./T*(TF1(V[0]+exc_aff, V[1]+inh_aff)-V[0])
    
    def A1(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                .5*V[2]*diff2_fe_fe(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fe_fi(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fi_fe(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[4]*diff2_fi_fi(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                TF2(V[0]+exc_aff, V[1]+inh_aff)-V[0])
    
    def A2(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                1./Ne*TF(V[0]+exc_aff, V[1]+inh_aff)*(1./T-TF(V[0]+exc_aff, V[1]+inh_aff))+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[0])**2+\
                2.*V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                2.*V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (-2.*V[2]))
    
    def A3(V, exc_aff=0, inh_aff=0): # mu, nu = e,i, then lbd = e then i
        return 1./T*(\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[0])*(TF(V[0]+exc_aff, V[1]+inh_aff)-V[1])+\
                (-2.*V[3]))
    
    def A4(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                1./Ni*TF(V[0]+exc_aff, V[1]+inh_aff)*(1./T-TF(V[0]+exc_aff, V[1]+inh_aff))+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[1])**2+\
                2.*V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                2.*V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (-2.*V[4]))
    
    def Diff_OP(V, exc_aff=0, inh_aff=0):
        return np.array([A0(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A1(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A2(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A3(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A4(V, exc_aff=exc_aff, inh_aff=inh_aff)])
    return Diff_OP
    

def build_up_differential_operator_for_sym_exc_inh(TF1, TF2, params, T=5e-3):
    """
    Implements Equation (3.16) in El BOustani & Destexhe 2009
    in the case of a network of two populations:
    one excitatory and one inhibitory
    
    Each neuronal population has the same transfer function
    this 2 order formalism computes the impact of finite size effects
    T : is the bin for the Markovian formalism to apply

    the time dependent vector vector is V=[fe,fi, sfe, sfi, sfefi]
    the function returns Diff_OP
    and d(V)/dt = Diff_OP(V)
    """
    
    # size of populations
    Ne, Ni = params['Ntot']*(1-params['gei']), params['Ntot']*params['gei']

    # we have the transfer function, now we also get its derivatives
    # TF, diff_fe, diff_fi, diff2_fe_fe, diff2_fe_fi, diff2_fi_fi, values = \
    #                         get_derivatives_of_TF(params)
    
    def A0(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                .5*V[2]*diff2_fe_fe(TF1, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fe_fi(TF1, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fi_fe(TF1, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[4]*diff2_fi_fi(TF1, V[0]+exc_aff, V[1]+inh_aff)+\
                TF1(V[0]+exc_aff, V[1]+inh_aff)-V[0])
    
    def A1(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                .5*V[2]*diff2_fe_fe(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fe_fi(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[3]*diff2_fi_fe(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                .5*V[4]*diff2_fi_fi(TF2, V[0]+exc_aff, V[1]+inh_aff)+\
                TF2(V[0]+exc_aff, V[1]+inh_aff)-V[0])
    
    def A2(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                1./Ne*TF(V[0]+exc_aff, V[1]+inh_aff)*(1./T-TF(V[0]+exc_aff, V[1]+inh_aff))+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[0])**2+\
                2.*V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                2.*V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (-2.*V[2]))
    
    def A3(V, exc_aff=0, inh_aff=0): # mu, nu = e,i, then lbd = e then i
        return 1./T*(\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[0])*(TF(V[0]+exc_aff, V[1]+inh_aff)-V[1])+\
                (-2.*V[3]))
    
    def A4(V, exc_aff=0, inh_aff=0):
        return 1./T*(\
                1./Ni*TF(V[0]+exc_aff, V[1]+inh_aff)*(1./T-TF(V[0]+exc_aff, V[1]+inh_aff))+\
                (TF(V[0]+exc_aff, V[1]+inh_aff)-V[1])**2+\
                2.*V[2]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fe(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                V[3]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                2.*V[4]*diff_fi(TF, V[0]+exc_aff, V[1]+inh_aff)+\
                (-2.*V[4]))
    
    def Diff_OP(V, exc_aff=0, inh_aff=0):
        return np.array([A0(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A1(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A2(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A3(V, exc_aff=exc_aff, inh_aff=inh_aff),\
                         A4(V, exc_aff=exc_aff, inh_aff=inh_aff)])
    return Diff_OP


def diff_fe(TF, fe, fi, df=1e-4):
    return (TF(fe+df/2., fi)-TF(fe-df/2.,fi))/df

def diff_fi(TF, fe, fi, df=1e-4):
    return (TF(fe, fi+df/2.)-TF(fe, fi-df/2.))/df

def diff2_fe_fe(TF, fe, fi, df=1e-4):
    return (diff_fe(TF, fe+df/2., fi)-diff_fe(TF,fe-df/2.,fi))/df

def diff2_fi_fe(TF, fe, fi, df=1e-4):
    return (diff_fi(TF, fe+df/2., fi)-diff_fi(TF,fe-df/2.,fi))/df

def diff2_fe_fi(TF, fe, fi, df=1e-4):
    return (diff_fe(TF, fe, fi+df/2.)-diff_fe(TF,fe, fi-df/2.))/df

def diff2_fi_fi(TF, fe, fi, df=1e-4):
    return (diff_fi(TF, fe, fi+df/2.)-diff_fi(TF,fe, fi-df/2.))/df



