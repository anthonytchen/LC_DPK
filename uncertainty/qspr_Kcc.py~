''' The module for QSPR model of partition coefficient
    between stratum corneum / corneocyte (and water) '''

import numpy as np

###########################################################
def qspr_lg10Ksc_nlh(paras, lg10Kow, lg10Ksc, sig2=1, retSig2=False):
    ''' Function to calculate the negative log likelihood given <paras>, <lg10Kow>, <lg10Ksc> and [sig2]
    Note that the data given in lg10Ksc is the VOLUMETRIC partition coefficient between stratum corneum and water
    '''

    import math
    import matplotlib.pyplot as plt

    n_dat = lg10Kow.shape[0]

    sse = qspr_lg10Ksc_sse(paras, lg10Kow, lg10Ksc)
    likelihood = -0.5*n_dat*np.log(sig2) - 0.5*n_dat*np.log(2*math.pi) - 0.5*sse/sig2
    nlh = -likelihood
        
    if (retSig2):
        return sse/n_dat
    else:
        return nlh

###########################################################
def qspr_lg10Ksc_sse(paras, lg10Kow, lg10Ksc, disp=False):
    ''' Function to calculate the SSE (sum of square error) given <paras>, <lg10Kow>, <lg10Ksc>
    Note that the data given in lg10Ksc is the VOLUMETRIC partition coefficient between stratum corneum and water
    '''

    n_dat = lg10Kow.shape[0]

    lg10Ksc_pred = qspr_lg10Ksc(paras, lg10Kow)
    err = lg10Ksc - lg10Ksc_pred
    sse =  np.sum( np.square(err) )

    if (disp):
        plt.plot( lg10Ksc, lg10Ksc_pred, 'ro')
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls='r-', c=".c")
        plt.show()       

    return sse

###########################################################
def qspr_lg10Ksc(paras, lg10Kow):
    ''' Function to predict the log10 partition coefficient between stratum corneum (and water)
    Note that the prediction is the VOLUMETRIC partition coefficient between stratum corneum and water
    '''
    
    Ksc_pred = qspr_Ksc(paras, lg10Kow)        
    return np.log10(Ksc_pred)


###########################################################
def qspr_Ksc(paras, lg10Kow):
    ''' Function to predict the partition coefficient between stratum corneum (and water)
    Note that the prediction is the VOLUMETRIC partition coefficient between stratum corneum and water 
    '''

    a = np.exp(paras[0])
    b = np.exp(paras[1])
    K_ow = 10**lg10Kow

    w_pro = 0.77
    w_lip = 0.23
    w_wat = 2.99

    #w_pro = 0.45*0.875; 
    #w_lip = 0.45*0.125; 
    #w_wat = 0.55; 

    rho_pro = 1.37
    v_pro = w_pro/rho_pro
    rho_lip = 0.90
    v_lip = w_lip/rho_lip
    rho_wat = 1.00
    v_wat = w_wat/rho_wat

    v_total = v_pro + v_lip + v_wat
    phi_pro = v_pro / v_total
    phi_lip = v_lip / v_total
    phi_wat = v_wat / v_total

    Ksc_pred = phi_pro*rho_pro/rho_wat* a*np.power(K_ow,b) + phi_lip*rho_lip/rho_wat* np.power(K_ow,0.69) + phi_wat
        
    return Ksc_pred


def qspr_lg10Ksc(paras, lg10Kow, lg10Z=None):
    ''' Function to predict the VOLUMETRIC partition coefficient between stratum corneum (and water)
    Args:
        paras, lg10Kow are not used if lg10Z is present
        lg10Z contains [lg10Kcc, lg10Klp]
    '''

    if lg10Z is None:
        return qspr_lg10Ksc(paras, lg10Kow)

    lg10Kcc = lg10Z[0]
    lg10Klp = lg10Z[1]
    
    w_pro = 0.77
    w_lip = 0.23
    w_wat = 2.99

    rho_pro = 1.37
    v_pro = w_pro/rho_pro
    rho_lip = 0.90
    v_lip = w_lip/rho_lip
    rho_wat = 1.00
    v_wat = w_wat/rho_wat

    v_total = v_pro + v_lip + v_wat
    phi_pro = v_pro / v_total
    phi_lip = v_lip / v_total
    phi_wat = v_wat / v_total

    Ksc_pred = phi_pro*rho_pro/rho_wat* (10**lg10Kcc) + phi_lip*rho_lip/rho_wat* (10**lg10Klp) + phi_wat
        
    return np.log10(Ksc_pred)


###########################################################
def qspr_lg10Kcc(paras, lg10Kow):
    ''' Function to predict the log10 volumetric partition coefficient between corneocyte (and water)'''
    
    Kcc_pred = qspr_Kcc(paras, lg10Kow)
    return np.log10(Kcc_pred)


###########################################################
def qspr_Kcc(paras, lg10Kow):
    ''' Function to predict the volumetric partition coefficient between corneocyte (and water)'''

    a = np.exp(paras[0])
    b = np.exp(paras[1])
    K_ow = 10**lg10Kow

    rho_pro = 1.37
    rho_wat = 1.00

    Kcc_pred = rho_pro/rho_wat* a*np.power(K_ow,b)
        
    return Kcc_pred
