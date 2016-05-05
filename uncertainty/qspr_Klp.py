# The module for QSPR model of partition coefficient between lipid (in stratum corneum) (and water).
# As a convention, the parameters should be in natural log, while Kow, Plp and Klp are in log10
#
# We follow Int. J. Pharm 398 (2010) 114 to use Klp for volumetric partition coefficient, 
#     and Plp for mass partition coefficient.

import numpy as np

###########################################################
def qspr_lg10Plp_nlh(paras, lg10Kow, lg10Plp, sig2=1, retSig2=False):
    ''' Function to calculate the negative log likelihood given <paras>, <lg10Kow>, <lg10Plp> and [sig2]
    Note that the data given in lg10Klp is the MASS partition coefficient between lipid and water '''

    import math
    import matplotlib.pyplot as plt

    n_dat = lg10Kow.shape[0]

    sse = qspr_lg10Plp_sse(paras, lg10Kow, lg10Plp)
    likelihood = -0.5*n_dat*np.log(sig2) - 0.5*n_dat*np.log(2*math.pi) - 0.5*sse/sig2
    nlh = -likelihood
        
    if (retSig2):
        return sse/n_dat
    else:
        return nlh

###########################################################
def qspr_lg10Plp_sse(paras, lg10Kow, lg10Plp, disp=False):
    ''' Function to calculate the SSE (sum of square error) given <paras>, <lg10Kow>, <lg10Plp>
    Note that the data given in lg10Plp is the MASS partition coefficient between lipid and water
    '''

    n_dat = lg10Kow.shape[0]

    lg10Plp_pred = qspr_lg10Plp(paras, lg10Kow)
    err = lg10Plp - lg10Plp_pred
    sse =  np.sum( np.square(err) )

    if (disp):
        plt.plot( lg10Plp, lg10Plp_pred, 'ro')
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls='r-', c=".c")
        plt.show()       

    return sse


###########################################################
def qspr_Plp(paras, lg10Kow):
    ''' Function to predict the partition coefficient between lipid (and water)
    Note that the prediction is the MASS partition coefficient between lipid and water
    '''

    Plp_pred = 10**qspr_lg10Plp(paras, lg10Kow)
    return Plp_pred


###########################################################
def qspr_lg10Plp(paras, lg10Kow):
    ''' Function to predict the log10 partition coefficient between lipid (and water)
    Note that the prediction is the MASS partition coefficient between lipid and water
    '''
    
    paras = np.exp(paras)
    lg10Plp_pred = paras*lg10Kow
    return lg10Plp_pred

###########################################################
def qspr_Klp(paras, lg10Kow):
    ''' Function to predict the partition coefficient between lipid (and water)
    Note that the prediction is the VOLUMETRIC partition coefficient between lipid and water
    '''

    rho_lip = 0.90
    rho_wat = 1.00

    Klp_pred = rho_lip/rho_wat * qspr_Plp(paras, lg10Kow)
    return Klp_pred


###########################################################
def qspr_lg10Klp(paras, lg10Kow):
    ''' Function to predict the log10 partition coefficient between lipid (and water)
    Note that the prediction is the MASS partition coefficient between lipid and water
    '''
    
    lg10Klp_pred = np.log10( qspr_Klp(paras, lg10Kow) )
    return lg10Klp_pred
