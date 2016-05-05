''' The module for QSPR model of diffusion coefficient
    in lipid (in stratum corneum)'''

import numpy as np
import math

###########################################################
def qspr_lnDlp_nlh(paras, MW, lnDlp, sig2=1, retSig2=False):
    ''' Function to calculate the negative log likelihood given <paras>, <MW>, <lnDlp> and [sig2]
    Args:
    - MW: molecular weight
    - lnDlp: measured diffusion coefficient in lipid (in natural log)
    '''

    import matplotlib.pyplot as plt

    n_dat = MW.shape[0]

    sse = qspr_lnDlp_sse(paras, MW, lnDlp)
    likelihood = -0.5*n_dat*np.log(sig2) - 0.5*n_dat*np.log(2*math.pi) - 0.5*sse/sig2
    nlh = -likelihood
        
    if (retSig2):
        return sse/n_dat
    else:
        return nlh

###########################################################
def qspr_lnDlp_sse(paras, MW, lnDlp, disp=False):
    ''' Function to calculate the SSE (sum of square error) given <paras>, <MW>, <lnDlp>'''

    n_dat = MW.shape[0]

    lnDlp_pred = qspr_lnDlp(paras, MW)
    err = lnDlp - lnDlp_pred
    sse =  np.sum( np.square(err) )

    if (disp):
        plt.plot( lnDlp, lnDlp_pred, 'ro')
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls='r-', c=".c")
        plt.show()       

    return sse

###########################################################
def qspr_Dlp(paras, MW):
    ''' Function to predict the diffusion coefficient in lipid
    D = a * exp(-b * r^2) where r can be calculated from MW  '''

    D_lp_pred = np.exp( qspr_lnDlp(paras, MW) )
    return D_lp_pred

###########################################################
def qspr_lnDlp(paras, MW):
    ''' Function to predict the ln diffusion coefficient in lipid '''

    lna = paras[0]
    b = np.exp(paras[1])
    r = np.power( 0.91*MW/4*3/math.pi, 0.333333333333333333333 )

    lnDlp_pred = lna - b*r*r
    return lnDlp_pred
