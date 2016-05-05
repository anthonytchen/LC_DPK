''' The module for QSPR model of partition coefficient
    between stratum corneum / corneocyte (and water) '''

import numpy as np

###########################################################
def qspr_lg10Ksc_nlh(paras, K_ow, K_sc, sig2=1, retSig2=False):
    ''' Function to calculate the negative log likelihood given <paras>, <K_ow>, <K_sc> and [sig2]'''

    import math
    import matplotlib.pyplot as plt

    n_dat = len(K_sc)

    sse = qspr_Ksc_sse(paras, K_ow, K_sc)
    likelihood = -0.5*n_dat*np.log(sig2) - 0.5*n_dat*np.log(2*math.pi) - 0.5*sse/sig2
    nlh = -likelihood
        
    if (retSig2):
        return sse/n_dat
    else:
        return nlh

###########################################################
def qspr_lg10Ksc_sse(paras, K_ow, K_sc, disp=False):
    ''' Function to calculate the SSE (sum of square error) given <paras>, <K_ow>, <K_sc>'''

    n_dat = len(K_sc)

    lg10Ksc_pred = qspr_lg10Ksc(paras, K_ow)
    err = np.log10(K_sc) - lg10Ksc_pred
    sse =  np.sum( np.square(err) )

    if (disp):
        plt.plot( np.log10(K_sc), lg10Ksc_pred, 'ro')
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls='r-', c=".c")
        plt.show()       

    return sse

###########################################################
def qspr_lg10Ksc(paras, K_ow):
    ''' Function to predict the log10 partition coefficient between stratum corneum (and water)'''
    
    K_sc_pred = qspr_Ksc(paras, K_ow)        
    return np.log10(K_sc_pred)


###########################################################
def qspr_Ksc(paras, K_ow):
    ''' Function to predict the partition coefficient between stratum corneum (and water)'''

    a = paras[0]
    b = paras[1]

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

    K_sc_pred = phi_pro*rho_pro/rho_wat* a*np.power(K_ow,b) + phi_lip*rho_lip/rho_wat* np.power(K_ow,0.69) + phi_wat
        
    return K_sc_pred

