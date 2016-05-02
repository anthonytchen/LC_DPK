def qspr():
    ''' Function to build a QSPR model and to calculate parameter uncertainty'''

    import numpy as np
    from scipy.optimize import minimize

    dat = np.loadtxt("K_ow_sc.txt")
    K_ow = dat[:,0]
    K_sc = dat[:,1]

    paras0 = [4.2, 0.31]
    bnds = ((0, None), (0, None))
    disp = 0

    res = minimize(qspr_Kcc, paras0, args=(dat,disp), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})
    print res.x

    disp = 1
    #qspr_Kcc(res.x, dat, disp)
    qspr_Kcc(paras0, dat, disp)


def qspr_Kcc(paras, dat, disp):
    ''' Function to build a QSPR model for the partition coefficient in corneocyte'''

    import numpy as np
    import matplotlib.pyplot as plt

    a = paras[0]
    b = paras[1]
    K_ow = dat[:,0] # partition coefficient octanol -- water
    K_sc = dat[:,1] # partition coefficient stratum corneum -- water

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

    K_sc_pred = phi_pro*rho_pro/rho_wat* np.power(a*K_ow,b) + phi_lip*rho_lip/rho_wat* np.power(K_ow,0.69) + phi_wat;
    err = np.log10(K_sc) - np.log10(K_sc_pred)

    rase = np.sqrt( np.mean( np.square(err) ) )

    if (disp):
        plt.plot( np.log10(K_sc), np.log10(K_sc_pred), 'ro');
        plt.show()
        #refline(1,0);    
        
    return rase

