def qspr(func_post):
    ''' Function to build a QSPR model and to calculate parameter uncertainty'''

    import numpy as np
    from scipy.optimize import minimize

    dat = np.loadtxt("K_ow_sc.txt")
    K_ow = np.copy(dat[:,0])
    K_sc = np.copy(dat[:,1])

    paras0 = [4.2, 0.31]
    bnds = ((0, None), (0, None))
    disp = 0

    res = minimize(func_post, paras0, args=(K_ow,K_sc), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})
    sig2 = func_post(res.x, K_ow, K_sc, retSig2=True)

    H = calcHess(func_post, res.x, K_ow, K_sc, sig2)
    # print H
    # qspr_Kcc(paras0, dat, disp)
    return (res.x, np.linalg.inv(H))

def calcHess(func_post, paras, X, Y, sig2=1):
    ''' Function to calculate the Hessian of negative
        log posterior w.r.t. model parameters '''

    import numpy as np

    n_paras = len(paras)
    H = np.zeros( (n_paras, n_paras) )
    
    delta_paras = np.fabs(paras) * 1e-3
    delta_paras[ delta_paras<1e-8 ] = 1e-8 # to avoid too small values

    for i in range(n_paras):
        for j in range(n_paras):

            if (i>j):
                H[i,j] = H[j,i]

            else:
                p1 = np.copy(paras)                
                p1[i] += delta_paras[i]
                p1[j] += delta_paras[j]
                t1 = func_post(p1, X, Y, sig2)

                p2 = np.copy(paras)                
                p2[i] += delta_paras[i]
                p2[j] -= delta_paras[j]
                t2 = func_post(p2, X, Y, sig2)

                p3 = np.copy(paras)                
                p3[i] -= delta_paras[i]
                p3[j] += delta_paras[j]
                t3 = func_post(p3, X, Y, sig2)

                p4 = np.copy(paras)                
                p4[i] -= delta_paras[i]
                p4[j] -= delta_paras[j]
                t4 = func_post(p4, X, Y, sig2)

                H[i,j] = (t1-t2-t3+t4) / (4*delta_paras[i]*delta_paras[j])            

    return H


