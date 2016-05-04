import numpy as np
import math
import scipy.sparse as sp
import scipy.sparse.linalg as spln
import matplotlib.pyplot as plt

# Example:
# qspr_unctn.qspr(qspr_Kcc.qspr_Ksc_nlh, qspr_Kcc.qspr_lg10Ksc)

# todo list:
# 1. add qspr for Klip, Kcc, Dlip, Dcc
# 2. propagate the uncertainty in those parameters to prediction of permeability

###########################################################
def qspr(func_post, func_pred):
    ''' Function to build a QSPR model and to calculate parameter uncertainty
    Args:
    - func_post: the function to calculate the negative log posterior (or likelihood) of the QSPR model
    - func_pred: the function to predict the property
    '''

    from scipy.optimize import minimize

    dat = np.loadtxt("K_ow_sc.txt")
    K_ow = np.copy(dat[:,0])
    K_sc = np.copy(dat[:,1])
    n_dat = len(K_ow)

    paras0 = [4.2, 0.31]
    bnds = ((0, None), (0, None))
    disp = 0

    res = minimize(func_post, paras0, args=(K_ow,K_sc), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})
    sig2 = func_post(res.x, K_ow, K_sc, retSig2=True)
    print sig2

    H = calcHess(func_post, res.x, K_ow, K_sc, sig2)

    norm_approx_mean = res.x
    norm_approx_cov = np.linalg.inv(H)

    if (0):
        N = 200
        Samples = ImpSam(func_post, res.x, K_ow, K_sc, sig2, (norm_approx_mean, norm_approx_cov), N)

        K_sc_pred = np.zeros((N,n_dat))

        fig, ax = plt.subplots()

        for i in range(N):
            K_sc_pred[i,:] = func_pred(Samples[0][i,:], K_ow)# + np.random.normal(0, np.sqrt(sig2))
            plt.plot( np.log10(K_ow), np.log10(K_sc_pred[i,:])+np.sqrt(sig2)*1.96, 'y.');
            plt.plot( np.log10(K_ow), np.log10(K_sc_pred[i,:])-np.sqrt(sig2)*1.96, 'y.');

            K_sc_pred_mean = func_pred(res.x, K_ow)
            plt.plot( np.log10(K_ow), np.log10(K_sc_pred_mean), 'bx')
            plt.plot( np.log10(K_ow), np.log10(K_sc), 'ro')
            
        plt.show()       

    fig, ax = plt.subplots()
    lg10Ksc_pred = np.zeros((n_dat,2)) # 2 columns: mean & std
    for i in range(n_dat):
        lg10Ksc_pred[i,0], lg10Ksc_pred[i,1] = TransUnctnKF(func_pred, norm_approx_mean, norm_approx_cov, K_ow[i])
        lg10Ksc_pred[i,1] += sig2
        lg10Ksc_pred[i,1] = np.sqrt(lg10Ksc_pred[i,1])

        plt.plot( np.log10(K_ow[i]), lg10Ksc_pred[i,0]+lg10Ksc_pred[i,1]*1.96, 'y.');
        plt.plot( np.log10(K_ow[i]), lg10Ksc_pred[i,0]-lg10Ksc_pred[i,1]*1.96, 'y.');

    plt.plot( np.log10(K_ow), lg10Ksc_pred[:,0], 'bx')
    plt.plot( np.log10(K_ow), np.log10(K_sc), 'ro')
    plt.show()       

    return (norm_approx_mean, norm_approx_cov, Samples, lg10Ksc_pred)


###########################################################
def ImpSam(func_post, paras, X, Y, sig2, norm_proposal, N):
    ''' Function for importance sampling
    Args:
    -- norm_proposal: the normal proposal function, [0] is mean, [1] is cov
    '''

    from random import uniform

    samples = np.random.multivariate_normal(norm_proposal[0], norm_proposal[1], N)
    d = len(norm_proposal[0])
    
    proposal_prob = np.zeros((N,1))
    post_prob = np.zeros((N,1))
    weights = np.zeros((N,1))

    for i in range(N):
        post_prob[i] = - func_post(samples[i,:], X, Y, sig2) # the function 'func_post' returns negative log posterior
        proposal_prob[i] = lognormpdf(samples[i,:], norm_proposal[0], norm_proposal[1])
        weights[i] = math.exp( post_prob[i] - proposal_prob[i] )

    weights /= np.sum(weights)
    #print weights

    # resampling -- Algorithm 2 in IEEE Trans Sig Proc (2002) 50:174.

    cdf = np.zeros((N,1))
    cdf[0] = weights[0]
    for i in range(1,N):
        cdf[i] = cdf[i-1] + weights[i]
    #print cdf

    re_sampled = np.zeros((N,d))
    u0 = uniform(0, 1.0/N)   
    for j in range(N):
        uj = u0 + float(j)/N
        #print ('uj=', uj)
        for i in range(N):
            if (uj<=cdf[i]):
                break
        re_sampled[j] = samples[i]
        #print i

    #return (samples, weights)
    return (re_sampled, samples)

###########################################################
def TransUnctnKF(func, paras_mean, paras_cov, X):
    ''' Function to transform the uncertainty, represented as a normal distribution (paras_mean, paras_cov)
    through a deterministic (func) to calculate the normal approximation of the output uncertainty
    X is the additional inputs that are needed for func
    '''

    y_mean = func(paras_mean, X)

    n_paras = len(paras_mean)
    
    grad = np.zeros( (n_paras, 1) )
    
    delta_paras = np.fabs(paras_mean) * 1e-3
    delta_paras[ delta_paras<1e-8 ] = 1e-8 # to avoid too small values

    # finite difference approximation of the gradient
    for i in range(n_paras):
        p1 = np.copy(paras_mean)
        p1[i] += delta_paras[i]
        t1 = func(p1, X)
        grad[i] = (t1-y_mean) / delta_paras[i]

    y_cov = np.transpose(grad).dot( paras_cov.dot(grad) )

    return (y_mean, y_cov)



###########################################################
def calcHess(func_post, paras, X, Y, sig2=1):
    ''' Function to calculate the Hessian of negative
        log posterior w.r.t. model parameters '''

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


###########################################################
def lognormpdf(x,mu,S):
    """ Calculate gaussian probability density of x, when x ~ N(mu,sigma) """

    nx = len(S)
    norm_coeff = nx*math.log(2*math.pi)+np.linalg.slogdet(S)[1]

    err = x-mu
    if (sp.issparse(S)):
        numerator = spln.spsolve(S, err).T.dot(err)
    else:
        numerator = np.linalg.solve(S, err).T.dot(err)

    return -0.5*(norm_coeff+numerator)
