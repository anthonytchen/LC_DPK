''' The module that contains computational routines for hybrid model
identification and uncertainty quantification
'''

import numpy as np
import math
import scipy.sparse as sp
import scipy.sparse.linalg as spln
#import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping, brute

# Example:

# todo list:


###########################################################
def Estep(func_top, func_low, theta, X, Y, sig2_y, sig2_z, N=100):
    ''' The E-step of the Monte Carlo EM algorithm
    Args:
    - func_top: top level function;  y = f(theta,x,z) + epsilon
    - func_low: low level function;  z = h(theta,x)   + zeta
    - theta: model parameters
    - X: data for input variables
    - Y: data for top level output variables
    - sig2_y: variance of top level function noise (epsilon)
    - sig2_z: variance of low level function noise (zeta)
    - N: number of MC samples
    Rtns:
    '''

    n_dat = X.shape[0]
    d_y = Y.shape[1]
    d_z = 2 #func_low(theta, np.array(X[0,:])).shape[0] # dimension of z

    Ymc = np.zeros((N, d_y))
    Z = np.zeros((N, d_z, n_dat)) # MC samples
    W = np.zeros((N, n_dat)) # weights

    sig_y = np.sqrt(sig2_y)
    sig_z = np.sqrt(sig2_z)

    for n in range(n_dat):

        z_n = func_low(theta, np.array(X[n,:]))        
        Z[:,:,n] = np.tile(z_n,(N,1))
        np.tile(z_n,(N,1))

        # MC sampling
        for d in range(d_z):
            rd = np.random.normal(0, sig_z[d], N)
            Z[:,d,n] += rd

        s = .0;
        for i in range(N):
            Ymc[i,:] = func_top( theta, np.array(X[n,:]), Z[i,:,n] )
            W[i,n] = np.exp( lognormpdf( Ymc[i,:], Y[n,:], np.diag(sig2_y) ) )
            s += W[i,n]
        W[:,n] /= s # normalising the weights

    return (Z, W)


###########################################################
def Mstep_main(func_top, func_low, theta0, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z):
    
    res = minimize(Mstep_theta_obj, theta0, args=(func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z), 
                   method='CG', options={'disp': True, 'maxiter': 100})

    theta1 = res.x      
    #theta1 = theta0
    #var = Mstep_var(theta1, func_top, func_low, X, Y, Z, W)                   
    var = 0
    
    return (theta1, var)

def Mstep_theta_obj(theta, func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z):
    ''' The objective function (negative log-likelihood) for optimising theta
    Terms that are not dependent on theta (thus constant as far as optimisation is concerned)
    are not calculated.
    '''

    n_dat_Xy = Xy.shape[0]
    Zmc = Zsamples[0]
    Wmc = Zsamples[1]    
    n_mc = Zmc.shape[0]
    
    beta = 1.0/sig2_z
    alpha = 1.0/sig2_y

    neg_lnlik = .0
    
    # For high-level X-Y data
    for n in range(n_dat_Xy):
        
        z_func = func_low(theta, np.array(Xy[n,:]))

        for i in range(n_mc):
            err = Zmc[i,:,n] - z_func
            neg_lnlik += Wmc[i,n]* 0.5* np.sum( beta* (np.array(err)**2) )
            y_func = func_top(theta, Xy[n,:], Zmc[i,:,n])
            err = Y[n,:] - y_func
            neg_lnlik += Wmc[i,n]* 0.5*(alpha*err*err)

    # For low-level X-Z data
    d_z = len(Xz)
    
    for i in range (d_z):
        Xtmp = Xz[i]
        Ztmp = Z[i]
        n_dat_Xz = Xtmp.shape[0]   
        for n in range(n_dat_Xz):       
            z_func = func_low(theta, np.array(Xtmp[n,:]))
            err = Ztmp[n] - z_func[:,i]
            neg_lnlik += 0.5* np.sum( beta[i]* (np.array(err)**2) )
            
            
    return neg_lnlik


def Mstep_theta_grad(func_top, func_low, theta, X, Y, sig2_y, sig2_z, Z, W):
    ''' The gradient of the objective function with respect to theta
    '''
    n_dat = X.shape[0]
    d_y = Y.shape[1]
    n_mc = Z.shape[0]
    d_z = Z.shape[1]

    beta = 1.0/sig2_z
    alpha = 1.0/sig2_y

    grad = np.zeros(theta.shape)
    for n in range(n_dat):
        
        z_func = func_low(theta, X[n,:])
        y_func = func_top(theta, X[n,:], Z[:,:,n])
        gd_h = calc_grad(func_low, theta, X[n,:])

        for i in range(n_mc):
            err = Z[i,:,n] - z_func
            gd1 = beta*err* gd_h
            err = Y[n,:] - y_func
            gd2 = alpha*err* calc_grad(func_top, theta, X[n,:], Z[i,:,n])
            grad += W[i,n]* (gd1+gd2)

    return grad

def calc_grad(func, theta, X, Z=None):
    '''Calculate the gradient of <func> w.r.t parameters <theta> with given <X> and <Z>
    using finite difference
    '''
    n_theta = theta.shape[0]
    theta1 = np.array(theta.shape)
    grad = np.array(theta.shape)

    for i in range(n_theta):
        delta = theta[i]*1e-4
        theta1 = theta
        if np.abs(delta) < 1e-5:
            delta = 1e-5
        theta1[i] += delta

        if Z is not None:
            f = func(theta, X, Z)
            f1 = func(theta1, X, Z)
        else:
            f = func(theta, X)
            f1 = func(theta1, X)

        grad[i] = (f1-f) / delta
    
# to be updated, and need to be optimised iteratively as well
def Mstep_var(theta, func_top, func_low, X, Y, Z, W):
    ''' The function to calculate the optimal values of the variance terms
    '''

    n_dat = X.shape[0]
    d_y = Y.shape[1]
    n_mc = Z.shape[0]
    d_z = Z.shape[1]

    beta = np.zeros([1,d_z])
    alpha = np.zeros([1,d_y])

    for n in range(n_dat):
        
        z_func = func_low(theta, np.array(X[n,:]))        

        for i in range(n_mc):
            beta += W[i,n]* np.array( Z[i,:,n] - z_func )**2
            y_func = func_top(theta, X[n,:], Z[i,:,n])
            alpha += W[i,n]* np.array( Y[n,:] - y_func )**2

    #beta = N/beta;
    #alpha = N/alpha;      

    return (alpha, beta)

def testFunc_top(theta, X, Z):
    ''' Function to predict the partition coefficient between stratum corneum (and water)
    Note that the prediction is the VOLUMETRIC partition coefficient between stratum corneum and water 
    Here used as a test function of the top-level model
    Args:
      theta -- model parameters
      X -- lg10Kow, np.mat, dim: [n_dat, 1]
      Z -- lg10Kcc & lg10Klp, np.mat, dim: [n_dat, 2]
    Rtns:
      Y -- Ksc_pred, predicted coefficient between stratum corneum (and water)
    '''

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

    Kcc = np.power(10, Z[:,0])
    Klp = np.power(10, Z[:,1])

    Ksc_pred = phi_pro*rho_pro/rho_wat* Kcc + phi_lip*rho_lip/rho_wat* Klp + phi_wat
    Y =  Ksc_pred
        
    return Y


def testFunc_low(theta, X):
    ''' Function to predict the volumetric partition coefficient between corneocyte (and water)
    and that between lipid (and water)
    Here used as a test function of the low-level model
    Args:
      theta -- model parameters
      X -- lg10Kow, dim: [n_dat, 2]
    Rtns:
      Z -- dim: [n_dat, 2]
    '''
    a = np.exp(paras[0])
    b = np.exp(paras[1])
    c = np.exp(paras[2])

    n_dat = X.shape[0]
    n_x = X.shape[1]

    Z = np.zeros(X.shape)
    Z = np.mat(Z)

    rho_pro = 1.37
    rho_lip = 0.90
    rho_wat = 1.00

    for i in range(n_dat):
        if ( X[i, 0] != None ):
            Z[i, 0] = rho_pro/rho_wat* a*np.power(K_ow,b) # stratum corneum
        elif ( X[i, 1] != None ):
            Z[i, 1] = rho_lip/rho_wat* np.power(K_ow,c) # lipid
        else:
            print "error"
        
    return Z


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
