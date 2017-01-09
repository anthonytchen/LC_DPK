''' The module that contains computational routines for hybrid model
identification and uncertainty quantification
'''

import warnings
import numpy as np
import math
import scipy.sparse as sp
import scipy.sparse.linalg as spln
#import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping, brute

# global variable
Nfeval = 1

# Example:

# todo list:
# 1. add a function to allow plug-in type model identification
# 2. simulate data to compare plug-in and full bayesian approach

###########################################################
def EMmain(func_top, func_low, Xy, Y, Xz, Z, theta0, sig2_y0, sig2_z0, Nmc=10, Niter=10, bnds=None):
    ''' The main function to run EM algorithm for parameter estimation
    Args:
    Rtns:
    '''

    theta = np.empty_like(theta0)
    sig2_y = np.empty_like(sig2_y0)
    sig2_z = np.empty_like(sig2_z0)
    np.copyto(theta, theta0)
    np.copyto(sig2_y, sig2_y0)
    np.copyto(sig2_z, sig2_z0)
    
    for i in range(Niter):
        Zsamples = Estep(func_top, func_low, theta, Xy, Y, sig2_y, sig2_z, Nmc)
        feval = Mstep_theta_obj(theta, func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z)
        print 'EM Iter {0:4d}, f = {1:.6f}'.format(i, float(feval[0]))
        
        Mrlt = Mstep_main(func_top, func_low, theta, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z, bnds)
        np.copyto(theta, Mrlt[0])
        np.copyto(sig2_y, Mrlt[1])
        np.copyto(sig2_z, Mrlt[2])

    # calculate the hessian to determine the variability of the parameter estimate
    H = calcHessVargs(Mstep_theta_obj_wrapper, theta, func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z)
    V = np.linalg.inv(H)
    
    return (theta, sig2_y, sig2_z, V, Zsamples)


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

        z_n = func_low(theta, X[n,:])        
        Z[:,:,n] = np.tile(z_n,(N,1))
        np.tile(z_n,(N,1))

        # MC sampling
        for d in range(d_z):
            rd = np.random.normal(0, sig_z[d], N)
            Z[:,d,n] += rd

        s = .0;
        for i in range(N):
            Ymc[i,:] = func_top( theta, np.array(X[n,:]), Z[i,:,n] )
            aa = np.exp( lognormpdf( Ymc[i,:], Y[n,:], np.diag(sig2_y) ) )
            if np.isfinite(aa):
                W[i,n] = aa
            else:
                W[i,n] = 0
            s += W[i,n]
        W[:,n] /= s # normalising the weights

    return (Z, W)


###########################################################
def Mstep_main(func_top, func_low, theta0, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z, bnds=None):
    ''' Function to run the M-step
    '''
    
    global Nfeval
    Nfeval = 1

    if bnds == None:
        res = minimize(Mstep_theta_obj, theta0, args=(func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z), 
                       method='BFGS', callback=callbackF, options={'disp': True, 'maxiter': 100})
    else:
        res = minimize(Mstep_theta_obj, theta0, args=(func_top, func_low, Xy, Y, Zsamples, Xz, Z, sig2_y, sig2_z), 
                       method='L-BFGS-B', bounds=bnds, callback=callbackF, options={'disp': True, 'maxiter': 100})
   
    theta1 = res.x
    var = Mstep_var(theta1, func_top, func_low, Xy, Y, Zsamples, Xz, Z)
    
    return (theta1, var[0], var[1])



def callbackF(Xi):
    ''' callback function to display information for the optimiser
    '''
    global Nfeval
    print '\t Iter {0:4d}:  Para values: {1:}'.format(Nfeval, Xi)
    Nfeval += 1

def Mstep_theta_obj_wrapper(theta, *args):
    ''' The wrapper to pass variable arguments to Mstep_theta_obj
    '''

    return Mstep_theta_obj(theta, *args)
    #return Mstep_theta_obj(theta, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])

    
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
        #print '\t theta = {0:}'.format(theta)
        z_func = func_low(theta, np.array(Xy[n,:]))

        for i in range(n_mc):
            err = Zmc[i,:,n] - z_func
            neg_lnlik += Wmc[i,n]* 0.5* np.sum( beta* (np.array(err)**2) )
            y_func = func_top(theta, Xy[n,:], Zmc[i,:,n])
            err = Y[n,:] - y_func
            neg_lnlik += Wmc[i,n]* 0.5*(alpha*err*err)

    # For low-level X-Z data
    d_z = len(Xz)

    # Data in Xz and Z are saved in lists, each item in the lists represent one
    #    Z-variable
    for i in range (d_z):
        Xtmp = Xz[i]
        Ztmp = Z[i]
        n_dat_Xz = Xtmp.shape[0]   
        for n in range(n_dat_Xz):
            z_func = func_low(theta, np.array(Xtmp[n,:]))
            err = Ztmp[n] - z_func[:,i]
            #print '\t err = {0:}'.format(err)
            neg_lnlik += 0.5* np.sum( beta[i]* (np.array(err)**2) )
            
    if np.isfinite(neg_lnlik):
        return neg_lnlik
    else:
        return 1e10

## not verified yet
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

def calc_grad_theta_Z(func, theta, X, Z):
    '''Calculate the gradient of <func> w.r.t parameters <theta> & <Z>, with given <X> and <Z>
    using finite difference
    '''

    grad_theta = calc_grad_theta(func, theta, X, Z)

    n_Z = Z.shape[1]
    n_dat = X.shape[0]

    f = func(theta, X, Z)
    if f.ndim == 1:
        d_f = 1
    else:
        d_f = f.shape[1]
    
    Z1 = np.zeros(Z[0,:].shape)
    grad = np.zeros((n_dat, d_f, n_Z))

    for i in range(n_dat):
        for j in range(n_Z):

            delta = Z[i,j]*1e-4
            np.copyto(Z1, Z[i,:])
            if np.abs(delta) < 1e-5:
                delta = 1e-5
            Z1[j] += delta

            f = func(theta, X[i,:], Z[i,:])
            f1 = func(theta, X[i,:], Z1)

            gd = (f1-f) / delta
            if np.isscalar(gd):
                grad[i,:,j] = gd
            else:
                grad[i,:,j] = gd[:,None]

    return (grad_theta, grad)

def calc_grad_theta(func, theta, X, Z=None):
    '''Calculate the gradient of <func> w.r.t parameters <theta> with given <X> and <Z>
    using finite difference
    '''
    n_theta = theta.shape[0]
    n_dat = X.shape[0]

    if Z is not None:
        f = func(theta, X, Z)
    else:
        f = func(theta, X)

    if f.ndim == 1:
        d_f = 1
    else:
        d_f = f.shape[1]
    
    theta1 = np.zeros(theta.shape)
    grad = np.zeros((n_dat, d_f, n_theta))

    for i in range(n_theta):
        delta = theta[i]*1e-4
        np.copyto(theta1,theta)
        if np.abs(delta) < 1e-5:
            delta = 1e-5
        theta1[i] += delta

        if Z is not None:
            f = func(theta, X, Z)
            f1 = func(theta1, X, Z)
        else:
            f = func(theta, X)
            f1 = func(theta1, X)

        gd = (f1-f) / delta
        if gd.ndim < 2:
            grad[:,:,i] = gd[:,None]
        else:
            grad[:,:,i] = gd

    return grad


def Mstep_var(theta, func_top, func_low, Xy, Y, Zsamples, Xz, Z):
    ''' The function to calculate the optimal values of the variance terms
    '''
    
    n_dat_Xy = Xy.shape[0]
    Zmc = Zsamples[0]
    Wmc = Zsamples[1]    
    n_mc = Zmc.shape[0]    

    sse_l = np.zeros(Zmc[0,:,0].shape)
    sse_h = np.zeros(Y[0,:].shape)
    sig2_z = np.zeros(sse_l.shape)
    sig2_y = np.zeros(sse_h.shape)
    
    # For high-level X-Y data
    for n in range(n_dat_Xy):
        z_func = func_low(theta, np.array(Xy[n,:]))

        for i in range(n_mc):
            err = Zmc[i,:,n] - z_func
            sse_l += Wmc[i,n]* np.squeeze(err)**2
            y_func = func_top(theta, Xy[n,:], Zmc[i,:,n])
            err = Y[n,:] - y_func
            sse_h += Wmc[i,n]* np.squeeze(err)**2

    sig2_y = sse_h / n_dat_Xy
    
    # For low-level X-Z data
    d_z = len(Xz)

    # Data in Xz and Z are saved in lists, each item in the lists represent one
    #    Z-variable
    for i in range (d_z):
        Xtmp = Xz[i]
        Ztmp = Z[i]
        n_dat_Xz = Xtmp.shape[0]   
        for n in range(n_dat_Xz):
            z_func = func_low(theta, np.array(Xtmp[n,:]))
            err = Ztmp[n] - z_func[:,i]
            sse_l[i] += err*err

        sig2_z[i] = sse_l[i] / (n_dat_Xz+n_dat_Xy)

    return (sig2_y, sig2_z)

def pred(func_top, func_low, X, theta, sig2_y, sig2_z, V):
    ''' Function to make prediction (both mean and variance
    Args:
    theta - parameter (estimated)
    V - the covariance matrix of theta
    '''

    n_dat = X.shape[0]
    n_theta = theta.shape[0]
    
    # step 1: predict for lower-level model
    
    Z_mean = func_low(theta, X)
    grad_low = calc_grad(func_low, theta, X)
    d_func_low = grad_low.shape[1] # grad_low in shape of [n_dat, d_func, n_paras]
    
    Z_cov = np.zeros((n_dat, d_func_low, d_func_low))
    for i in range(n_dat): 
        gd = np.mat(grad_low[i,:,:])
        #print gd
        Z = gd * V * np.matrix.transpose(gd) + np.diag(sig2_z)
        #print Z
        Z_cov[i,:,:] = np.array(Z)

        
    # step 2: predict for higher-level model
    
    Y_mean = func_top(theta, X, Z_mean)
    grad_high = calc_grad_theta_Z(func_top, theta, X, Z_mean)
    grad_high_theta = grad_high[0] # grad_high_theta in shape of [n_dat, d_func, n_paras]
    grad_high_Z = grad_high[1]     # grad_high_Z in shape of [n_dat, d_func, d_Z]
    d_func_high = grad_high_theta.shape[1] 

    Y_cov = np.zeros((n_dat, d_func_high, d_func_high))
    zero_mat = np.mat(np.zeros((V.shape[0], Z_cov.shape[1])))
    zero_mat_trans = np.matrix.transpose(zero_mat)
    for i in range(n_dat):
        Z = np.mat(Z_cov[i,:,:])
        block_mat = np.bmat([ [V, zero_mat], [zero_mat_trans, Z] ])
        gd = np.bmat( [ grad_high_theta[i,:,:], grad_high_Z[i,:,:] ] )
        Y = gd * block_mat * np.matrix.transpose(gd) + np.diag(sig2_y)
        print Y
        Y_cov[i,:,:] = np.array(Y)

    return (Z_mean, Z_cov, Y_mean, Y_cov)
    
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

    if Z.ndim == 1:
        Kcc = 10**Z[0]
        Klp = 10**Z[1]
    else:
        Kcc = 10**Z[:,0]
        Klp = 10**Z[:,1]

    Ksc_pred = phi_pro*rho_pro/rho_wat* Kcc + phi_lip*rho_lip/rho_wat* Klp + phi_wat
    Y =  Ksc_pred
        
    return np.log10(Y)


def testFunc_low(theta, X):
    ''' Function to predict the volumetric partition coefficient between corneocyte (and water)
    and that between lipid (and water)
    Here used as a test function of the low-level model
    Args:
      theta -- model parameters
      X -- lg10Kow, dim: [n_dat, 1]
    Rtns:
      Z -- dim: [n_dat, 2]
    '''

    #if theta.shape[0] < 3:
    #    print 'problem'

    theta = np.squeeze(np.asarray(theta))
    # to avoid some mysterious conversion of np array to np matrix by the optimiser
        
    a = np.exp(theta[0])
    b = np.exp(theta[1])
    c = np.exp(theta[2])

    n_dat = X.shape[0]

    Z = np.zeros((n_dat, 2))
    #Z = np.mat(Z)

    rho_pro = 1.37
    rho_lip = 0.90
    rho_wat = 1.00

#    for i in range(n_dat):
#        K_ow = 10**X[i,0]
#        Z[i, 0] = rho_pro/rho_wat* a*np.power(K_ow,b) # stratum corneum
#        if np.isfinite( np.power(K_ow,b) ) == False:
#            print 'problem power'
#        Z[i, 1] = rho_lip/rho_wat* np.power(K_ow,c) # lipid
#        
#    return np.log10(Z)

    for i in range(n_dat):
        lg10_K_ow = X[i,0]
        Z[i, 0] = np.log10(rho_pro/rho_wat) + np.log10(a) + b*lg10_K_ow # stratum corneum
        Z[i, 1] = np.log10(rho_lip/rho_wat) + c*lg10_K_ow # lipid
        
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
def calcHessVargs(func_post, paras, *args):
    ''' Function to calculate the Hessian of negative
        log posterior w.r.t. model parameters with variable arguments
    '''

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
                t1 = func_post(p1, *args)

                p2 = np.copy(paras)                
                p2[i] += delta_paras[i]
                p2[j] -= delta_paras[j]
                t2 = func_post(p2, *args)

                p3 = np.copy(paras)                
                p3[i] -= delta_paras[i]
                p3[j] += delta_paras[j]
                t3 = func_post(p3, *args)

                p4 = np.copy(paras)                
                p4[i] -= delta_paras[i]
                p4[j] -= delta_paras[j]
                t4 = func_post(p4, *args)

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
