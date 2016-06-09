# The module that contains computational routines for hybrid model
#   identification and uncertainty quantification

import numpy as np
import math
import scipy.sparse as sp
import scipy.sparse.linalg as spln
import matplotlib.pyplot as plt

# Example:

# todo list:


###########################################################
def Estep(func_top, func_low, X, Y, theta, sig2_y, sig2_z, N=100):
    ''' The E-step of the Monte Carlo EM algorithm
    Args:
    - func_top: top level function;  y = f(x,z,theta) + epsilon
    - func_low: low level function;  z = h(x,theta)   + zeta
    - X: data for input variables
    - Y: data for top level output variables
    - theta: model parameters
    - sig2_y: variance of top level function noise (epsilon)
    - sig2_z: variance of low level function noise (zeta)
    - N: number of MC samples
    Rtns:
    '''

    n_dat = X.shape[0]
    d_y = Y.shape[1]
    d_z = func_low(X[0,:], theta).shape[0] # dimension of z

    Ymc = np.zeros((N, d_y))
    Z = np.zeros((N, d_z, n_dat)) # MC samples
    W = np.zeros((N, n_dat)) # weights

    sig_y = np.sqrt(sig2_y)
    sig_z = np.sqrt(sig2_z)

    for n in range(n_dat):

        z_n = func_low(X[n,:], theta)
        Z[:,:,n] = np.tile(z_n,(N,1))

        # MC sampling
        for d in range(d_z):
            rd = np.random.normal(0, sig_z[d], N)
            Z[:,d,n] += rd

        Ymc = func_top( X[n,:], Z[:,:,n], theta )
        s = .0;
        for i in range(N):
            W[i,n] = np.exp( lognormpdf( Ymc[i,:], Y[n,:], np.diag(sig2_y) ) )
            s += W[i,n]
        W[:,n] /= s # normalising the weights

    return (Z, W)


###########################################################
def Mstep_main():


def Mstep_theta_obj(func_top, func_low, X, Y, theta, sig2_y, sig2_z, Z, W):
    ''' The objective function (negative log-likelihood) for optimising theta
    Terms that are not dependent on theta (thus constant as far as optimisation is concerned)
    are not calculated.
    '''

    n_dat = X.shape[0]
    d_y = Y.shape[1]
    n_mc = Z.shape[0]
    d_z = Z.shape[1]

    beta = 1.0/sig2_z
    alpha = 1.0/sig2_y

    neg_lnlik = .0
    for n in range(n_dat):
        
        z_func = func_low(X[n,:], theta)
        y_func = func_top(X[n,:], Z[:,:,n], theta )

        for i in range(n_mc):
            err = Z[i,:,n] - z_func
            neg_lnlik += W[i,n]* 0.5*(beta*err*err)
            err = Y[n,:] - y_func
            neg_lnlik += W[i,n]* 0.5*(alpha*err*err)

    return neg_lnlik


def Mstep_theta_grad(func_top, func_low, X, Y, theta, sig2_y, sig2_z, Z, W):
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
        
        z_func = func_low(X[n,:], theta)
        y_func = func_top(X[n,:], Z[:,:,n], theta )
        gd_h = calc_grad(func_low, X[n,:], theta)

        for i in range(n_mc):
            err = Z[i,:,n] - z_func
            gd1 = beta*err* gd_h
            err = Y[n,:] - y_func
            gd2 = alpha*err* calc_grad(func_top, X[n,:], Z[i,:,n], theta)
            grad += W[i,n]* (gd1+gd2)

    return neg_lnlik

def calc_grad(func, X, Z=None, thetea):
    '''Calculate the gradient of <func> w.r.t parameters <theta> with given <X> and <Z>
    '''
    

def Mstep_var(func_top, func_low, X, Y, theta, Z, W):
    ''' The function to calculate the optimal values of the variance terms
    '''

    n_dat = X.shape[0]
    d_y = Y.shape[1]
    n_mc = Z.shape[0]
    d_z = Z.shape[1]

    beta = np.zeros((d_z,1))
    alpha = np.zeros((d_y,1))

    for n in range(n_dat):
        
        z_func = func_low(X[n,:], theta)
        y_func = func_top(X[n,:], Z[:,:,n], theta )

        for i in range(n_mc):
            beta += W[i,n]* ( Z[i,:,n] - z_func )**2
            alpha += W[i,n]* ( Y[n,:] - f_func )**2

    beta = N/beta;
    alpha = N/alpha;      

    return (alpha, beta)



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
