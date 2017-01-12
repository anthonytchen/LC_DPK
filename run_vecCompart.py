import os
import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing

#from multiprocessing import Pool

import Config as conf
import Chemical as chem
import Skin_Setup as skin

# example:
#    For caffeine:
#        dpk_perm( np.array([-1, -1, -1, -1]), np.array([194.19, 10**-0.07]) )

def main():
    import imp
    hybmdl = imp.load_source('hybmdl', './uncertainty/hybmdl.py')

    dat_Plp = np.matrix( np.loadtxt("./uncertainty/Kow_Plp_lg10.txt") )
    dat_Ksc = np.matrix( np.loadtxt("./uncertainty/Kow_Ksc_lg10.txt") )
    dat_Kcc = np.matrix( np.loadtxt("./uncertainty/Kow_Kcc_lg10.txt") )
    paras0 = np.log( np.array([4.2, 0.31, 0.69]) )
    bnds = ((-10, 10), (-10, 10), (-10, 10))

    sig2_y = np.array([0.05])
    sig2_z = np.array([0.05, 0.05])

    Xy = dat_Ksc[:,0].reshape((-1, 1))
    Y = dat_Ksc[:,1].reshape((-1, 1))
    Xz = ( dat_Plp[:,0].reshape((-1, 1)), dat_Kcc[:,0].reshape((-1, 1)) )
    Z = ( dat_Plp[:,1].reshape((-1, 1)), dat_Kcc[:,1].reshape((-1, 1)) )

    paras = np.empty_like (paras0)
    np.copyto(paras, paras0)

    mdl = hybmdl.PluginMain(hybmdl.testFunc_top_plugin, hybmdl.testFunc_low, Xy, Y, Xz, Z, paras0, sig2_y, sig2_z, 10, bnds)

    theta = mdl[0]
    sig2_z = mdl[2]
    V = mdl[3]
    
    # predict for permeability

    # caffeine, nicotine
    mw = [194.19, 162.23]
    lg10kow = [-0.07, 1.17]
    X = np.array([mw, lg10kow]).T

    results = Parallel(n_jobs=2)(delayed(
        hybmdl.pred)(perm_Kw_cc_lip, Kw_cc_lip, X[i,:].reshape((1,-1)), theta, np.array([0]), sig2_z, V) for i in range(2) )
    
    #pool = Pool() # multiprocessing
    #prd1 = pool.apply_async( hybmdl.pred(perm_Kw_cc_lip, Kw_cc_lip, X[0,:].reshape((1,-1)), theta, np.array([0]), sig2_z, V) )
    #prd2 = pool.apply_async( hybmdl.pred(perm_Kw_cc_lip, Kw_cc_lip, X[1,:].reshape((1,-1)), theta, np.array([0]), sig2_z, V) )
    
    # prd = hybmdl.pred(perm_Kw_cc_lip, Kw_cc_lip, X, theta, np.array([0]), sig2_z, V)

    #todo: finally check the transform-based uncertainty results against simple sensitivity analysis!
    return (mdl, results)


def dpk_perm(sc_ptys, chem_ptys):
    ''' Function to compute the steady-state permeability
    Args:
       sc_ptys -- the properties of the stratum corneum: [lip_Kw, lip_D, cc_Kw, cc_D]
                  negative value means to calcualate by using the built-in QSPR
       chem_ptys -- the properties of the chemical: [MW, Kow]
                  *** NOTE: this is Kow, not log10_Kow! ***
    Rtn:
       perm -- permeability if steady-state reached; 
               otherwise return False and issue warning message
    '''

    # The following setting for folders is only meant to fulfill
    #   the requirement of the underlying C++ classes.
    #   For simple permeability calculation they are not really used.
    cfn = 'data/CommonPermeability.cfg'
    
    dir = 'Dir'
    os.system(' mkdir -p ./' + dir)

    fn_conc = dir + '/conc'
    fn_coord_x = dir + '/coord_x'
    fn_coord_y = dir + '/coord_y'
    fn_amount = dir + '/amount'
    fn_flux = dir + '/flux'
   
    _conf = conf.Config()
    _conf.ReadConfigFile(cfn)
    _conf.m_mw = chem_ptys[0]
    _conf.m_K_ow = chem_ptys[1]

    _chem = chem.Chemical()
    _chem.InitConfig(_conf)

    _skin = skin.Skin_Setup()
    _skin.InitConfig(_chem, _conf)
    _skin.setScProperties(sc_ptys[0], sc_ptys[1], sc_ptys[2], sc_ptys[3])


    t_end = 3600*100
    t_inv = 3600*4
    
    b_1st_save = True;

    _skin.saveCoord( fn_coord_x, fn_coord_y )
    _skin.saveGrids(b_1st_save, fn_conc)
    _skin.compCompartAmount()
    _skin.saveAmount(b_1st_save, fn_amount)
    b_1st_save = not b_1st_save
      
    t_simu = .0
    perm = -1e10
    while t_simu < t_end:
        start = time.clock()

        _skin.diffuseMoL(t_simu, t_simu+t_inv)
        _skin.saveGrids(b_1st_save, fn_conc)
        _skin.compCompartAmount()
        _skin.saveAmount(b_1st_save, fn_amount)

        end = time.clock()
        cpu_time_used = end - start
        print 'Simulation time {0:.2f} hr, cpu time  {1:.2f} s'.format( (t_simu+t_inv)/3600, cpu_time_used)
                                                                    
        t_simu += t_inv

        # When flux reaches steady state, break
        #  only applicable for permeability calculation
        fluxP = skin.new_doubleP() # create a swig pointer to double to be used in _skin.saveFlux
        if  _skin.saveFlux(b_1st_save, fn_flux, fluxP) and _conf.m_bInfSrc:
            flux = skin.doubleP_value(fluxP)
            perm = flux/_conf.m_conc_vehicle
            print 'Steady-state flux = {0:.5e}, permeability = {1:.5e}'.format( flux, perm )
            break
  
    _skin.Release();


    if perm > 0:
        return perm
    else:
        warnings.warn("Steady state not reached; consider increasing the simulation time")
        return False

def perm_Kw_cc_lip(theta, X, Kw_skin):
    '''
    '''
    mw = X[0,0]
    lg10Kow = X[0,1]

    Kw = np.squeeze(np.asarray(Kw_skin))
    #Kw = Kw_cc_lip(theta, lg10Kow)    
    #Kw = np.squeeze(Kw)
    Kw_cc = 10**Kw[0]
    Kw_lip = 10**Kw[1]

    perm = dpk_perm( np.array([Kw_lip, -1, Kw_cc, -1]), np.array([mw, 10**lg10Kow]) )
    return np.log10(perm).reshape((1,-1))
                     
    
def Kw_cc_lip(theta, X):
    ''' Function to predict the volumetric partition coefficient between corneocyte (and water)
    and that between lipid (and water)
    Here used as a test function of the low-level model
    Args:
      theta -- model parameters
      X -- [mw, lg10Kow], dim: [n_dat, 2]
    Rtns:
      Z -- dim: [n_dat, 2]
    '''

    theta = np.squeeze(np.asarray(theta))
    # to avoid some mysterious conversion of np array to np matrix by the optimiser
        
    a = np.exp(theta[0])
    b = np.exp(theta[1])
    c = np.exp(theta[2])

    #print X
    #X = np.mat(X)
    n_dat = X.shape[0]

    Z = np.zeros((n_dat, 2))
    #Z = np.mat(Z)

    rho_pro = 1.37
    rho_lip = 0.90
    rho_wat = 1.00
    #print X
    #print Z

    for i in range(n_dat):
        lg10_K_ow = X[i,1]
        Z[i,0] = np.log10(rho_pro/rho_wat) + np.log10(a) + b*lg10_K_ow # corneocyte
        Z[i,1] = np.log10(rho_lip/rho_wat) + c*lg10_K_ow # lipid
        
    #print Z
    return Z


