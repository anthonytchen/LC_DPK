import os
import time
import warnings
import numpy as np

import Config as conf
import Chemical as chem
import Skin_Setup as skin

def dpk_permeability(cfn = 'data/Caffeine.cfg'):
    ''' Function to compute the steady-state permeability
    Args:
       cfn -- the file for configuration
    Rtn:
       perm -- permeability if steady-state reached; 
               otherwise return False and issue warning message
    '''

    # The following setting for folders is only meant to fulfill
    #   the requirement of the underlying C++ classes.
    #   For simple permeability calculation they are not really used.
    
    dir = 'Dir'
    os.system(' mkdir -p ./' + dir)

    fn_conc = dir + '/conc'
    fn_coord_x = dir + '/coord_x'
    fn_coord_y = dir + '/coord_y'
    fn_amount = dir + '/amount'
    fn_flux = dir + '/flux'
   
    _conf = conf.Config()
    _conf.ReadConfigFile(cfn)

    _chem = chem.Chemical()
    _chem.InitConfig(_conf)

    _skin = skin.Skin_Setup()
    _skin.InitConfig(_chem, _conf)

    # todo: test StraCorn::setGridsProperties()

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

