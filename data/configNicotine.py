def configDPK():
    '''Function to prepare the input parameters for dermatokinetic simulation of NICOTINE'''

    from DPKInputParas import *  
    cfg = DPKInputParas()

    # skin geometric parameters, only change those are not default values
    cfg.n_layer_x_sc = 16+7

    # solute parameters
    cfg.MW = 162.23
    cfg.K_ow = 10**1.17
    cfg.pKa = 3.12 
    cfg.AorB = 'B' # (B)ase
    cfg.frac_non_ion = 0.31
    cfg.frac_unbound = 0.95
    cfg.par_dermis2blood = 1.0/(10**0.04) # log_P blood:skin is 0.04
    cfg.k_clear_blood = 23.3e-6 # 1400 ml/min = 23.3e-6 m3/s

    # vehicle parameters

    return cfg
