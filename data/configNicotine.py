def configDPK():
    '''Function to prepare the input parameters for dermatokinetic simulation of NICOTINE'''

    from DPKInputParas import *  
    cfg = DPKInputParas()

    # skin geometric parameters, only change those are not default values
    cfg.n_layer_x_sc = 23 # 22.6 \mu m / 0.0875
    cfg.x_len_viaepd = 100e-6
    cfg.x_len_dermis = 1200e-6
    cfg.n_grids_x_ve = 10
    cfg.n_grids_x_de = 10
    cfg.n_layer_y_sc = 1
    cfg.offset_y_sc = 40.04e-6
    cfg.b_inf_src = 0

    # solute parameters
    cfg.MW = 162.23
    cfg.K_ow = 10**1.17
    cfg.pKa = 3.12 
    cfg.AorB = 'B' # (B)ase
    cfg.frac_non_ion = 0.31
    cfg.frac_unbound = 0.95
    cfg.par_dermis2blood = 1.0/(10**0.04) # log_P blood:skin is 0.04
    cfg.k_clear_blood = 23.3e-6 # 1400 ml/min = 23.3e-6 m3/s, 1250 = 20.8e-6, 1540 = 25.7e-6

    # vehicle parameters

    return cfg
