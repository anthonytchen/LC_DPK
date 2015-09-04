def configDPK():
    '''Function to prepare the input parameters for dermatokinetic simulation of CAFFEINE'''

    from DPKInputParas import *  
    cfg = DPKInputParas()

    # skin geometric parameters, only change those are not default values
    cfg.n_layer_x_sc = 16+5
    cfg.x_len_viaepd = 100e-6
    cfg.x_len_dermis = 1200e-6
    cfg.n_grids_x_ve = 10
    cfg.n_grids_x_de = 10
    cfg.n_layer_y_sc = 1
    cfg.offset_y_sc = 40.04e-6
    cfg.b_inf_src = 0

    # solute parameters
    cfg.MW = 194.19
    cfg.K_ow = 10**-0.07
    cfg.pKa = 14.0 # to do
    cfg.AorB = 'B' # (B)ase to be checked
    cfg.frac_non_ion = 0.31 # to be checked
    cfg.frac_unbound = 0.65 # Ref.: Blanchard J (1982). Protein binding of caffeine in young and elderly males, Journal of Pharmaceutical Sciences, 71: 1415–1418.
    cfg.par_dermis2blood = 1.0/(10**0.04) # log_P blood:skin is 0.04
    cfg.k_clear_blood = 2.69e-6 # 9.7 L/hr = 2.69e-6 m3/s, Ref.: Carbo M, et al. (1989). Effect of quinolones on caffeine disposition, Clinical Pharmacology & Therapeutics, 45: 234-240.

    # vehicle parameters

    return cfg
