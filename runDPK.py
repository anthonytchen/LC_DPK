class strdConc:
    ''' Class definition for structured concentration and mass profiles at a particular time point'''
    def __init__(self):

        # concentration values
        self.vehicle = []
        self.stracorn2d_coordx = []
        self.stracorn2d_coordy = []
        self.stracorn2d_data = []
        self.stracorn1d = []
        self.viaepd = []
        self.dermis = []
        self.blood = []

        # mass values
        self.mass_applied = [] # mass dose transported to skin so far
        self.mass_stracorn_total = []
        self.mass_stracorn_lipid = []
        self.mass_stracorn_corneocyte = []
        self.mass_viaepd = []
        self.mass_dermis = []
        self.mass_blood = []
        self.mass_cleared = []


def plot_stracorn(dat, fn_str):

    import numpy as np
    import matplotlib.pyplot as plt

    nx_layers = 23

    nx_grids_lipid = 2; # of x-grids for lipid layer, 2
    nx_grids_cc = 4; # of x-grids for corneocyte layer, 4
    ny_grids_lipid = 2; # of y-grids for lipid layer, 2
    ny_grids_cc_dn = 2; # of y-grids for dn-part of the offset corneocyte layer, 2

    plt.subplots(figsize=(8,8))

    # Plot the lines that mark the boundary of corneocytes
    for i in range(nx_layers): # starting from 0
        
        if np.remainder(i, 2) == 0: # if i is even number
            x_st = ny_grids_lipid - 1 - 0.5
            x_ed = ny_grids_lipid - 1 + ny_grids_cc_dn + ny_grids_lipid + ny_grids_cc_dn*8 - 0.5
            y_st = i*(nx_grids_lipid+nx_grids_cc) + nx_grids_lipid - 0.5
            y_ed = y_st + nx_grids_cc
            plt.plot([x_st, x_ed], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st, x_ed], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st, x_st], [y_st, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_ed, x_ed], [y_st, y_ed], color='k', linestyle='-', linewidth=1)

        else:
            x_ed1 = ny_grids_lipid - 1 + ny_grids_cc_dn - 0.5
            x_st2 = ny_grids_lipid - 1 + ny_grids_cc_dn + ny_grids_lipid - 0.5
            x_ed2 = x_st2 + ny_grids_cc_dn*8 + ny_grids_lipid - 1
            y_st = i*(nx_grids_lipid+nx_grids_cc) + nx_grids_lipid - 0.5
            y_ed = y_st + nx_grids_cc
            plt.plot([-0.5, x_ed1], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_ed2], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([-0.5, x_ed1], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_ed2], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            
            plt.plot([x_ed1, x_ed1], [y_st, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_st2], [y_st, y_ed], color='k', linestyle='-', linewidth=1)

    plt.imshow(dat, interpolation='none', aspect='auto')
    cb = plt.colorbar(shrink=0.98, pad=0.05)
    cb.set_clim(0, 2900)
    plt.tight_layout()
    plt.axis('off')
    plt.savefig(fn_str, bbox_inches='tight')
    plt.close()
    #plt.show()

def calConcProfile(conc_obj, skin_obj, conc_raw, coord_x, coord_y):
    ''' The function to work out concentration at different depth    '''

    import numpy as np

    # calulcate concentration in vehicle
    idx = 0
    conc_obj.vehicle = np.matrix( [coord_x[ idx : idx+skin_obj.m_dim_vh ], conc_raw[ idx : idx+skin_obj.m_dim_vh ]] )
    idx += skin_obj.m_dim_vh


    # calculate concentration in stratum corneum
    n_layer_x_sc = skin_obj.getNLayerXSc(); # layer of cells in x direction in SC
    a1 = np.zeros(n_layer_x_sc);
    skin_obj.get1DCoordSC(a1);
    a2 = np.zeros(n_layer_x_sc);
    skin_obj.get1DConcSC(a2);
    conc_obj.stracorn1d = np.matrix( [a1, a2] );

    conc_obj.stracorn2d_coordx = np.array( coord_x[ idx : idx+skin_obj.m_dim_sc ] )
    conc_obj.stracorn2d_coordy = np.array( coord_y[ idx : idx+skin_obj.m_dim_sc ] )
    n_x_sc = skin_obj.getNGridsXSc()
    n_y_sc = skin_obj.getNGridsYSc()
    conc_obj.stracorn2d_data = np.zeros( (n_x_sc, n_y_sc) );
    for i in range(n_x_sc):
        st = idx + i*n_y_sc;
        ed = idx + (i+1)*n_y_sc;
        conc_obj.stracorn2d_data[i,:] = conc_raw[st:ed];

    idx += skin_obj.m_dim_sc;

    # calculate concentration in viable epidermis
    coord = coord_x[ idx : idx+skin_obj.m_dim_ve ] + coord_x[ idx+1 : idx+skin_obj.m_dim_ve+1 ]
    coord /= 2; # use the middle point as the coordinate
    conc_obj.viaepd = np.matrix( [coord, conc_raw[ idx : idx+skin_obj.m_dim_ve ]] )
    idx += skin_obj.m_dim_ve;

    # calculate concentration in dermis
    coord = coord_x[ idx : idx+skin_obj.m_dim_de ] + coord_x[ idx+1 : idx+skin_obj.m_dim_de+1 ]
    coord /= 2; # use the middle point as the coordinate
    conc_obj.dermis = np.matrix( [coord, conc_raw[ idx : idx+skin_obj.m_dim_de ]] )
    idx += skin_obj.m_dim_de;

    # calculate concentration in blood
    conc_obj.blood = np.matrix( [coord_x[ idx : idx+skin_obj.m_dim_bd ], conc_raw[ idx : idx+skin_obj.m_dim_bd ]] )

    # calculate amount of solute in differnet layers
    
    n_amount = 9
    dat = np.zeros(n_amount);
    skin_obj.getLayersAmount(dat);
    
    conc_obj.mass_applied = dat[0]-dat[1]
    conc_obj.mass_stracorn_total = dat[2]
    conc_obj.mass_stracorn_lipid = dat[3]
    conc_obj.mass_stracorn_corneocyte = dat[4]
    conc_obj.mass_viaepd = dat[5]
    conc_obj.mass_dermis = dat[6]
    conc_obj.mass_blood = dat[7]
    conc_obj.mass_cleared = dat[8]

def calcAUC(t, conc):
    '''Calculate AUC, t_max, C_max from time-profile of concentration values
    '''
    import numpy as np

    auc = .0
    for i in range(len(t)-1):
        auc += (conc[i] + conc[i+1])/2 * (t[i+1]-t[i])
    
    t_max = np.amax(t)
    t_interp = np.linspace(0, t_max, 1000)
    conc_interp = np.interp(t_interp, t, conc)

    t_max = t_interp[np.argmax(conc_interp)]
    C_max = np.amax(conc_interp)

    return (auc, t_max, C_max)

def runDPKfunc(partition_vehicle, diffu_vehicle, dx_vehicle, t_range, dose_factor=1, clear_blood=None):
    '''Run DPK using provided partition and diffusivity of

    partition_vehicle = 1 : partition coefficient of vehicle to water
    diffu_vehicle = 9.727e-10 : diffusivity of solute in vehicle
    dx_vehicle = 100e-6 : depth of vehicle in meter
    t_range = range(0, t_end, t_inv) : the time points for which simulation is needed
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    import Skin as skin
    import Chemical as chem

    from data import configNicotine
    reload(configNicotine)
    cfg = configNicotine.configDPK()
    if clear_blood != None:
        cfg.k_clear_blood = clear_blood	

    area_vehicle = dose_factor*3.5*1e-4; # cm2 patch, represented in m2 
    # dx_vehicle = 100e-6  # thickness of vehicle in meter
    conc_vehicle = dose_factor*15*1e-6/(area_vehicle*dx_vehicle); # mg in cm2 patch, in kg/m3
 
    _chem = chem.Chemical();
    _skin = skin.Skin();

    _chem.Init(cfg.MW, cfg.K_ow, cfg.pKa, cfg.frac_non_ion, cfg.frac_unbound, cfg.AorB);
    _skin.Init(_chem, conc_vehicle, diffu_vehicle, partition_vehicle, 
               cfg.par_dermis2blood, cfg.k_clear_blood,
               dx_vehicle, area_vehicle, cfg.x_len_viaepd, cfg.x_len_dermis,
               cfg.n_layer_x_sc, cfg.n_layer_y_sc, cfg.n_grids_x_ve, cfg.n_grids_x_de, 
               cfg.offset_y_sc, cfg.b_inf_src);

    conc = np.zeros(_skin.m_dim_all);
    _skin.getGridsConc(conc);

    coord_x = np.zeros(_skin.m_dim_all);
    coord_y = np.zeros(_skin.m_dim_all);
    _skin.getXCoord(coord_x);
    _skin.getYCoord(coord_y);

    stdConc = strdConc();

    blood_conc = np.zeros(len(t_range))
    vehicle_conc = np.zeros(len(t_range))

    mass_applied = np.zeros(len(t_range))
    mass_stracorn_total = np.zeros(len(t_range))
    mass_stracorn_lipid = np.zeros(len(t_range))
    mass_stracorn_corneocyte = np.zeros(len(t_range))
    mass_viaepd = np.zeros(len(t_range))
    mass_dermis = np.zeros(len(t_range))
    mass_blood = np.zeros(len(t_range))
    mass_cleared = np.zeros(len(t_range))

    calConcProfile(stdConc, _skin, conc, coord_x, coord_y)
    blood_conc[0] = stdConc.blood[1,:];
    vehicle_conc[0] = stdConc.vehicle[1,:];

    mass_applied[0] = stdConc.mass_applied
    mass_stracorn_total[0] = stdConc.mass_stracorn_total
    mass_stracorn_lipid[0] = stdConc.mass_stracorn_lipid
    mass_stracorn_corneocyte[0] = stdConc.mass_stracorn_corneocyte
    mass_viaepd[0] = stdConc.mass_viaepd
    mass_dermis[0] = stdConc.mass_dermis
    mass_blood[0] = stdConc.mass_blood
    mass_cleared[0] = stdConc.mass_cleared

    b_1st_save = 1

    for i in range(len(t_range)-1):

        t = t_range[i]
        t_next = t_range[i+1]

        #t1 = t/3600
        #print( 'time = %.3e' % t1 )

        _skin.diffuseMoL(t, t_next)
        _skin.saveGrids(b_1st_save, 'python_rlt')
        if b_1st_save:
             b_1st_save = 0

        _skin.getGridsConc(conc);
        calConcProfile(stdConc, _skin, conc, coord_x, coord_y)

        #flux1 = _skin.compFlux_2sc();
        #flux2 = _skin.compFlux_sc2down();
        #flux3 = _skin.compFlux_ve2down();
        #flux4 =  _skin.compFlux_de2down();

        #print( '\t flux_2sc = %.3e, flux2ve = %e, flux_2de = %.3e, flux_2bd = %.3e\n' % (flux1, flux2, flux3, flux4) )

        fn_str = 'stracorn_' + str(t_next) + '.png'
	plot_stracorn( stdConc.stracorn2d_data, fn_str )

        # plt.plot( np.hstack( (stdConc.stracorn1d[0,:], stdConc.viaepd[0,:], stdConc.dermis[0,:]) ), 
        #          np.hstack( (stdConc.stracorn1d[1,:], stdConc.viaepd[1,:], stdConc.dermis[1,:]) ), 'r--o' )
        blood_conc[i+1] = stdConc.blood[1,:];
        vehicle_conc[i+1] = stdConc.vehicle[1,:]

        mass_applied[i+1] = stdConc.mass_applied
        mass_stracorn_total[i+1] = stdConc.mass_stracorn_total
        mass_stracorn_lipid[i+1] = stdConc.mass_stracorn_lipid
        mass_stracorn_corneocyte[i+1] = stdConc.mass_stracorn_corneocyte
        mass_viaepd[i+1] = stdConc.mass_viaepd
        mass_dermis[i+1] = stdConc.mass_dermis
        mass_blood[i+1] = stdConc.mass_blood
        mass_cleared[i+1] = stdConc.mass_cleared

        if abs(t_next-24*3600) < 0.1 : # remove vehicle
            _skin.removeVehicle()
        #if abs(t_next % (24*3600)) < 0.1 : # reset vehicle
        #    _skin.resetVehicle(conc_vehicle, partition_vehicle, diffu_vehicle)

    return (blood_conc, vehicle_conc, mass_applied, mass_stracorn_total, mass_stracorn_lipid, mass_stracorn_corneocyte, mass_viaepd, mass_dermis, mass_blood, mass_cleared)


def objFunCalib(paras, t_range, blood_conc_data, dose_factor=1.0):
    ''' Function to calculate the objective function to be minimised for use with model calibration'''

    import numpy as np

    partition_vehicle = 10**paras[0]
    #diffu_vehicle = 10**paras[1]
    #dx_vehicle = 100*1e-6

    #partition_vehicle = 1
    diffu_vehicle = 9.727e-10
    dx_vehicle = paras[1]*1e-6

    conc = runDPKfunc( partition_vehicle, diffu_vehicle, dx_vehicle, t_range, dose_factor )
    blood_conc = conc[0]*1e6 # converting from kg/m^3 to ng/ml

    # calculate the linear trapesoidal AUC
    auc_simu = np.zeros(t_range.size-1)
    auc_data = np.zeros(t_range.size-1)
    for i in range(len(t_range)-1):
        tt = ( t_range[i+1] - t_range[i] ) / 3600
        auc_simu[i] = 0.5*tt* ( blood_conc[i] + blood_conc[i+1] )
        auc_data[i] = 0.5*tt* ( blood_conc_data[i] + blood_conc_data[i+1] )

    #print blood_conc
    #print auc_simu
    #print auc_data
    #err = np.sqrt( np.mean( np.square( blood_conc_data - blood_conc ) ) )
    err = np.sum( np.square( auc_data - auc_simu ) )
    #print err
    return err
