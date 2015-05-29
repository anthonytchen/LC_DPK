class strdConc:
    ''' Class definition for structured concentration profiles
    '''
    def __init__(self):
        self.vehicle = []
        self.stracorn2d_coordx = []
        self.stracorn2d_coordy = []
        self.stracorn2d_data = []
        self.stracorn1d = []
        self.viaepd = []
        self.dermis = []
        self.blood = []

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
    for i in range(0, n_x_sc):
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

def runDPKfunc(partition_vehicle, diffu_vehicle, t_range):
    '''Run DPK using provided thickness of vehicle and diffusivity
    partition_vehicle = 1 : partition coefficient of vehicle to water
    diffu_vehicle = 5e-14 : diffusivity of solute in vehicle
    t_range = range(0, t_end, t_inv) : the time points for which simulation is needed
    '''
    import numpy as np

    import Skin as skin
    import Chemical as chem

    n_layer_x_sc = 16+5;
    n_layer_y_sc = 1;
    n_grids_x_ve = 10;
    n_grids_x_de = 10;
    offset_y_sc = 0;
    b_inf_src = 0;

    # Nicotine
    MW = 162.23;
    K_ow = pow(10, 1.17);
    pKa = 3.12; # a base

    area_vehicle = 1*3.5*1e-4; # cm2 patch, represented in m2 
    dx_vehicle = 100e-6  # thickness of vehicle in meter
    conc_vehicle = 1*15*1e-6/(area_vehicle*dx_vehicle); # mg in cm2 patch, in kg/m3
 
    _chem = chem.Chemical();
    _skin = skin.Skin();

    _chem.Init(MW, K_ow, pKa, 'B'); # the last letter denotes acid (A) or base (B)
    _skin.Init(_chem, conc_vehicle, diffu_vehicle, partition_vehicle, dx_vehicle, area_vehicle, 
               n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de, offset_y_sc, b_inf_src);

    conc = np.zeros(_skin.m_dim_all);
    _skin.getGridsConc(conc);

    coord_x = np.zeros(_skin.m_dim_all);
    coord_y = np.zeros(_skin.m_dim_all);
    _skin.getXCoord(coord_x);
    _skin.getYCoord(coord_y);

    stdConc = strdConc();

    blood_conc = np.zeros(t_range.size)
    vehicle_conc = np.zeros(t_range.size)

    calConcProfile(stdConc, _skin, conc, coord_x, coord_y)
    blood_conc[0] = stdConc.blood[1,:];
    vehicle_conc[0] = stdConc.vehicle[1,:];

    for i in range(0, t_range.size-1, 1):

        t = t_range[i]
        t_next = t_range[i+1]

        #print( 'time = %.3e' % t )

        _skin.diffuseMoL(t, t_next)
        _skin.getGridsConc(conc);
        calConcProfile(stdConc, _skin, conc, coord_x, coord_y)

        #flux1 = _skin.compFlux_2sc();
        #flux2 = _skin.compFlux_sc2down();
        #flux3 = _skin.compFlux_ve2down();
        #flux4 =  _skin.compFlux_de2down();

        #print( '\t flux_2sc = %.3e, flux2ve = %e, flux_2de = %.3e, flux_2bd = %.3e\n' % (flux1, flux2, flux3, flux4) )

        # plt.imshow( stdConc.stracorn2d_data )
        # plt.plot( np.hstack( (stdConc.stracorn1d[0,:], stdConc.viaepd[0,:], stdConc.dermis[0,:]) ), 
        #          np.hstack( (stdConc.stracorn1d[1,:], stdConc.viaepd[1,:], stdConc.dermis[1,:]) ), 'r--o' )
        blood_conc[i+1] = stdConc.blood[1,:];
        vehicle_conc[i+1] = stdConc.vehicle[1,:]

        if abs(t_next-24*3600) < 0.1 : # remove vehicle
            _skin.removeVehicle()

    return (blood_conc, vehicle_conc)


def objFunCalib(paras, t_range, blood_conc_data):
    ''' Function to calculate the objective function to be minimised for use with model calibration'''

    import numpy as np

    partition_vehicle = 10**paras[0]
    diffu_vehicle = 10**paras[1]

    conc = runDPKfunc( partition_vehicle, diffu_vehicle, t_range )
    blood_conc = conc[0]*1e6 # converting from kg/m^3 to ng/ml

    # calculate the linear trapesoidal AUC
    auc_simu = np.zeros(t_range.size-1)
    auc_data = np.zeros(t_range.size-1)

    for i in range(0, t_range.size-1, 1):
        tt = ( t_range[i+1] - t_range[i] ) / 3600
        auc_simu[i] = 0.5*tt* ( blood_conc[i] + blood_conc[i+1] )
        auc_data[i] = 0.5*tt* ( blood_conc_data[i] + blood_conc_data[i+1] )

    #print blood_conc
    #print auc_simu
    #print auc_data
    #err = np.sqrt( np.mean( np.square( blood_conc_data - blood_conc ) ) )
    err = np.sqrt( np.mean( np.square( auc_data - auc_simu ) ) )
    #print err
    return err