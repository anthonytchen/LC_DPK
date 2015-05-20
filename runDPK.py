import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import Skin as skin
import Chemical as chem

n_layer_x_sc = 16; 
n_layer_y_sc = 2;
n_grids_x_ve = 10;
n_grids_x_de = 10;
offset_y_sc = 0;
b_inf_src = 0;

# Nicotine
MW = 162.23;
K_ow = pow(10, 1.17);
pKa = 3.12; # a base

dx_vehicle = 1e-4; # thickness of vehicle in meter
area_vehicle = 3.5*1e-4; # cm2 patch, represented in m2
conc_vehicle = 30.0*1e-6/(area_vehicle*dx_vehicle); # mg in cm2 patch, with patch thickness 0.1cm; in kg/m3
diffu_vehicle = 5e-14; # diffusivity of solute in vehicle
  
 
_chem = chem.Chemical();
_skin = skin.Skin();
  
_chem.Init(MW, K_ow, pKa, 'B'); # the last letter denotes acid (A) or base (B)
_skin.Init(_chem, conc_vehicle, diffu_vehicle, dx_vehicle, area_vehicle, 
           n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de, offset_y_sc, b_inf_src);
  
conc = np.zeros(_skin.m_dim_all);
_skin.getGridsConc(conc);

coord_x = np.zeros(_skin.m_dim_all);
coord_y = np.zeros(_skin.m_dim_all);
_skin.getXCoord(coord_x);
_skin.getYCoord(coord_y);

''' Class definition for structured concentration profiles
'''
class strdConc:
    def __init__(self):
        self.vehicle = []
        self.stracorn2d_coordx = []
        self.stracorn2d_coordy = []
        self.stracorn2d_data = []
        self.stracorn1d = []
        self.viaepd = []
        self.dermis = []
        self.blood = []

stdConc = strdConc();


''' The function to work out concentration at different depth
'''
def calConcProfile(conc_obj, skin_obj, conc_raw, coord_x, coord_y):

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

t_end = 3600*240
t_inv = 3600
for t in range(0, t_end, t_inv):

    print( 'time = %.3e\t' % t )

    _skin.diffuseMoL(t, t+t_inv)
    _skin.getGridsConc(conc);
    calConcProfile(stdConc, _skin, conc, coord_x, coord_y)

    flux1 = _skin.compFlux_2sc();
    flux2 = _skin.compFlux_sc2down();
    flux3 = _skin.compFlux_ve2down();
    flux4 =  _skin.compFlux_de2down();

    print( 'flux_2sc = %.3e, flux2ve = %e, flux_2de = %.3e, flux_2bd = %.3e\n' % (flux1, flux2, flux3, flux4) )

    plt.imshow( stdConc.stracorn2d_data )
    plt.plot( np.hstack( (stdConc.stracorn1d[0,:], stdConc.viaepd[0,:], stdConc.dermis[0,:]) ), 
              np.hstack( (stdConc.stracorn1d[1,:], stdConc.viaepd[1,:], stdConc.dermis[1,:]) ), 'r--o' )
