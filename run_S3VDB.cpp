// The console application for predicting dermal and systemic kinetics
//

#include "stdafx.h"
#include "arg.h"
#include "Chemical.h"
#include "Skin_S3VDB.h"


int main (int argc, char* argv[])
{
  double MW, log_K_ow, log_K_vh, K_vh, pKa,
    x_len_sb_sur, y_len_sb_har, x_len_viaepd, x_len_dermis,
    conc_vehicle, diffu_vehicle, partition_dermis2blood, k_clear_blood, area_vehicle,
    t_simu, t_end, t_inv, t_remove, offset_y_sc;
  bool b_1st_save = true;
  int b_inf_src = 0;
  int i, n_layer_x_sc, n_layer_y_sc, n_grids_x_sb_sur, n_grids_y_sb_har, n_grids_x_ve, n_grids_x_de;
  char fn_coord_x[1024], fn_coord_y[1024];
	
  static int nDis = 1;
  static char *fn = "chemical_name";
  char fn_conc[1024] = "conc";

	
  t_end = 900; t_inv = 10; t_remove = t_end/3600; // simulation time and interval in seconds
  n_layer_x_sc = 23; //16
  n_layer_y_sc = 1; // 2
  n_grids_x_ve = 10;
  n_grids_x_de = 10;
  offset_y_sc = 40.03751e-6;
  x_len_viaepd = 100e-6; // depth of viable epidermis
  x_len_dermis = 1200e-6; // depth of dermis
  n_grids_x_sb_sur = 10;  
  n_grids_y_sb_har = 5;

  // Nicotine
  
  MW = 162.23;
  log_K_ow = 1.17;
  pKa = 3.12; // a base

  /* vanakoski  paper */
  x_len_sb_sur = 20e-6; // thickness of surface sebum
  y_len_sb_har = 20e-6; // width of hair follicle sebum
  //area_vehicle = 22.5*1e-4; // cm2 patch, represented in m2
  //conc_vehicle = 39.37*1e-6/(area_vehicle*dx_vehicle); // mg in cm2 patch, with patch thickness 0.1cm; in kg/m3
  diffu_vehicle = -1; // using diffusivity in water
  log_K_vh = log10(0.7);
  partition_dermis2blood = 1.0/pow(10,0.04);
  k_clear_blood = 23.3e-6; // 1400 ml/min = 23.3e-6 m3/s, 1250 = 20.8e-6, 1540 = 25.7e-6

	
  // Provide a command line user interface
  static Config_t params[] = {

    "tinv", "Simulation time interval (s)",
    "-tinv", DOUBLE, (caddr_t)&t_inv,

    "tend", "Simulation end time (s)",
    "-tend", DOUBLE, (caddr_t)&t_end,  

    "MW", "molecular weight",
    "-MW", DOUBLE, (caddr_t)&MW,

    "log_Kow", "log10 of partition octanol:water",
    "-Kow", DOUBLE, (caddr_t)&log_K_ow,

    "log_K_vh", "log10 of partition vehicle:water",
    "-Kvh", DOUBLE, (caddr_t)&log_K_vh,

    "n_layer_x_sc", "Number of layers (x-axis, verticle) of cells in stratum corneum",
    "-nxsc", INT, (caddr_t)&n_layer_x_sc,
    
    "n_layer_y_sc", "Number of layers (y-axis, lateral) of cells in stratum corneum",
    "-nysc", INT, (caddr_t)&n_layer_y_sc,

    "x_len_sb_sur", "Depth of surface sebum (m)",
    "-xsbs", DOUBLE, (caddr_t)&x_len_sb_sur,

    "y_len_sb_har", "Width of hair sebum (m)",
    "-ysbh", DOUBLE, (caddr_t)&y_len_sb_har,

    "x_len_viaepd", "Depth of viable epidermis (m)",
    "-xve", DOUBLE, (caddr_t)&x_len_viaepd,

    "x_len_dermis", "Depth of dermis (m)",
    "-xde", DOUBLE, (caddr_t)&x_len_dermis,

    "fn", "File name prefix to store concentration, coordinates, etc",
    "-fn", STRING, (caddr_t)&fn,
   
    "dis", "Display options; 0 - most parsimonious; 3 - most verbose",
    "-dis", INT, (caddr_t)&nDis,
	  
    0, 0, 0, NOTYPE, 0

  };
  if ( argc < 2) {
    pusage ( argv[0], params );
    exit (0);
  }
  if ( ppconf (argc, argv, params, true) ) {
    pusage (argv[0], params);
    exit (1);
  }
  
  sprintf(fn_conc, "%s.conc", fn);
  sprintf(fn_coord_x, "%s.coord_x", fn);
  sprintf(fn_coord_y, "%s.coord_y", fn);

  clock_t start, end;
  double cpu_time_used;
  
  Chemical _chem;
  _chem.Init(MW, pow(10, log_K_ow), pKa, 0.31, 0.95, 'B'); // the last letter denotes acid (A) or base (B)

  Skin_S3VDB _skin;
  K_vh = pow(10, log_K_vh);
  _skin.Init(&_chem, 1, &conc_vehicle, &K_vh, &diffu_vehicle,
	     x_len_sb_sur, n_grids_x_sb_sur, y_len_sb_har, n_grids_y_sb_har,
	     n_layer_x_sc, n_layer_y_sc, offset_y_sc,
	     x_len_viaepd, n_grids_x_ve, x_len_dermis, n_grids_x_de,
	     &partition_dermis2blood, &k_clear_blood);
  _skin.saveCoord( fn_coord_x, fn_coord_y );
	
  if ( nDis > 1 )
    _skin.displayGrids();	


  double flux1, flux2, flux3, flux4;
  
  for ( t_simu=.0; t_simu<t_end; t_simu+=t_inv ){
    start = clock();
   
    _skin.diffuseMoL(t_simu, t_simu+t_inv);	
    _skin.saveGrids(b_1st_save, fn_conc);

    if ( b_1st_save )
      b_1st_save = !b_1st_save;
		
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		
    if ( nDis > 0 ) {
      printf("Simulation time is %e, cpu time = %e s \n", t_simu+t_inv, cpu_time_used);
      fflush(stdout);
    }

    if ( t_simu+t_inv-3600*t_remove > -1 )
      _skin.removeVehicle();
  }

  _skin.Release();
  return 0;
}

