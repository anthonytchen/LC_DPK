// LCP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "arg.h"
#include "Chemical.h"
#include "Skin.h"


int main (int argc, char* argv[])
{
  double MW, K_ow, pKa,
    conc_vehicle, diffu_vehicle, partition_vehicle, partition_dermis2blood, k_clear_blood, area_vehicle,
    t_simu, t_end, t_inv, offset_y_sc, dx_vehicle, x_len_viaepd, x_len_dermis;
  bool b_1st_save = true;
  int b_inf_src = 0;
  int i, n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de;
  char fn_coord_x[20], fn_coord_y[20];
	
  static int nDis = 1;
  static char *fn_conc = "conc";
  static char *pre_coord = "coord";
	
  t_end = 900; t_inv = 10; // simulation time and interval ()in seconds
  n_layer_x_sc = 16+5; //16
  n_layer_y_sc = 2; // 2
  n_grids_x_ve = 10;
  n_grids_x_de = 10;
  offset_y_sc = 0;
  x_len_viaepd = 100e-6; // depth of viable epidermis
  x_len_dermis = 1200e-6; // depth of dermis
  
  // 4-Cyanophenol
  /*
  MW = 119.12;  // Da, i.e. g/mol
  K_ow = pow(10,1.6); // partition coefficient between octanol and water
  pKa = 7.8; 
  diffu_vehicle = 9.12e-10; // diffusivity of solute in vehicle
  */

  // DEET
  /*
  MW = 191.3;
  K_ow = pow(10, 2.02);
  pKa = 3.12; // not used

  dx_vehicle = ;
  area_vehicle = 4*6*1e-4; // cm2 application area, represented in m2
  conc_vehicle = 15.0*1e-6/(area_vehicle*dx_vehicle); // mg in cm2 patch, with patch thickness 0.1cm; in kg/m3
  diffu_vehicle = 1; // diffusivity of solute in vehicle
  */
  // Nicotine
  
  MW = 162.23;
  K_ow = pow(10, 1.17);
  pKa = 3.12; // a base

  dx_vehicle = 1e-4; // thickness of vehicle in meter
  area_vehicle = 3.5*1e-4; // cm2 patch, represented in m2
  conc_vehicle = 15.0*1e-6/(area_vehicle*dx_vehicle); // mg in cm2 patch, with patch thickness 0.1cm; in kg/m3
  diffu_vehicle = 1e-13; // diffusivity of solute in vehicle
  partition_vehicle = 1;
  partition_dermis2blood = 1.0/pow(10,0.04);
  k_clear_blood = 23.3e-6; // 1400 ml/min = 23.3e-6 m3/s, 1250 = 20.8e-6, 1540 = 25.7e-6
  
	
  // Provide a command line user interface
  static Config_t params[] = {

    "tinv", "Simulation time interval (s)",
    "-tinv", DOUBLE, (caddr_t)&t_inv,

    "tend", "Simulation end time (s)",
    "-tend", DOUBLE, (caddr_t)&t_end,  

    "offset_y", "Offset on the left boundary",
    "-ofy", DOUBLE, (caddr_t)&offset_y_sc,
    
    "n_layer_y", "Number of y (lateral) layers",
    "-ny", INT, (caddr_t)&n_layer_y_sc,

    "fn_conc", "File name to store concentration",
    "-fn_conc", STRING, (caddr_t)&fn_conc,

    "inf_src", "Whether the vehicle is an infinite source",
    "-inf_src", INT, (caddr_t)&b_inf_src,
	  
    "pre_coord", "Prefix (max 10 chars) to store coordinates to <fn_coord>_x.txt and <fn_coord>_y.txt",
    "-pre_coord", STRING, (caddr_t)&pre_coord,
	  
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
  
  clock_t start, end;
  double cpu_time_used;
  
  Chemical _chem;
  Skin _skin;
  
  //  _chem.Init(MW, K_ow, pKa, 'A'); // the last letter denotes acid (A) or base (B)
  _chem.Init(MW, K_ow, pKa, 0.31, 0.95, 'B'); // the last letter denotes acid (A) or base (B)

  _skin.Init(_chem, conc_vehicle, diffu_vehicle, partition_vehicle, partition_dermis2blood, k_clear_blood,
	     dx_vehicle, area_vehicle, x_len_viaepd, x_len_dermis,
	     n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de, offset_y_sc, b_inf_src);
  
  strcpy(fn_coord_x, pre_coord); strcat(fn_coord_x, "_x");
  strcpy(fn_coord_y, pre_coord); strcat(fn_coord_y, "_y");
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

    flux1 = _skin.compFlux_2sc();
    flux2 = _skin.compFlux_sc2down();
    flux3 = _skin.compFlux_ve2down();
    flux4 =  _skin.compFlux_de2down();
    printf("Time %e flux = %e %e %e %e\n", t_simu+t_inv, flux1, flux2, flux3, flux4);

    // _skin.resetVehicle(conc_vehicle, partition_vehicle, diffu_vehicle);
  }

  _skin.Release();
  return 0;
}

