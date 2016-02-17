// The console application for predicting stratum corneum permeability
//    of infinite dose (default) and sink conditions
//

#include "stdafx.h"
#include "arg.h"
#include "Chemical.h"
#include "Skin_VS.h"


int main (int argc, char* argv[])
{
  double MW, log_K_ow, log_K_vh, K_vh, pKa,
    conc_vehicle, diffu_vehicle, area_vehicle,
    t_simu, t_end, t_inv, dx_vehicle;
  bool b_1st_save = true;
  int b_inf_src = 1;
  int i, n_layer_x_sc, n_layer_y_sc;
  char fn_coord_x[1024], fn_coord_y[1024], fn_flux[1024];
	
  static int nDis = 1;
  static char *fn = "chemical_name";
  char fn_conc[1024] = "conc";
	
  t_end = 900; t_inv = 10; // simulation time and interval ()in seconds
  n_layer_x_sc = 12;
  n_layer_y_sc = 1;
  
  // Nicotine
  
  MW = 162.23;
  log_K_ow = 1.17;
  pKa = 3.12; // a base

  dx_vehicle = 1e-4; // thickness of vehicle in meter
  area_vehicle = 3.5*1e-4; // cm2 patch, represented in m2
  conc_vehicle = 1; // in kg/m3, same as mg/ml
  diffu_vehicle = -1; // diffusivity of solute in vehicle; negative value means using diffusivity in water
  log_K_vh = 0; // essentially water vehicle
	
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

    "C_vh", "concentration (kg/m3 = mg/ml) in vehicle",
    "-Cvh", DOUBLE, (caddr_t)&conc_vehicle,
    
    "n_layer_x_sc", "Number of (verticle) cell layers in stratum corneum",
    "-nx", INT, (caddr_t)&n_layer_x_sc,

    "fn", "File name to store concentration, coordinates, and flux history",
    "-fn", STRING, (caddr_t)&fn,

    "inf_src", "Whether the vehicle is an infinite source",
    "-inf_src", INT, (caddr_t)&b_inf_src,
	  
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
  sprintf(fn_flux, "%s.flux", fn);

  clock_t start, end;
  double cpu_time_used;
  
  Chemical _chem;
  _chem.Init(MW, pow(10, log_K_ow), pKa, 1, 1, 'B'); // the last 3 inputs are not used when calculating permeability through SC

  Skin_VS _skin;

  K_vh = pow(10, log_K_vh);
  _skin.Init(&_chem, 1, &conc_vehicle, &K_vh, &diffu_vehicle,
	     dx_vehicle, area_vehicle, n_layer_x_sc, n_layer_y_sc, 40.03751e-6, b_inf_src);
  _skin.saveCoord( fn_coord_x, fn_coord_y );
  if ( nDis > 1 )
    _skin.displayGrids();	

  double flux1, flux2;

  FILE *file = fopen(fn_flux, "w");

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

    _skin.compFlux_2sc(&flux1);
    _skin.compFlux_sc2down(&flux2);
    fprintf(file, "Time %e flux = %e %e\n", t_simu+t_inv, flux1, flux2);
    if ( fabs(flux1-flux2)/flux1 < 1e-4 )
      break;
  }
  printf("Steady-state flux = %e, permeability = %e\n", flux1, flux1/conc_vehicle);

  fclose(file);
  _skin.Release();
  return 0;
}

