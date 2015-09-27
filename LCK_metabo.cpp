// LCP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "arg.h"
#include "Chemical.h"
#include "Skin.h"


int main (int argc, char* argv[])
{
  double MW, K_ow, pKa, f_non, f_u,
    conc_vehicle[2], diffu_vehicle[2], partition_vehicle[2], partition_dermis2blood[2], k_clear_blood[2], area_vehicle,
    t_simu, t_end, t_inv, offset_y_sc, dx_vehicle, x_len_viaepd, x_len_dermis;
  bool b_1st_save = true;
  int b_inf_src = 1;
  int i, n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de;
  char fn_coord_x[20], fn_coord_y[20];
	
  static int nDis = 1;
  static char *fn_conc = "conc";
  static char *pre_coord = "coord";
	
  // rat abdominal skin
  t_end = 900; t_inv = 10; // simulation time and interval ()in seconds
  n_layer_x_sc = 1; //16
  n_layer_y_sc = 1; // 2
  n_grids_x_ve = 10;
  n_grids_x_de = 10;
  offset_y_sc = 0;
  x_len_viaepd = 7e-6; // depth of viable epidermis
  x_len_dermis = 834e-6; // depth of dermis
  

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
  
  Chemical _chem[2];
  Skin _skin;

  area_vehicle = 0.95*1e-4; // cm2 effective diffusion are, represented in m2
  dx_vehicle = 1e-4; // thickness of vehicle in meter, assumed
  diffu_vehicle[0] = diffu_vehicle[1] = 1e-9; // diffusivity of solute in vehicle, assumed
  partition_dermis2blood[0] = partition_dermis2blood[1] = 1.0; // not used
  k_clear_blood[0] = k_clear_blood[1] = 23.3e-6; // m3/s, not used  

  
  // 1. Parent chemical: ethyl nicotinate (EN)
  MW = 151.16;  K_ow = pow(10, 1.32);  pKa = 3.24; // a base
  f_non = 0.99; f_u = 0.31; // f_u need verification
  _chem[0].Init(MW, K_ow, pKa, f_non, f_u, 'B'); // the last letter denotes acid (A) or base (B)
  conc_vehicle[0] = 145; // mol/m3
  partition_vehicle[0] = 1; // assuming aqueous solution

  // 2. Metabolite: nicotinic acid (NA)
  MW = 123.11;  K_ow = pow(10, 0.36);  pKa = 4.75; // an acid
  f_non = 0.02; f_u = 0.31; // f_u need verification
  _chem[1].Init(MW, K_ow, pKa, f_non, f_u, 'A'); // the last letter denotes acid (A) or base (B)
  conc_vehicle[1] = 0;
  partition_vehicle[0] = 1; // assuming aqueous solution

  bool b_has_blood = false;
  _skin.Init(_chem, 2, b_has_blood, conc_vehicle, diffu_vehicle, partition_vehicle, partition_dermis2blood, k_clear_blood,
	     dx_vehicle, area_vehicle, x_len_viaepd, x_len_dermis,
	     n_layer_x_sc, n_layer_y_sc, n_grids_x_ve, n_grids_x_de, offset_y_sc, b_inf_src);
  _skin.InitReaction(0, 1, 5.7/3600, 0.837); // 5.7 mol m-3 h-1 converted to mol m-3 s-1
  
  strcpy(fn_coord_x, pre_coord); strcat(fn_coord_x, "_x");
  strcpy(fn_coord_y, pre_coord); strcat(fn_coord_y, "_y");
  _skin.saveCoord( fn_coord_x, fn_coord_y );
	
  if ( nDis > 1 )
    _skin.displayGrids();	


  double flux1[2], flux2[2], flux3[2], flux4[2];
  
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

    _skin.compFlux_2sc(flux1);
    _skin.compFlux_sc2down(flux2);
    _skin.compFlux_ve2down(flux3);
     _skin.compFlux_de2down(flux4);
    printf("  Time %e flux = %e %e %e %e\n", t_simu+t_inv, flux1[0], flux2[0], flux3[0], flux4[0]);
    printf("  Time %e flux = %e %e %e %e\n", t_simu+t_inv, flux1[1], flux2[1], flux3[1], flux4[1]);

    // _skin.resetVehicle(conc_vehicle, partition_vehicle, diffu_vehicle);
  }

  _skin.Release();
  return 0;
}

