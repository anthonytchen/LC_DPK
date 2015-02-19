// LCP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "arg.h"
#include "Skin.h"


int main (int argc, char* argv[])
{
  double g, d, s, t, concSource, *K_ow, *MW, *perm_data, *perm_simu, DSource,
    t_simu, t_end, t_inv, offset_y, temp1, temp2;
  int i, n_layer_x, n_layer_y, n_data;
	
  static int nDis = 1;
  static char *fn_data = "sc_perm.dat";
	
  t_end = 60*60*48; t_inv = 60*60*2; // simulation time and interval in seconds
  g=0.075E-6; d=40E-6; s=0.075E-6; t=0.8E-6;
  n_layer_x = 16; //16
  n_layer_y = 1; // 2
  offset_y = 0;

  concSource = .1 * 1e3; // in mol/m^3
  DSource = 9.12e-10; // diffusivity of solute in vehicle
	
  // Provide a command line user interface
  static Config_t params[] = {

    "tinv", "Simulation time interval (s)",
    "-tinv", DOUBLE, (caddr_t)&t_inv,

    "tend", "Simulation end time (s)",
    "-tend", DOUBLE, (caddr_t)&t_end,  
  
    "dis", "Display options; 0 - most parsimonious; 3 - most verbose",
    "-dis", INT, (caddr_t)&nDis,

    "fn_data", "File name to read data from",
    "-fn_data", STRING, (caddr_t)&fn_data,

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

  Skin _skin;
  
  FILE *fp = fopen(fn_data, "r");
  if ( fp == NULL ){
    fclose(fp);
    exit(0);
  }
  fscanf(fp, "%d", &n_data);

  K_ow = new double [n_data];
  MW = new double [n_data];
  perm_data = new double[n_data];
  perm_simu = new double[n_data];

  for ( i = 0; i < n_data; i++ ) {
    // read data
    fscanf(fp, "%lf", K_ow+i);
    K_ow[i] = pow(10, K_ow[i]);
    fscanf(fp, "%lf", MW+i);
    fscanf(fp, "%lf", perm_data+i);

    printf("K_ow=%.2lf, MW=%.2lf\n", log10(K_ow[i]), MW[i]);
    
    _skin.Init( g, d, s, t, K_ow[i], MW[i], concSource, DSource, n_layer_x, n_layer_y, t_inv, offset_y );
    _skin.createGrids();

    // simulate for each of the data points
    perm_simu[i] = -100;
    for ( t_simu=.0; t_simu<t_end; t_simu+=t_inv ){

      _skin.diffuseMoL_cv(t_simu, t_simu+t_inv);

      temp1 = _skin.compFlux_2sc()/concSource; temp1 = log10(temp1);
      temp2 = _skin.compFlux_sc2ve()/concSource; temp2 = log10(temp2);
      printf("\t %.2lf %.2lf %.2lf\n", t_simu+t_inv, temp1, temp2);

      if ( fabs(temp1-temp2) < 1e-3 )
	break;
    }

    perm_simu[i] = temp1;
  }

  for ( i = 0; i < n_data; i++ )
    printf("%.2f  %.2f\n", perm_data[i], perm_simu[i]);    

  fclose(fp); 
  
  delete [] K_ow;
  delete [] MW;
  delete [] perm_data;
  delete [] perm_simu;
  _skin.Release();

  return 0;
}
