/*! 
  The console application for predicting stratum corneum permeability
    of infinite dose (default) and sink conditions
*/

#include "stdafx.h"
#include "arg.h"
#include "Config.h"
#include "Chemical.h"
#include "Skin_Setup.h"


int main (int argc, char* argv[])
{
  double t_simu, t_end, t_inv;
  bool b_1st_save = true;
  int i;
  char fn_coord_x[1024], fn_coord_y[1024], fn_amount[1024], fn_flux[1024];
	
  static int nDis = 1;
  static char *dir = "Dir";
  static char *cfn = "Nicotine.cfg";
  char fn_conc[1024] = "conc";
	
  t_end = 900; t_inv = 10; // simulation time and interval in seconds

  // Provide a command line user interface
  static Config_t params[] = {

    "tinv", "Simulation time interval (s)",
    "-tinv", DOUBLE, (caddr_t)&t_inv,

    "tend", "Simulation end time (s)",
    "-tend", DOUBLE, (caddr_t)&t_end,  

    "cfn", "File name of the configuration file",
    "-cfn", STRING, (caddr_t)&cfn,
    
    "dir", "Directory to store simulation data",
    "-dir", STRING, (caddr_t)&dir,

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

  sprintf(fn_conc, "mkdir -p ./%s", dir);
  system(fn_conc);
  sprintf(fn_conc, "./%s/conc", dir);
  sprintf(fn_coord_x, "./%s/coord_x", dir);
  sprintf(fn_coord_y, "./%s/coord_y", dir);
  sprintf(fn_amount, "./%s/amount", dir);
  sprintf(fn_flux, "./%s/flux", dir);

  clock_t start, end;
  double cpu_time_used;
  
  Config _conf;
  _conf.ReadConfigFile(cfn);

  Chemical _chem;
  _chem.InitConfig(_conf);


  Skin_Setup _skin;
  _skin.InitConfig(&_chem, _conf);

  _skin.saveCoord( fn_coord_x, fn_coord_y );
  if ( nDis > 1 )
    _skin.displayGrids();	

  double flux;

  _skin.saveGrids(b_1st_save, fn_conc);
  _skin.compCompartAmount();
  _skin.saveAmount(b_1st_save, fn_amount);
  b_1st_save = !b_1st_save;
      
  FILE *file = fopen(fn_flux, "w");

  for ( t_simu=.0; t_simu<t_end; t_simu+=t_inv ){
    start = clock();
   
    _skin.diffuseMoL(t_simu, t_simu+t_inv);	
    _skin.saveGrids(b_1st_save, fn_conc);
    _skin.compCompartAmount();
    _skin.saveAmount(b_1st_save, fn_amount);
	
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		
    if ( nDis > 0 ) {
      printf("Simulation time is %e, cpu time = %e s \n", t_simu+t_inv, cpu_time_used);
      fflush(stdout);
    }

    // when flux reaches steady state, break
    // only applicable for permeability calculation
    if ( _skin.saveFlux(b_1st_save, fn_flux, &flux) && _conf.m_bInfSrc ) {
      printf("Steady-state flux = %e, permeability = %e\n", flux, flux/_conf.m_conc_vehicle);
      break;
    }
  }

  fclose(file);
  _skin.Release();
  return 0;
}

